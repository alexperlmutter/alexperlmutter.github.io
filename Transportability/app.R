library(shiny)
library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)

expit <- function(x) 1 / (1 + exp(-x))

# Helper function: Standardized Mean Difference (SMD)
calc_smd <- function(val, group, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(val))

  w1 <- weights[group == 1]
  w0 <- weights[group == 0]
  v1 <- val[group == 1]
  v0 <- val[group == 0]

  if (sum(w1) == 0 || sum(w0) == 0) return(NA)

  m1 <- sum(v1 * w1) / sum(w1)
  m0 <- sum(v0 * w0) / sum(w0)

  var1 <- sum(w1 * (v1 - m1)^2) / sum(w1)
  var0 <- sum(w0 * (v0 - m0)^2) / sum(w0)

  denom <- sqrt((var1 + var0) / 2)
  if (is.na(denom) || denom == 0) return(0)

  return((m1 - m0) / denom)
}

# Monte Carlo single iteration worker
run_single_sim <- function(input) {
  N <- input$n_pop
  p_unlinked <- input$pct_unlinked / 100

  # 1. Assign Linkage Status: S=1 (Unlinked) vs S=0 (Linked)
  S <- rbinom(N, 1, p_unlinked)

  # 2. Simulate Claims Covariates (L) conditional on Linkage Status (S)
  white <- ifelse(S == 0, rbinom(N, 1, input$p_white_link/100), rbinom(N, 1, input$p_white_unlink/100))
  diabetes <- ifelse(S == 0, rbinom(N, 1, input$p_diab_link/100), rbinom(N, 1, input$p_diab_unlink/100))
  cv_hosp <- ifelse(S == 0, rbinom(N, 1, input$p_cv_link/100), rbinom(N, 1, input$p_cv_unlink/100))
  endo_visit <- ifelse(S == 0, rbinom(N, 1, input$p_endo_link/100), rbinom(N, 1, input$p_endo_unlink/100))
  fall_hist <- ifelse(S == 0, rbinom(N, 1, input$p_fall_link/100), rbinom(N, 1, input$p_fall_unlink/100))

  # 3. Treatment Assignment (A) driven by Claims L
  lin_a <- -0.5 + 0.5 * white + 0.9 * diabetes + 0.7 * cv_hosp + 0.8 * endo_visit + 0.4 * fall_hist
  A <- rbinom(N, 1, expit(lin_a))

  # 4. Selective Test Ordering (O) driven by L and A
  beta_diab <- log(input$or_ord_diab)
  beta_endo <- log(input$or_ord_endo)
  beta_tx   <- log(input$or_ord_tx)
  target_p  <- input$p_hba1c_order / 100
  base_offset <- log(target_p / (1 - target_p)) - 1.0

  lin_o_hba1c <- base_offset + beta_diab * diabetes + beta_endo * endo_visit + beta_tx * A
  O_hba1c <- ifelse(S == 0, rbinom(N, 1, expit(lin_o_hba1c)), 0)

  # 5. Generate Biometrics U driven by Claims L and Ordering O
  order_effect <- input$order_on_u_effect
  hba1c <- 5.5 + 1.8 * diabetes + 0.8 * endo_visit + 0.3 * cv_hosp + order_effect * O_hba1c + rnorm(N, 0, 0.6)
  ldl <- 100 + 25 * cv_hosp + 15 * diabetes + 10 * endo_visit + rnorm(N, 0, 18)

  df <- data.frame(id = 1:N, S = S, A = A, white = white, diabetes = diabetes,
                   cv_hosp = cv_hosp, endo_visit = endo_visit, fall_hist = fall_hist,
                   hba1c = hba1c, ldl = ldl, O_hba1c = O_hba1c)

  # Filter complete cases in linked data (S=0 & O_hba1c=1)
  df_cc <- df %>% filter(S == 0 & O_hba1c == 1)

  # --------------------------------------------------------------------------
  # WEIGHT ESTIMATION
  # --------------------------------------------------------------------------
  # 1. IPTW Model (Targets Tested Subpopulation: S=0, O=1)
  ps_mod <- glm(A ~ white + diabetes + cv_hosp + endo_visit + fall_hist, data = df_cc, family = binomial())
  df_cc$p_a <- pmin(pmax(predict(ps_mod, type = "response"), 0.01), 0.99)
  df_cc$W_IPTW <- ifelse(df_cc$A == 1, 1 / df_cc$p_a, 1 / (1 - df_cc$p_a))

  # 2. SMR Weight Model (Transports to Unlinked: S=1)
  smr_mod <- glm(S ~ white + diabetes + cv_hosp + endo_visit + fall_hist, data = df, family = binomial())
  df$p_s <- pmin(pmax(predict(smr_mod, type = "response"), 0.01), 0.99)
  df$W_SMR <- ifelse(df$S == 0, df$p_s / (1 - df$p_s), 1.0)
  df_cc <- df_cc %>% left_join(df %>% select(id, W_SMR), by = "id")

  # 3. IPSW Weight Model (Transports to Full Linked: S=0)
  df_linked <- df %>% filter(S == 0)
  ord_mod <- glm(O_hba1c ~ diabetes + endo_visit + A, data = df_linked, family = binomial())
  df_linked$p_o <- pmin(pmax(predict(ord_mod, type = "response"), 0.01), 0.99)
  df_linked$W_IPSW <- 1 / df_linked$p_o
  df_cc <- df_cc %>% left_join(df_linked %>% select(id, W_IPSW), by = "id")

  # Composite Weight Definitions
  df_cc <- df_cc %>%
    mutate(
      W_Unadjusted = 1.0,
      W_IPTW_Only  = W_IPTW,
      W_IPTW_SMR   = W_IPTW * W_SMR,
      W_IPTW_IPSW  = W_IPTW * W_IPSW,
      W_Fully_Adj  = W_IPTW * W_SMR * W_IPSW
    )

  scenarios <- c(
    "1. Unweighted (Observed Sample)",
    "2. IPTW Only [Target: Tested S=0, O=1]",
    "3. IPTW + SMR [Target: Unlinked Tested S=1, O=1]",
    "4. IPTW + IPSW [Target: Full Linked S=0]",
    "5. IPTW + SMR + IPSW [Target: Total Population N]"
  )
  w_vars <- c("W_Unadjusted", "W_IPTW_Only", "W_IPTW_SMR", "W_IPTW_IPSW", "W_Fully_Adj")

  smd_res <- lapply(1:5, function(i) {
    w_vec <- df_cc[[w_vars[i]]]
    data.frame(
      Scenario = scenarios[i],
      HbA1c_SMD = calc_smd(df_cc$hba1c, df_cc$A, w_vec),
      LDL_SMD   = calc_smd(df_cc$ldl, df_cc$A, w_vec)
    )
  })

  bind_rows(smd_res)
}

library(shiny)

# ==============================================================================
# APP OVERVIEW & THEORETICAL SUMMARY UI PANEL
# ==============================================================================
library(shiny)

# ==============================================================================
# APP OVERVIEW & THEORETICAL SUMMARY UI PANEL (ROBUST MATHJAX)
# ==============================================================================
overview_panel <- wellPanel(
  style = "background-color: #f8f9fa; border-left: 5px solid #0056b3; margin-bottom: 20px;",
  withMathJax(
    tags$div(
      tags$h2("App Purpose & Theoretical Overview"),
      tags$p(
        "The primary purpose of this Shiny application is to demonstrate how proxy-adjustment strategies can control unmeasured confounding in a target population (unlinked claims, \\(S=1\\)) using clinical proxy variables (\\(L\\)) measured in an analytical study population (linked claims–EHR, \\(S=0\\))."
      ),
      tags$p(
        "When an unmeasured confounder \\(U_k\\) (e.g., baseline HbA1c categorized as \\(U_k=1\\) [\\(\\ge 6.5\\%\\)] vs. \\(U_k=0\\) [< 6.5%]) is missing in claims data (\\(S=1\\)), observed claims covariates \\(L\\) can serve as proxies for \\(U_k\\). Standard Inverse Probability of Treatment Weighting (IPTW) on \\(L\\) balances treatment groups (\\(A\\)) with respect to \\(L\\) and, by proxy, \\(U_k\\) within \\(S=0\\). However, because the goal is to address confounding due to \\(U_k\\) in \\(S=1\\) rather than \\(S=0\\), proxy adjustment in the study population alone is insufficient if the joint distribution of \\(L\\) and \\(A\\) differs substantially between \\(S=0\\) and \\(S=1\\)."
      ),
      tags$p(
        "Under proxy sufficiency, proxy adjustment remains consistent across populations if the underlying confounding and proxy mechanisms operate similarly (e.g., \\(L\\) increases the probability of \\(U_k\\) and \\(A\\) in both \\(S=0\\) and \\(S=1\\)). In effect, population linkage may act as a weak effect modifier of the treatment–confounder relationship:"
      ),
      tags$p(
        "\\[ P(A=1 \\mid U_k=1, L, S=0) - P(A=1 \\mid U_k=0, L, S=0) \\neq P(A=1 \\mid U_k=1, L, S=1) - P(A=1 \\mid U_k=0, L, S=1) \\]"
      ),
      tags$p(
        "Additionally, when \\(U_k\\) is a selectively ordered biomarker (available only when test order \\(O_k=1\\)), IPTW on \\(L\\) balances \\(U_k\\) across treatment arms among tested patients (\\(S=0, O_k=1\\)). However, this balance does not automatically extend to untested or unlinked populations if \\(S=0\\) and \\(O_k=1\\) systematically differ in \\(L\\)."
      ),
      tags$p(
        "To transport proxy-adjusted inferences to the target population (\\(S=1\\)), a multi-stage weighting pipeline is required:"
      ),
      tags$ul(
        tags$li(tags$strong("IPTW (\\(W_{\\text{IPTW}}\\)): "), "Adjusts for confounding (differences in \\(U\\) and \\(L\\) between \\(A\\)) based only on measured proxies \\(L\\)."),
        tags$li(tags$strong("IPSW (\\(W_{\\text{IPSW}}\\)): "), "Adjusts for selective test ordering \\(O_k\\) to transport inferences from tested patients \\(O_k=1\\) to the full linked population \\(S=0\\)."),
        tags$li(tags$strong("SMR (\\(W_{\\text{SMR}}\\)): "), "Adjusts for linkage selection \\(S\\) to transport inferences from the linked population \\(S=0\\) to the unlinked target population \\(S=1\\).")
      ),
      tags$hr(),
      tags$h3("Simulation Mechanics & Key Insights"),
      tags$p(
        "This application uses Monte Carlo simulation to evaluate how combining these weighting strategies affects proxy performance and covariate balance across different target populations:"
      ),
      tags$ul(
        tags$li(
          tags$strong("Tested Study Sample (\\(S=0, O_k=1\\)): "),
          "IPTW alone achieves optimal proxy adjustment for \\(U_k\\) in the observed complete-case sample because test ordering \\(O_k\\) is held constant."
        ),
        tags$li(
          tags$strong("Target Population Transportability (\\(S=1\\) or Total \\(N\\)): "),
          "When transporting inferences to populations that differ from \\(S=0, O_k=1\\), joint weighting schemes (e.g., IPTW + IPSW + SMR) account for selection differences in \\(L\\), \\(S\\), and \\(O_k\\). Proxy adjustment is less effective."
        ),
        tags$li(
          tags$strong("Interactive Assumption Testing: "),
          "Users can dynamically manipulate data-generating parameters—including the unlinked population proportion, proxy strength (\\(L \\rightarrow U_k\\)), linkage selection (\\(L \\rightarrow S\\)), test ordering drivers (\\(L, A \\rightarrow O_k\\)), and causal shifts (\\(O_k \\rightarrow U_k\\))—to evaluate estimator performance under varying observational conditions."
        )
      )
    )
  )
)

# ==============================================================================
# UI DEFINITION
# ==============================================================================
ui <- fluidPage(
  titlePanel("Proxy Assessment & Target Population Transportability"),

  # Insert the overview panel here at the top of the page
  overview_panel,

  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Monte Carlo Settings"),
      numericInput("n_sims", "Number of Simulations:", value = 1000, min = 10, max = 1000, step = 50),
      numericInput("n_pop", "Sample Size (N per sim):", value = 5000, min = 1000, max = 200000, step = 2000),
      sliderInput("pct_unlinked", "Unlinked Population (%):", min = 10, max = 100, value = 90, step = 5),

      hr(),
      h4("Claims Covariates L (Linked vs Unlinked %)"),
      fluidRow(
        column(6, numericInput("p_white_link", "White (L):", 52)),
        column(6, numericInput("p_white_unlink", "White (U):", 30))
      ),
      fluidRow(
        column(6, numericInput("p_diab_link", "Diabetes (L):", 35)),
        column(6, numericInput("p_diab_unlink", "Diabetes (U):", 15))
      ),
      fluidRow(
        column(6, numericInput("p_cv_link", "CV Hosp (L):", 25)),
        column(6, numericInput("p_cv_unlink", "CV Hosp (U):", 10))
      ),
      fluidRow(
        column(6, numericInput("p_endo_link", "Endo Visit (L):", 30)),
        column(6, numericInput("p_endo_unlink", "Endo Visit (U):", 8))
      ),
      fluidRow(
        column(6, numericInput("p_fall_link", "Fall Hist (L):", 20)),
        column(6, numericInput("p_fall_unlink", "Fall Hist (U):", 6))
      ),

      hr(),
      h4("Biometric Ordering Dynamics (L + A -> O)"),
      sliderInput("p_hba1c_order", "HbA1c Ordering Rate (%):", min = 20, max = 90, value = 60),
      sliderInput("or_ord_diab", "OR: Diabetes on Ordering:", min = 1.0, max = 5.0, value = 3.0, step = 0.25),
      sliderInput("or_ord_endo", "OR: Endo Visit on Ordering:", min = 1.0, max = 5.0, value = 2.5, step = 0.25),
      sliderInput("or_ord_tx", "OR: Treatment A on Ordering:", min = 1.0, max = 5.0, value = 2.0, step = 0.25),

      hr(),
      h4("Ordering Effect on Biomarker (O -> U)"),
      sliderInput("order_on_u_effect", "Shift in HbA1c Due to Test Ordering (O -> U):",
                  min = 0.0, max = 2.0, value = 0.6, step = 0.1),

      actionButton("run_mc", "Run Monte Carlo Simulation", class = "btn-primary", style = "width: 100%; margin-top: 15px;")
    ),

    mainPanel(
      width = 8,
      wellPanel(
        h4("Target Population Legend"),
        p(strong("IPTW Only:"), " Adjusts for treatment confounding solely within the ", em("tested subpopulation (S=0, O=1)"), ". Because ordering is constant in this sample, IPTW yields the lowest sample SMD."),
        p(strong("IPTW + IPSW:"), " Transports inferences from tested patients back to the ", em("Full Linked Population (S=0)"), " by adjusting for selective test ordering."),
        p(strong("IPTW + SMR + IPSW:"), " Transports inferences to the ", em("Total Combined Target Population (N)"), " by adjusting for both selective ordering and unlinked data transportability.")
      ),

      tabsetPanel(
        tabPanel("SMD Imbalance Rates (|SMD| > 0.1)",
                 p("Proportion of simulation runs exceeding absolute SMD > 0.1 across weighting strategies."),
                 plotOutput("imbalance_plot", height = "480px"),
                 hr(),
                 tableOutput("imbalance_table")),

        tabPanel("SMD Distribution across Runs",
                 p("Distribution of SMDs across Monte Carlo iterations."),
                 plotOutput("smd_box_plot", height = "480px"))
      )
    )
  )
)

# ==============================================================================
# SERVER LOGIC
# ==============================================================================
server <- function(input, output, session) {

  mc_results <- reactive({
    input$run_mc

    n_sims <- isolate(input$n_sims)

    withProgress(message = 'Running Monte Carlo Simulations...', value = 0, {
      sim_list <- vector("list", n_sims)
      for (i in 1:n_sims) {
        sim_list[[i]] <- run_single_sim(isolate(reactiveValuesToList(input)))
        incProgress(1 / n_sims)
      }
      bind_rows(sim_list, .id = "sim_id")
    })
  })

  output$imbalance_plot <- renderPlot({
    res <- mc_results()

    summary_df <- res %>%
      group_by(Scenario) %>%
      summarize(
        HbA1c_Imbalanced = mean(abs(HbA1c_SMD) > 0.1, na.rm = TRUE),
        LDL_Imbalanced   = mean(abs(LDL_SMD) > 0.1, na.rm = TRUE)
      ) %>%
      pivot_longer(cols = c("HbA1c_Imbalanced", "LDL_Imbalanced"), names_to = "Biometric", values_to = "Prop_Imbalanced") %>%
      mutate(
        Biometric = ifelse(Biometric == "HbA1c_Imbalanced", "HbA1c (Selectively Ordered)", "LDL (Routinely Measured)"),
        Pct_Label = sprintf("%.1f%%", Prop_Imbalanced * 100)
      )

    ggplot(summary_df, aes(x = Prop_Imbalanced, y = reorder(Scenario, desc(Scenario)), fill = Biometric)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
      geom_text(aes(label = Pct_Label), position = position_dodge(width = 0.9), hjust = -0.1, size = 4) +
      scale_x_continuous(labels = scales::percent, limits = c(0, 1.15)) +
      scale_fill_manual(values = c("HbA1c (Selectively Ordered)" = "firebrick", "LDL (Routinely Measured)" = "steelblue")) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom") +
      labs(
        title = paste("Proportion of Simulations with |SMD| > 0.1 Across", isolate(input$n_sims), "Runs"),
        subtitle = "Explicit Target Populations Labeled on Y-Axis",
        x = "Proportion of Runs Imbalanced (|SMD| > 0.1)",
        y = NULL,
        fill = "Biometric:"
      )
  })

  output$imbalance_table <- renderTable({
    mc_results() %>%
      group_by(Scenario) %>%
      summarize(
        `HbA1c Imbalance Rate (|SMD|>0.1)` = sprintf("%.1f%%", mean(abs(HbA1c_SMD) > 0.1, na.rm = TRUE) * 100),
        `HbA1c Mean SMD` = sprintf("%.3f", mean(HbA1c_SMD, na.rm = TRUE)),
        `LDL Imbalance Rate (|SMD|>0.1)` = sprintf("%.1f%%", mean(abs(LDL_SMD) > 0.1, na.rm = TRUE) * 100),
        `LDL Mean SMD` = sprintf("%.3f", mean(LDL_SMD, na.rm = TRUE))
      )
  })

  output$smd_box_plot <- renderPlot({
    res <- mc_results() %>%
      pivot_longer(cols = c("HbA1c_SMD", "LDL_SMD"), names_to = "Biometric", values_to = "SMD") %>%
      mutate(Biometric = ifelse(Biometric == "HbA1c_SMD", "HbA1c (Selectively Ordered)", "LDL (Routinely Measured)"))

    ggplot(res, aes(x = SMD, y = reorder(Scenario, desc(Scenario)), fill = Biometric)) +
      geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "red", size = 0.8) +
      geom_vline(xintercept = 0, linetype = "solid", color = "gray50") +
      geom_boxplot(alpha = 0.7, outlier.size = 1) +
      scale_fill_manual(values = c("HbA1c (Selectively Ordered)" = "firebrick", "LDL (Routinely Measured)" = "steelblue")) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom") +
      labs(
        title = "Distribution of Standardized Mean Differences (SMD)",
        subtitle = "Red dashed lines mark the [-0.1, 0.1] balance corridor",
        x = "Standardized Mean Difference (SMD)",
        y = NULL,
        fill = "Biometric:"
      )
  })
}

shinyApp(ui = ui, server = server)
