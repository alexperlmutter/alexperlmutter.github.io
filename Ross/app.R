library(shiny)
library(dplyr)
library(ggplot2)
library(broom)

expit <- function(x) 1 / (1 + exp(-x))

ui <- fluidPage(
  titlePanel("Ross Framework: Causal Effect Estimation under Selective Testing & Transportability"),

  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("1. Population & True Causal Effect"),
      numericInput("n_pop", "Total Population Size (N):", value = 6000, min = 2000, max = 20000, step = 1000),
      sliderInput("true_theta_a", "True Treatment Effect on Outcome Y (Theta_A):", min = -2.0, max = 2.0, value = 0.8, step = 0.1),

      hr(),
      h4("2. Transportability Bias (S=1 vs S=0)"),
      sliderInput("or_s_x", "How many times MORE LIKELY are high-risk patients to be in UNLINKED Target Claims vs LINKED EHR? [OR]:",
                  min = 0.2, max = 5.0, value = 3.0, step = 0.2),
      checkboxInput("misspec_s", "Misspecify SMR Transport Model (Omit non-linear X^2)?", value = FALSE),

      hr(),
      h4("3. Selective Test Ordering (O=1 vs O=0)"),
      sliderInput("or_o_a", "How many times MORE LIKELY are TREATED patients (A=1) to get tested vs untreated? [OR]:",
                  min = 1.0, max = 6.0, value = 3.5, step = 0.5),
      sliderInput("or_o_z", "How many times MORE LIKELY are patients with worse Biometrics (Z) to get tested? [OR]:",
                  min = 1.0, max = 5.0, value = 2.5, step = 0.5),

      hr(),
      h4("4. Biometric Confounding Strength"),
      sliderInput("z_conf_y", "Impact of Biometric (Z) on Outcome Y per SD increase:", min = 0.0, max = 2.0, value = 1.0, step = 0.1),

      actionButton("resimulate", "Re-run Causal Pipeline", class = "btn-primary", style = "width: 100%; margin-top: 10px;")
    ),

    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Treatment Effect Recovery",
                 plotOutput("forest_plot", height = "450px"),
                 hr(),
                 tableOutput("summary_table")),

        tabPanel("Weight Diagnostics & ESS",
                 plotOutput("weight_dist_plot", height = "350px"),
                 hr(),
                 tableOutput("ess_table"))
      )
    )
  )
)

server <- function(input, output, session) {

  sim_data <- reactive({
    input$resimulate

    N <- input$n_pop

    # Convert intuitive ORs to log-odds parameters
    alpha_x <- log(input$or_s_x)
    delta_a <- log(input$or_o_a)
    delta_z <- log(input$or_o_z)

    # 1. Claims Covariate X
    X <- rnorm(N, 0, 1)

    # 2. Selection into Target Claims (S=1) vs Linked EHR (S=0)
    S_prob <- expit(-0.5 + alpha_x * X + 0.4 * (X^2))
    S <- rbinom(N, 1, S_prob)

    # 3. Biometric Confounder Z
    Z <- 0.6 * X + rnorm(N, 0, 1)

    # 4. Treatment Assignment A
    A_prob <- expit(-0.2 + 0.7 * X + 0.9 * Z)
    A <- rbinom(N, 1, A_prob)

    # 5. Selective Test Ordering O
    O_prob <- expit(-0.5 + 0.4 * X + delta_a * A + delta_z * Z)
    O <- rbinom(N, 1, O_prob)

    # 6. Continuous Outcome Y
    Y <- 1.0 + input$true_theta_a * A + 0.7 * X + input$z_conf_y * Z + rnorm(N, 0, 1)

    data.frame(id = 1:N, X = X, Z = Z, A = A, S = S, O = O, Y = Y)
  })

  results <- reactive({
    df <- sim_data()

    # --------------------------------------------------------------------------
    # Stage 1: Transportability Model (W_SMR)
    # --------------------------------------------------------------------------
    if (input$misspec_s) {
      smr_mod <- glm(S ~ X, data = df, family = binomial())
    } else {
      smr_mod <- glm(S ~ X + I(X^2), data = df, family = binomial())
    }

    df$ps_smr <- pmin(pmax(predict(smr_mod, type = "response"), 0.01), 0.99)
    df$W_SMR <- ifelse(df$S == 0, df$ps_smr / (1 - df$ps_smr), 1.0)

    # --------------------------------------------------------------------------
    # Stage 2: Test Ordering Model (W_ORD) on Linked Cohort S=0
    # --------------------------------------------------------------------------
    df_linked <- df %>% filter(S == 0)
    ord_mod <- glm(O ~ X + A, data = df_linked, family = binomial())

    df_linked$ps_ord <- pmin(pmax(predict(ord_mod, type = "response"), 0.01), 0.99)
    df_linked$W_ORD <- 1 / df_linked$ps_ord

    df <- df %>% left_join(df_linked %>% select(id, W_ORD), by = "id")
    df_cc <- df %>% filter(S == 0 & O == 1)

    # --------------------------------------------------------------------------
    # Stage 3: Biometric Treatment Propensity Models
    # --------------------------------------------------------------------------
    # Unweighted Naive Model
    tx_mod_naive <- glm(A ~ X + Z, data = df_cc, family = binomial())
    df_cc$ps_tx_naive <- predict(tx_mod_naive, type = "response")
    df_cc$W_TX_naive <- ifelse(df_cc$A == 1, 1 / df_cc$ps_tx_naive, 1 / (1 - df_cc$ps_tx_naive))

    # SMR Weighted Only (Ignores W_ORD)
    tx_mod_smr <- glm(A ~ X + Z, data = df_cc, weights = W_SMR, family = binomial())
    df_cc$ps_tx_smr <- predict(tx_mod_smr, type = "response")
    df_cc$W_TX_smr <- ifelse(df_cc$A == 1, 1 / df_cc$ps_tx_smr, 1 / (1 - df_cc$ps_tx_smr))

    # Full Ross Approach A (Weighted by W_SMR * W_ORD)
    df_cc$W_OBS <- df_cc$W_SMR * df_cc$W_ORD
    tx_mod_full <- glm(A ~ X + Z, data = df_cc, weights = W_OBS, family = binomial())
    df_cc$ps_tx_full <- predict(tx_mod_full, type = "response")
    df_cc$W_TX_full <- ifelse(df_cc$A == 1, 1 / df_cc$ps_tx_full, 1 / (1 - df_cc$ps_tx_full))

    # Composite Final Weights
    df_cc <- df_cc %>%
      mutate(
        W_final_naive = W_TX_naive,
        W_final_smr_only = W_SMR * W_TX_smr,
        W_final_ord_only = W_ORD * W_TX_full,
        W_final_full = W_SMR * W_ORD * W_TX_full
      )

    # --------------------------------------------------------------------------
    # Stage 4: Outcome Models
    # --------------------------------------------------------------------------
    m_oracle <- lm(Y ~ A + X + Z, data = df)
    m_naive  <- lm(Y ~ A, data = df_cc)
    m_no_ord <- lm(Y ~ A, data = df_cc, weights = W_final_smr_only)
    m_no_smr <- lm(Y ~ A, data = df_cc, weights = W_final_ord_only)
    m_full   <- lm(Y ~ A, data = df_cc, weights = W_final_full)

    res_table <- bind_rows(
      tidy(m_oracle, conf.int = TRUE) %>% filter(term == "A") %>% mutate(Model = "1. Oracle Target (Full Data)"),
      tidy(m_naive, conf.int = TRUE) %>% filter(term == "A") %>% mutate(Model = "2. Naive Complete-Case (Unweighted)"),
      tidy(m_no_ord, conf.int = TRUE) %>% filter(term == "A") %>% mutate(Model = "3. Ignores Ordering Bias (SMR Weight Only)"),
      tidy(m_no_smr, conf.int = TRUE) %>% filter(term == "A") %>% mutate(Model = "4. Ignores Transportability (ORD Weight Only)"),
      tidy(m_full, conf.int = TRUE) %>% filter(term == "A") %>% mutate(Model = "5. Full 3-Stage Pipeline (Ross Approach A)")
    ) %>%
      select(Model, estimate, std.error, conf.low, conf.high)

    list(estimates = res_table, weights_df = df_cc)
  })

  output$forest_plot <- renderPlot({
    res <- results()$estimates
    true_val <- input$true_theta_a

    ggplot(res, aes(x = estimate, y = reorder(Model, estimate))) +
      geom_vline(xintercept = true_val, linetype = "dashed", color = "red", size = 1) +
      geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.8, color = "navy") +
      theme_minimal(base_size = 14) +
      labs(
        title = "Treatment Effect Recovery Across Pipelines",
        subtitle = "Red dashed line indicates the True Causal Effect (Theta_A)",
        x = "Estimated Treatment Effect (Beta_A)",
        y = NULL
      )
  })

  output$summary_table <- renderTable({
    results()$estimates %>%
      mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
      rename(`Model Pipeline` = Model, Estimate = estimate, `Std Error` = std.error, `95% CI Lower` = conf.low, `95% CI Upper` = conf.high)
  })

  output$weight_dist_plot <- renderPlot({
    w_df <- results()$weights_df
    ggplot(w_df, aes(x = W_final_full)) +
      geom_histogram(bins = 40, fill = "teal", color = "black", alpha = 0.7) +
      theme_minimal() +
      labs(title = "Distribution of Composite Final Weights (W_SMR * W_ORD * W_TX)", x = "Composite Weight Value", y = "Frequency")
  })

  output$ess_table <- renderTable({
    w_df <- results()$weights_df
    calc_ess <- function(w) sum(w)^2 / sum(w^2)

    data.frame(
      Weight_Type = c("Raw SMR Weight (W_SMR)", "Raw Ordering Weight (W_ORD)", "Full Composite Weight"),
      N_Complete_Cases = nrow(w_df),
      Effective_Sample_Size = c(calc_ess(w_df$W_SMR), calc_ess(w_df$W_ORD), calc_ess(w_df$W_final_full)),
      Max_Weight = c(max(w_df$W_SMR), max(w_df$W_ORD), max(w_df$W_final_full))
    ) %>%
      mutate(ESS_Percent = paste0(round((Effective_Sample_Size / N_Complete_Cases) * 100, 1), "%"))
  })
}

shinyApp(ui, server)
