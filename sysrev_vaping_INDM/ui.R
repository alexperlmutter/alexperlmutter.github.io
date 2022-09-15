#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

ui <- fluidPage(
  titlePanel("Correction for independent nondifferential misclassification of nicotine vaping"),

  sidebarLayout(
    sidebarPanel(
      helpText("Change sensitivity and specificity of nicotine vaping to get corrected OR below"),

      sliderInput("Se",
                  label = "Sensitivity:",
                  min = min(dat$Se), max = max(dat$Se), value = 0.8, step = 0.05),

      sliderInput("Sp",
                  label = "Specificity:",
                  min = min(dat$Sp), max = max(dat$Sp), value = 0.99, step = 0.01),
    ),

    mainPanel(
      tableOutput("DataTable")
    )
  )
)
