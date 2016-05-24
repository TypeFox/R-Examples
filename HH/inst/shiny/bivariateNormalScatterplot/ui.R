library(shiny)
library(lattice)

## Define UI for application
shinyUI(pageWithSidebar(
  headerPanel("Bivariate Normal at Various Correlations"),
  ## Sidebar with a slider inputs
  sidebarPanel(
    sliderInput("rho", "rho",     -1, 1, .35, .05, animate=list(loop=TRUE, interval=500)),
    numericInput("seed","seed", 1234)
  ),
  mainPanel(
    plotOutput("correlationPlot")
  )
  ))
