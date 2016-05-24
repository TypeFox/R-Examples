library(shiny)
library(HH)

## Define UI for application
shinyUI(pageWithSidebar(
headerPanel("Bivariate Normal Density"),
    ## Sidebar with a slider inputs
  sidebarPanel(
    sliderInput("rho", "rho",     -.80, .80, -.80, .05, animate=list(loop=TRUE, interval=500)),
    sliderInput("angle", "angle in degrees",  0, 360, 112.5, 22.5, animate=list(loop=TRUE))),
  mainPanel(
    plotOutput("densityPlot")
  )
  ))
