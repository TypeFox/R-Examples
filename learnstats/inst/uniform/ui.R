#UI for a Uniform distribution

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("Interactive Uniform Density"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("range","Range of the uniform density",min=-5,max=10,value=c(0,1),step=1)
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("unifPlot")
    )
  )
))
