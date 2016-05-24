#This is the user interface for the F distribution shiny app.

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("Interactive F Density"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("df1","DF 1",min=1,max=50,value=1),
      sliderInput("df2","DF 2",min=1,max=50,value=1)
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("fPlot")
    )
  )
))
