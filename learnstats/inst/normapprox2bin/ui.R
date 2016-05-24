library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("The normal approximation to the binomial."),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("n",
                  "Number of Experiments:",
                  min = 1,
                  max = 750,
                  value = 30),
      sliderInput("p",
                  "Probability of Success:",
                  min = 0.01,
                  max = 0.99,
                  value = 0.5)
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("normbinPlot")
    )
  )
))
