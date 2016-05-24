#This is the user-interface portion of the interactive binomial distribution.

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("Introduction to the Binomial density"),

  # Sidebar with a slider input for number of bins
     mainPanel(
      p("The binomial density takes two parameters, the n parameter, which is the number of trials, and the p parameter,
        which is the probability of success for each trial. If we flipped a fair coin (which has the same probability of heads as tails) ten times,
        we would use a binomial(10,0.5) distribution, which is what is currently displayed to the right."),
     p("Try changing the n and p 
        parameters to see how the probabilities change."),
      plotOutput("pbinomialPlot"),
     sliderInput("singlep","p",min=0,max=1,value=0.5,step=0.01), 
     textOutput("singletext"),
     plotOutput("binomialPlot"),
     br(),
     br(),
     sliderInput("bothn","n",min=1,max=50,value=10),
     sliderInput("bothp","p",min=0,max=1,value=0.5,step=0.01),
     textOutput("bothtext")
    )
  
))
