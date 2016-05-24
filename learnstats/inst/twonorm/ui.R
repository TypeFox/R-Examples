
library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("An interactive introduction to the Normal Density"),

  # Sidebar with a slider input for number of bins
    p("The normal distribution is governed by two parameters. The first is the mean, or average
      of the distribution. Changing the mean changes the location of the probability density. 
      Below is a normal distribution that changes based on the mean you choose."),  
    plotOutput("meanPlot"),
    sliderInput("meanonly","Mean",min=-50,max=50,value=0,width="45%"),

    p("The second parameter is called the standard deviation. The standard deviation is dependent on measure 
      the squared distance from the mean. You can think of it as a measure of how close or spread out the 
      density is. Try changing the standard deviation to see how the density changes."),
   plotOutput("sdPlot"),
   sliderInput("sdonly","Standard Deviation",min=5,max=25,value=10,step=0.5,width="45%"),

    p("In real life, though, both parameters can be changed. Here are two normal distributions 
      that you can change both parameters for."),
        # Show a plot of the generated distribution
      plotOutput("twonormPlot"),
      fluidRow(column(5,offset=1,sliderInput("onemean","Red Mean",min=-50,max=50,value=-3)),
               column(5,sliderInput("twomean","Blue Mean",min=-50,max=50,value=4))),
      fluidRow(column(5,offset=1,sliderInput("onesdev","Red Standard deviation",min=4,max=25,value=5,step=0.5)),
               column(5,sliderInput("twosdev","Blue Standard deviation",min=4,max=25,value=10,step=0.5)))  
  
))
