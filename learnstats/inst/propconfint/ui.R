
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("Confidence Interval for Population Proportion"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("level",
                  "Confidence level: choose a value from 50 to 99",
                  min = 50,
                  max = 99,
                  value = 95),
      sliderInput("samsize",
                  "How large would you like your samples to be?",
                  min=10,
                  max=1000,
                  value=40),
      sliderInput("repnum",
                  "How many confidence intervals would you like to generate?",
                  min=10,
                  max=200,
                  value=20)
    ),
    # Show a plot of the generated distribution
    mainPanel(
      p("In this case, we're trying to estimate the population proportion. 
        You and I know that the true value is 0.6, and all of the data are 
        drawn from Binomial(n,p) samples, where p = 0.6.")
      ,p("But the people 
        conducting the experiment don't know this, so they're trying 
        to estimate the true population parameter from the samples they have. 
         They create confidence intervals from these estimates. "),
      plotOutput("confintPlot")
    )
  )
))
