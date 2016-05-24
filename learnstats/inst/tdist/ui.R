
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("Interactive Student's T Density"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      #sliderInput("bins",
       #           "Number of bins:",
        #          min = 1,
         #         max = 50,
          #        value = 30),
      sliderInput("degrees","Move the slider to change the degrees of freedom",min=1,max=100,value=1)
    ),

    # Show a plot of the generated distribution
    mainPanel(
      helpText("Though the t with one degree of freedom looks different from the Normal(0,1), as the degrees of freedom increase, it becomes very close to a standard Normal"),
      plotOutput("thePlot"),
      helpText("In all circumstances, the tails of the T distribution are \"fatter,\" which means outliers are more likely under the T distribution."),
      plotOutput("secondPlot"),
      br(),
      br(),
      helpText("© Daniel Walter 2014")
      )
  )
))
