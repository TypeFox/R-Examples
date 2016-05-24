library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Catenary"),
  
  
  sidebarPanel(
    numericInput("c1","Choose a value for C1:",1),
    numericInput("c2","Choose a value for C2:",1),
    numericInput("lambda","Choose a value for lambda:",1),
    numericInput("x0","Choose a value for left endpoint:",-2),
    numericInput("x1","Choose a value for right endpoint:",2)
  ),

  
  mainPanel(
    plotOutput("catPlot"),
    verbatimTextOutput("summary")
  )
))