library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Fitted Catenary"),
  
  
  sidebarPanel(
    
    numericInput("c1","C1:",0.4),
    numericInput("c2","C2:",2),
    numericInput("lambda","lambda:",3),
    numericInput("x0","left endpoint:",0),
    numericInput("x1","right endpoint:",4),
    numericInput("sd","standard deviation for simulation:",0.1),
    numericInput("R","Number of iterations:",10),
    checkboxInput("para_sim","Sim parabola",FALSE),
    conditionalPanel(
      condition = "input.para_sim == true",
      numericInput("a","a:",1),
      numericInput("b","b:",1),
      numericInput("c","c:",1)
    ),
    checkboxInput('aspect',"Fix aspect ratio",FALSE)
  ),

  
  mainPanel(
    plotOutput("catPlot"),
    h3("Fitted lines"),
    checkboxInput("cat","Catenary"),
    checkboxInput("para","Parabola"),
    verbatimTextOutput("summary")
  )
))