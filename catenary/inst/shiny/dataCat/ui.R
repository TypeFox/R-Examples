library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Ctesiphon catenary"),
  
  
  sidebarPanel(
    h4("Left endpoint"),
    selectInput("internal","Select observations",
                c("Internal"="int",
                  "External"='ext'))),

  
  mainPanel(
    plotOutput("catPlot"),
    verbatimTextOutput("summary")
  )
))