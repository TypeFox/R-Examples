library(shiny)
library(DT)


shinyUI(fluidPage(
  
  sidebarLayout(
    sidebarPanel(
      includeMarkdown("readme.md"),
      pre(includeText("prots.txt")),
      uiOutput("dynamic_ui")
    ),
    
    mainPanel(
      uiOutput("dynamic_tabset")    
    )
  )))
