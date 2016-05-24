library(shiny)
library(cbsodata)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("CBS OData"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      selectizeInput("table"
                    , "Table:"
                    , c("Loading...")
                    , options = list()
                    )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      dataTableOutput("table_list")
    )
  )
))