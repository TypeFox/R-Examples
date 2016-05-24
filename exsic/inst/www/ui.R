
# Define UI for dataset viewer application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Exsic DEMO: formatting exsiccatae"),
  
  # Sidebar with controls to select a dataset and specify the number
  # of observations to view
  sidebarPanel(
    helpText(""),
    selectInput("format", "Choose a format:", 
                choices = c("SBMG", "ASPT", "NYBG", "PK")),
    
    helpText(""),
    
    helpText("Any suggestion or questions should be directed to:"),
    helpText("r.simon@cgiar.org")
    
  ),
  
  # Show a summary of the dataset and an HTML table with the requested
  # number of observations
  # for progress timer see: https://gist.github.com/jcheng5/4495659 
  mainPanel(
    #verbatimTextOutput("summary"),
    tabsetPanel(
      tabPanel("Help", ""
        #list(
        #div("help", #includeHTML("vignette.html")
        #    ))
               ),
      tabPanel("Exsiccatae", htmlOutput("exsic")),
        
      tabPanel("Data table",{
        tableOutput("data") 
      })
    )
  )
  
))
