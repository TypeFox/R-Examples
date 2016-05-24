library(shiny)

# Define UI for dataset viewer application
shinyUI(pageWithSidebar(
  
  # Application title
  titlePanel("Data check: profiling a table against a set of rules"),
  
  # Sidebar with controls to select a dataset and specify the number
  # of observations to view
  sidebarPanel(
#     fileInput("dataset", "Input dataset file", multiple = FALSE, accept = NULL),
#     uiOutput("recLabels"),
#     fileInput("ruleset", "Input rule set file", multiple = FALSE, accept = NULL),
#     downloadButton('downloadData', 'Download score data'),
#     uiOutput("recScores"),
#     plotOutput("scoreSums"),
#     
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
      tabPanel("Vignette", list(
        div(id = "vignette", includeMarkdown(system.file("doc/index.md", package = "datacheck")))
      ))
#       ,
#       tabPanel("Table",list(
#         div(id='progress',includeHTML("js/timer.js")),
#         tableOutput("view")))
#       ,
#       tabPanel("Scores",list(
#         div(id='progress',includeHTML("js/timer.js")),
#         tableOutput("scores")))
#       ,
#       
#       tabPanel("Summaries",{
#         tableOutput("descriptive")  
#       })
#       ,      
#       tabPanel("Heatmap", p(
#         plotOutput("heatmap", height = "1000px"))
#         ),
#       tabPanel("Coverage", plotOutput("coverage", height = "1500px")),
#       tabPanel("Profile",tableOutput("profile"))
    )
  )
  
))

