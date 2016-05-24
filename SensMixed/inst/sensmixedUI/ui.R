shinyUI(fluidPage(
  ## for different themes install shinythemes
  #theme = shinytheme("spacelab"),
  #theme = "united.min.css",
  #theme = "style.css",
  # tags$head(
  #   tags$link(rel = "stylesheet", type = "text/css", href = "united.min.css")
  # ),

  div(titlePanel("Analysis of Sensory and Consumer data within a mixed effects model framework"), style = "color:#191970"),
  helpText("This application is a user-friendly interface for
      the R-package SensMixed"),
  #sidebarPanel(
  fluidRow(
    column(4,
#          wellPanel(
#           h4("Choose type of analysis"),
#           radioButtons("analysis", "Analysis of", 
#                        choices = c("Sensory data", "Consumer data"),
#                        inline = TRUE)
#          ),
           uiOutput("antypeUI"),
#          bsTooltip("analysis", "title", placement = "bottom", trigger = "hover"),
          uiOutput("AttrUI"),
          #submitButton("Run Analysis")
         # actionButton("goButton", "Run Analysis"),
        bsButton("goButton", label = "Run Analysis", type = "action", 
                 style = "primary")
    ),
    column(8,
    #       tabsetPanel(
      uiOutput("theTabset")
#       tabPanel("Data",
#                h4("Choose data"),
#                selectInput('uploaddata', '', 
#                            c("Read CSV file from local drive" = 1, 
#                              "TVbo data" = 2, "Ham data" = 3)),
#                uiOutput("UploadUI")
#       ),
#       uiOutput("plotUI"),
#       tabPanel("Plot output",
#                helpText("Note: This output is only dedicated for analysis of sensory data"),
#                conditionalPanel(
#                condition =  "input.analysis == 'Sensory data'",
#                sidebarLayout(
#                  sidebarPanel(
#                    selectInput('typeEffs', 'Plot effects', 
#                                c("random" = 1, "fixed" = 2, "scaling" = 3)),
#                    selectInput('typePlot', 'Plot type', 
#                                c("F" = FALSE, "d-prime" = TRUE)),
#                    selectInput('representPlot', 'Layout', 
#                                c("single" = FALSE, "multiple" = TRUE)),
#                    numericInput('scalePlot', label = "Scale plot", value = 1),
#                    downloadButton('downloadPlot', label = "Download Plot")                   
#                  ),                  
#                  mainPanel(
#                    plotOutput("plotsSensMixed")
#                  )
#                 
#                ),               
#                value = 1)),
#       tabPanel("Table output",
#                helpText("Note: This output is only dedicated for analysis of sensory data"),
#                sidebarLayout(
#                  sidebarPanel(
#                    selectInput('typeEffsTable', 'Type of effects', 
#                                c("random" = 1, "fixed" = 2, "scaling" = 3, 
#                                  "all" = 4)),
#                    selectInput("typetable2", "Type", c("html", "latex")),
#                    downloadButton('downloadTable', label = "Download Table")
#                  ),
#                  mainPanel(
#                    htmlOutput("tablesSensMixed")
#                    )
#                ), 
#                value = 2),
#       tabPanel("Step output",               
#                sidebarLayout(
# #                  conditionalPanel(condition =  "input.analysis == 'Sensory data'",
# #                  sidebarPanel(
# #                    uiOutput("AttrStepUI"))), 
#                  sidebarPanel(
#                    uiOutput("AttrStepUI")),
#                  mainPanel(
#                    htmlOutput("stepRand"), 
#                    br(),
#                    htmlOutput("stepFixed")
#                  )
#                ),              
#                value = 3),
#       tabPanel("Post-hoc",
#                sidebarLayout(
#                  sidebarPanel(
#                    uiOutput("AttrPosthocUI"),
#                    uiOutput("EffsPosthocUI")),                   
#                  mainPanel(
#                    plotOutput("posthocPlot"),
#                    htmlOutput("posthocTable")
#                  )
#                ),              
#                value = 4),
#       id="tabs1")
  )
  )
  ))
