shinyUI(fluidPage(
  titlePanel(h3("Output Master code"))
, navlistPanel(selected="Definition File", id='mynavlist'
             , tabPanel(title = "Definition File"
                      , value = "definitionUpload"
                      , fileInput(inputId = 'definitionFile'
                                , label = 'RDS file:'
                                , multiple = FALSE
                                , accept = c(".RDS", ".rds"))

                      , actionButton(inputId = "uploadDefinition"
                                   , label = "Upload Definition")

                      , conditionalPanel(condition = "input.uploadDefinition != 0"
                                       , verbatimTextOutput(outputId = 'definitionSummary')
                                       , br()
                                       , h5(textOutput(outputId = 'definitionUploaded')))
                        )
             , tabPanel(title = "Specify Sites"
                      , value = "specifySiteURLs"
                      , textInput(inputId = "siteName"
                                , label = "Site Name"
                                , value = "Site[0-9]+")
                      , textInput(inputId = "ocpuURL"
                                , label = "OpenCPU URL"
                                , value = "http://localhost:")
                      , actionButton(inputId = "addSite"
                                   , label = "Add site")
                      , conditionalPanel(condition = "input.addSite != 0"
                                       , verbatimTextOutput(outputId = 'siteList')
                                       , h5(textOutput(outputId = 'writeRCode')))
                        )
             , tabPanel(title = "Output R code"
                      , value = "outputResult"
                      , textInput(inputId = "outputFile"
                                , label = "Output File Name"
                                , value = "master.R")
                      , actionButton(inputId = "saveCode"
                                   , label = "Save Code")
                      , conditionalPanel(condition = "input.saveCode != 0"
                                       , h5(textOutput(outputId = 'codeSaved'))
                                       , br()
                                       , actionButton(inputId = "exitApp", label = "Exit"))
                        )
             , "-----"
             , tabPanel("Help-FAQ"
                      , h5("here in help")
                        )
               )
))


