shinyUI(fluidPage(
  titlePanel(h3("Upload the Computation Definition"))
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
                                       , actionButton("gotoDataInputs", "Continue"))
                        )
             , "-----"
             , tabPanel("Help-FAQ"
                      , h5("here in help")
                        )
               )
))


