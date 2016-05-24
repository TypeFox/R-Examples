shinyUI(fluidPage(
  titlePanel("Define the Stratified Cox Model")
, navlistPanel(selected = "Data Load"
             , id = 'navigationList' ## keeps track of our state for navigation panel on left
             , tabPanel(title = "Data Load"
                      , value = "dataLoad"
                      , selectInput(inputId = "input_type", label = "Data Input type",
                                    choices = c("CSV File", "Redcap API", "Postgres"))
                      , conditionalPanel(condition = "input.input_type == 'CSV File'"
                                       , textInput(inputId = "missingIndicators"
                                                 , label = "Optional missing values indicator(s), comma separated"
                                                 , value = "NA")

                                       , fileInput(inputId = 'dataFile'
                                                 , label = '(CSV data file)'
                                                 , multiple = FALSE
                                                 , accept = NULL))

                      , conditionalPanel(condition = "input.input_type == 'Redcap API'"
                                       , textInput(inputId = "redcapURL"
                                                 , label = "Redcap URL"
                                                 , value = "")

                                       , passwordInput(inputId = "redcapToken"
                                                     , label = "Redcap API Token"))

                      , conditionalPanel(condition = "input.input_type == 'Postgres'"
                                       , textInput(inputId = "dbName"
                                                 , label = "Database name"
                                                 , value = "")

                                       , textInput(inputId = "dbHost"
                                                 , label = "Database host"
                                                 , value = "")

                                       , textInput(inputId = "dbPort"
                                                 , label = "Database Port"
                                                 , value = "")

                                       , textInput(inputId = "dbUser"
                                                 , label = "Database User"
                                                 , value = "")

                                       , passwordInput(inputId = "dbPassword"
                                                     , label = "Database Password")

                                       , textInput(inputId = "dbTable"
                                                 , label = "Database Table Name"
                                                 , value = ""))

                      , actionButton(inputId = "loadData"
                                   , label = "Load Data")

                      , conditionalPanel(condition = "input.loadData != 0"
                                       , h5("Summary")
                                       , verbatimTextOutput(outputId = 'dataFileContentSummary')
                                       , br()
                                       , h3(textOutput(outputId = 'dataLoaded')))
                        )

             , tabPanel(title = "Formula Check"
                      , value = "formulaCheck"
                      , textInput(inputId = "formula"
                                , label = "Formula"
                                , value = "")
                      , actionButton(inputId = "checkFormula",
                                     label = "Check Formula")
                      , br()
                      , conditionalPanel(condition = "input.checkFormula != 0"
                                       , h5("Summary")
                                       , verbatimTextOutput(outputId = 'checkFormulaResult')
                                       , br()
                                       , h3(textOutput(outputId = 'formulaChecked')))
                        )

             , tabPanel(title = "Output Result"
                      , value = "outputResult"
                      , textInput(inputId = "outputFile"
                                , label = "Output File Name"
                                , value = "defn.rds")
                      , actionButton(inputId = "saveDefinition"
                                   , label = "Save Definition")
                      , conditionalPanel(condition = "input.saveDefinition != 0"
                                       , verbatimTextOutput(outputId = 'outputResult')
                                       , br()
                                       , h5(textOutput(outputId = 'definitionSaved'))
                                       , br()
                                       , actionButton(inputId = "exitApp", label = "Exit"))
                        )
             , "-----"
             , tabPanel("Help-FAQ"
                      , h5("here in help")
                        )
               )
))

