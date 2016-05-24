library(shiny)

source("chooser.R")

shinyUI(navbarPage(
  header=includeScript("www/js/jquery-ui.custom.min.js"),
  "OpenMx IFA Model Builder",
  tabPanel(
    "Data",
    sidebarLayout(
      sidebarPanel(
        fileInput('file1', 'Choose CSV File',
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
        tags$hr(),
        checkboxInput('dataHeader', 'Header?', TRUE),
        checkboxInput('dataRowNames', 'Row names?', TRUE),
        radioButtons('dataSep', 'Separator',
                     c(Comma=',',
                       Semicolon=';',
                       Tab='\t',
                       Whitespace=' '),
                     ','),
        radioButtons('dataQuote', 'Quote',
                     c(None='',
                       'Double Quote'='"',
                       'Single Quote'="'"),
                     '"'),
        tags$hr(),
        helpText("Load example data"),
        actionButton("exampleDataKCT", label = "KCT"),
        actionButton("exampleDataLSAT6", label = "LSAT6"),
        actionButton("exampleDataScience", label = "Science")
      ),
      mainPanel(
        tabsetPanel(id="dataPreviewTabset",
                    tabPanel("First 6 rows", value="front",
                             tags$p("Data file:"),
                             verbatimTextOutput("nameOfDataFile"),
                             tags$p("Number of rows:"),
                             verbatimTextOutput("numberOfDataRows"),
                             hr(),
                             helpText("Unparsed content (first 6 lines):"),
                             verbatimTextOutput('unparsedDataContents'),
                             helpText("Parsed content (first 6 lines):"),
                             textOutput("parseFileFeedback"),
                             tableOutput('dataContents'),
                             helpText("Only the first 6 rows are shown to give",
                                      "you an idea of whether the data loaded",
                                      "correctly.")),
                    tabPanel("Item summary",
                             selectInput("freqColumnName", label = "Row frequency column:",
                                         choices="No data loaded"),
                             hr(),
                             tableOutput('dataSummary'),
                             helpText("The number of outcomes listed here",
                                      "are derived solely from the data and",
                                      "are not affected by subsequent recoding.")))
      )
    )
  ),
  tabPanel("Outcomes",
           sidebarLayout(
             sidebarPanel(
               selectInput("focusedOutcomeSet", label = "Outcome set:",
                           choices="No data loaded"),
               selectInput("focusedOutcomeItem", label = "Item:",
                           choices="No data loaded"),
               tags$hr(),
               textInput("newOutcomeName", label = "Add outcome"),
               actionButton("addNewOutcomeAction", label = "Add It"),
               textOutput("addNewOutcomeActionFeedback"),
               tags$hr(),
               selectInput("focusedOutcomeMapFrom", label="Recode from",
                           choices="No outcomes loaded"),
               selectInput("focusedOutcomeMapTo", label="to",
                           choices="No outcomes loaded"),
               conditionalPanel('input.focusedOutcomeMapTo == "<Rename>"',
                                textInput("focusedOutcomeRenameTo", label = "New name")),
               actionButton("focusedOutcomeMapAction", label = "Add mapping"),
               textOutput("focusedOutcomeMapActionFeedback"),
               tags$hr(),
               numericInput("focusedRecodeRule", "Recode Rule", value=1, step=1),
               actionButton("resetRecodeAction", label = "Discard"),
               textOutput("resetRecodeActionFeedback")
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Recode",
                          helpText("Currently selected outcomes:"),
                          tableOutput("focusedOutcomeTable"),
                          helpText("Recode Table:"),
                          tableOutput("recodeTable")
                 ),
                 tabPanel("Reorder",
                          wellPanel(helpText(
                            "Drag to reorder.",
                            "The standard order is from incorrect (upper) ",
                            "to correct (lower) or from least (upper) to",
                            "most (lower). Exceptional items that use the",
                            "opposite order can be marked as reversed",
                            "on the next tab."),
                            uiOutput('reorderOutcomesSorterUI')),
                          helpText("Permutation Table:"),
                          tableOutput("permuteTable")
                 ),
                 tabPanel("Reverse",
                          helpText("Items on the right side will be reverse scored."),
                          uiOutput("reversePicker")
                 )
               )
             ))),
  tabPanel("Model",
           sidebarLayout(
             sidebarPanel(
               selectInput("focusedItemStart", label = "Edit items from:",
                           choices="No data loaded"),
               selectInput("focusedItemEnd", label = "to:",
                           choices="No data loaded"),
               actionButton("selectAllItemsAction", label = "Focus all items"),
               tags$hr(),
               selectInput("focusedItemModel", label = "Model:",
                           choices=c('drm', 'grm', 'nrm')),
               selectInput("focusedItemModelTc", label = "Nominal Tc:",
                           choices=c('as is')),
               textOutput("focusedItemModelTcFeedback"),
               selectInput("focusedItemParameter", label = "Parameter:",
                           choices="No item selected"),
               selectInput("focusedParameterFree", label = "Free",
                           choices="No parameter selected"),
               textInput("focusedParameterLabel", label = "Label"),
               actionButton("changeLabelAction", label = "Set Label"),
               sliderInput("focusedParameterPrior", label = "Prior mode",
                           min=2, max=25, value=3),
               actionButton("focusedParameterPriorSetAction", label = "Set"),
               actionButton("focusedParameterPriorClearAction", label = "Clear"),
               textOutput("focusedParameterPriorFeedback")
             ),
             mainPanel(tabsetPanel(
               tabPanel("Factors",
                        sliderInput("numFactors", "Number of factors:",
                                    min=0, max=5, value=1),
                        textInput("nameOfFactor1", "Factor 1"),
                        textInput("nameOfFactor2", "Factor 2"),
                        textInput("nameOfFactor3", "Factor 3"),
                        textInput("nameOfFactor4", "Factor 4"),
                        textInput("nameOfFactor5", "Factor 5")),
               tabPanel("Reorder",
                        uiOutput('reorderItemsSorterUI')),
               tabPanel("Parameters",
                        helpText("Starting values"),
                        tableOutput('itemStartingValuesTable'),
                        helpText("Is free?"),
                        tableOutput('itemFreeTable'),
                        helpText("Labels"),
                        tableOutput('itemLabelTable'),
                        helpText("Bayesian prior mode"),
                        tableOutput('itemPriorTable')),
               tabPanel("Exclude",
                        helpText("Choose which items to exclude (if any).",
                                 "Items on the right side will be excluded."),
                        uiOutput("excludePicker")),
               tabPanel("Summary", tableOutput('itemModelAssignment'))
             ))
           )),
  tabPanel("Settings",
           sidebarPanel(
             helpText("You can store all the settings in a separate file",
                      "to reproduce the current analysis or reuse with similar analyses."),
             actionButton("refreshSettingsAction", label = "Refresh!"),
             tags$hr(),
             downloadButton('downloadCoding', 'Download'),
             tags$hr(),
             fileInput('codingFile', paste('Choose a setting save file from',
                                           'which to restore the settings'),
                       accept=c('text/plain', '.R')),
             verbatimTextOutput("codingFileFeedback")),
           mainPanel(
             verbatimTextOutput("debugSettingsOutput")
           )
  ),
  tabPanel("Analysis", sidebarLayout(
    sidebarPanel(
      selectInput("boundPriorForm", "Functional form for dichotomous bound prior density",
                  c("Logit-normal", "Beta")),
      checkboxInput("showFitProgress", label = "Show model fitting progress", value = TRUE),
      checkboxInput("fitReferenceModels", label = "Fit reference models (for more fit statistics)",
                    value = FALSE),
      selectInput("infoMethod", "Information matrix method:", 
                  choices = c("Oakes", "Meat", "Agile SEM", "*none*")),
      actionButton("debugScriptAction", label = "Refresh!"),
      tags$hr(),
      downloadButton('downloadScript', 'Download'),
      withTags(ol(
        li('Download your analysis script.'),
        li('Open it in RStudio.'),
        li('Update the pathname to your data (if necessary).'),
        li('Click the Knit/HTML button at the top of your document.'),
        li('Note: RStudio does not render numbers or equations correctly. Open the HTML in a regular browser.')
      ))
    ),
    mainPanel(
      verbatimTextOutput("debugScriptOutput")
    )
  ))

))
