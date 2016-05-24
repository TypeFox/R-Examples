#####  shinyElicit ui

shinyUI(fluidPage(
  tags$head(tags$script(HTML('
      Shiny.addCustomMessageHandler("jsCode",
        function(message) {
          console.log(message)
          eval(message.code);
        }
      );
    '))),
  uiOutput("debugTools"),
  h1("Biomarker validation study design support"),
  a(
    href="Using_the_NNTbiomarker_package.htm", rel="help", target="_blank",
    #href="../doc/Using_the_NNTbiomarker_package.html", rel="help", target="_blank",
    #href="http://www.github.org/professorbeautiful/NNTbiomarkerHome/man/Using_the_NNTbiomarker_package.html",

    fluidRow(column(1, offset=2, strong(em("click for information:",
                                           style="color:lightgreen")))
             ,
             column(3, actionButton(inputId = "Info", label="",
                                    style="background:lightgreen",
                                    icon=icon("info-sign", lib="glyphicon")))
    )
  ),
  hr(),
  a(href="Steps.html", target="_blank",
    actionButton(inputId = "reportButton",
                 label = "When all steps are Done, you can click here for a report. (In Progress)")
  ),
  div(style="background:darkGrey",
      span(checkboxInput(inputId='stepTableCheckbox', value=TRUE,
                    label=em(strong("View NNT design table of stepping stones"))),
           actionButton(inputId = "autoFill", label="autoFill for testing")
      ),
      conditionalPanel('input.stepTableCheckbox',
                       tableOutput("steps")
      )
  ),
  hr(),
  textInput(inputId = "biomarkerReportTitle", label = "Report Title", width = "100%"),
  div(style="position:relative;overflow:scroll;height:1200px;background:lightgrey",
      ##The sectionHeader function provides the text from stepsTableInitial
      sectionHeader(  #1
        div(style="vertical-align:middle;font-size:150%",
            HTML(stringr::str_dup("&nbsp;", 15)),
            tags$textarea(id = "ClinicalScenario", style="width:100%; rows:4"),
            br(),
            "A description of the patients for whom it would be better to treat",
            tags$textarea(id = "BestToTreatDescription", style="width:100%; rows:4")
        )),
      sectionHeader(  #2
        fluidRow(
          column(1, ""),
          column(5, h3("Intended beneficiaries"), tags$textarea(id = "who", style="width:100%; rows:4")),
          column(5,
                 #HTML(stringr::str_dup("&nbsp;", 15)),
                 textInput(inputId = "Option_Treat", width="100%",
                           label = "Name of the more active decision choice"),
                 HTML(stringr::str_dup("&nbsp;", 15)),
                 textInput(inputId = "Option_Wait", width="100%",
                           label = "Name of the more conservative decision choice")
        )
        )),
      sectionHeader(  #3
        div(
          HTML(stringr::str_dup("&nbsp;", 15)),
          h3("Set the NNT discomfort range (NNTlower and NNTupper)."),
          HTML(stringr::str_dup("&nbsp;", 15)),
          h3("The current NNT is 1/prevalence, often between NNTlower and NNTupper."),
          fluidRow(
            #column(0, HTML("&nbsp;")),
            #           column(4, sliderInput("NNTlower", label = "NNTlower",
            #                                  value=7, min = 1, max=100, step=1)),
            #           column(4, sliderInput("prevalence", label = "prevalence = Pr(BestToAct) = 1/NNT",
            #                                  value=0.1, min = 0, max=1, step=0.05)),
            #           column(4, sliderInput("NNTupper", label = "NNTupper",
            #                                  value=17, min = 10, max=100, step=1))
            column(4, numericInput("NNTlower", label = "NNTlower",
                                   value=7, min = 1,  step=1)),
            column(4, numericInput("prevalence", label = "prevalence = Pr(BestToAct) = 1/NNT",
                                   value=0.1, min = 0, max=1, step=0.01)),
            column(4, numericInput("NNTupper", label = "NNTupper",
                                   value=17, min = 2, step=1))
          ),
          plotOutput(outputId = "plotDiscomfort", height='200px'),
          fluidRow(
            column(2, HTML("&nbsp;")),
            column(5, numericInput("NNTpos", label = "NNTpos, should be smaller than NNTlower",
                                   value=2, min = 1, step=0.5)),
            column(5, numericInput("NNTneg", label = "NNTneg, should be larger than NNTupper",
                                   value=30, min = 1, step=0.5))
          ),
          fluidRow(
            column(2, HTML("&nbsp;")),
            column(5, "Positive predictive value = 1/NNTpos = ",
                   textOutput("PPVderived") ),
            column(5, "Negative predictive value = 1 - 1/NNTneg = ",
                   textOutput("NPVderived") )
          )
         # plotOutput(outputId = "plotNNTgoals", height='250px')
        )  ### End of div()
      ),  ### End of sectionHeader()
      sectionHeader(  #4
        div(style="vertical-align:middle;font-size:150%",
            HTML(stringr::str_dup("&nbsp;", 15)),
            h3("Specific clinical benefit hoped for:")
            , tags$textarea(id = "SpecificBenefit",
                            style="width: 1000px; height: 150px;", rows=4)
        )),
      sectionHeader(    #5
        div(
          fluidRow(
            column(6,
                   numericInput("NpatientsProspective", label = "Prospective study sample size",
                                value=30, min = 10, max=1000, step = 1) ),
            column(6,
                   numericInput("percentPositive", label = "Expected percent with a positive test",
                                value=50, min = 5, max=95, step = 5) )
          ), ### End of fluidRow()
          fluidRow(
            column(6, h3("Follow-up"),
                   tags$textarea("follow_up", value = "(follow-up period & plan)")
            ),
            #textOutput
            column(6, h3("Notes"),
                   tags$textarea("ProspectiveStudyNotes",
                                 value = "(any notes on the prospective study design parameters)")
            )
          )  ### End of fluidRow()
        )  ### end of div()
      ),  ### End of sectionHeader()
      sectionHeader( #6
        div(
          h3("Required sensitivity and specificity for a retrospective study."),
          h2("Select the prevalence; hover mouse on plot for calculations."),
#           numericInput("prevalence", label = "prevalence",
#                        value=0.5, min = 0, max=1, step = 0.05),
          fluidRow(column(6,
                          plotOutput("contraBayesPlot",
                                     click="contraBayesPlot_click",
                                     hover="contraBayesPlot_hover",
                                     width='100%')),
                   column(6, tableOutput("selectedNNTPosNeg"),
                          numericInput("samplesizeCases",
                                       label = "Retrospective study #cases",
                                       value=30),
                          numericInput("samplesizeControls",
                                       label = "Retrospective study #controls",
                                       value=30),
                          h3("Anticipated results (TODO)")
                   )
          )
        ))
  ) # end scroll pane
))
