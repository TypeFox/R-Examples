#####  shinyCombined ui

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
  fluidRow(
    column(7,
           h1("Biomarker validation study design support")),
    column(5, a(
      href="Using_the_NNTbiomarker_package.htm", rel="help", target="_blank",
      #href="../doc/Using_the_NNTbiomarker_package.html", rel="help", target="_blank",
      #href="http://www.github.org/professorbeautiful/NNTbiomarkerHome/man/Using_the_NNTbiomarker_package.html",
      fluidRow(
        column(3,
               style="background:yellow",
               strong(em("Click for information:",
                         style="color:darkgreen")))
        ,
        column(1, style="background:yellow",
               actionButton(inputId = "Info", label="",
                            style="background:lightgreen",
                            icon=icon("info-sign", lib="glyphicon")))
      )
    )
    )
  ),
  hr(),
  fluidRow(style="text-align:center",
           #           column(4, sliderInput("NNTlower", label = "NNTlower",
           #                                  value=7, min = 1, max=100, step=1)),
           #           column(4, sliderInput("prevalence", label = "prevalence = Pr(BestToAct) = 1/NNT",
           #                                  value=0.1, min = 0, max=1, step=0.05)),
           #           column(4, sliderInput("NNTupper", label = "NNTupper",
           #                                  value=17, min = 10, max=100, step=1))
           column(1, HTML("&nbsp;")),
           column(2, numericInput("NNTpos", label = "",
                                  value=2, min = 1, step=0.1)),
           column(2, numericInput("NNTlower", label = "",
                                  value=7, min = 1,  step=0.1)),
           column(2, numericInput("prevalence", label = "",
                                  value=0.1, min = 0, max=1, step=0.01)),
           column(2, numericInput("NNTupper", label = "",
                                  value=17, min = 2, step=0.1)),
           column(2, numericInput("NNTneg", value=30, label = "", min = 1, step=0.1))
  ),
  fluidRow(style="font-weight:bold; text-align:center",
           column(1, HTML("&nbsp;")),
           column(2, "NNTpos"),
           column(2, "NNTlower"),
           column(2, "prevalence"),
           column(2, "NNTupper"),
           column(2, "NNTneg")
  ),
  fluidRow(style="font-weight:bold; ; text-align:center; font-style:oblique",
           column(1, HTML("&nbsp;")),
           column(2, "preferably smaller than NNTlower"),
           column(2, "below this 'treat' is acceptable"),
           column(2, "Pr(BestToAct) = 1/NNT"),
           column(2, "above this 'treat' is unacceptable"),
           column(2, "preferably larger than NNTupper")
  ),
  plotOutput(outputId = "plotNNTgoals",
             height='250px'),
  hr(),

  #   fluidRow(
  #     "Positive predictive value = 1/NNTpos = ",
  #     textOutput("PPVderived"),
  #     br(),
  #     "Negative predictive value = 1 - 1/NNTneg = ",
  #     textOutput("NPVderived")
  #   ),
  #   sectionHeader(6,
  #                 numericInput("samplesize", label = "Prospective study sample size",
  #                              value=30, min = 10, max=1000, step = 1)
  #   ),
  #   sectionHeader(7, div(
  # h3("Required sensitivity and specificity for a retrospective study."),
  div(style="background:lightgrey",
      fluidRow(
        column(6,
               column(12, offset=2,
                      h3("Contra-Bayes plot", style="center"),
                      h4("Hover mouse on plot for calculations; click to set."),
                      tableOutput("parameterTable")
               ),
               column(10, offset=1,
                      plotOutput("contraBayesPlot",
                                 click="contraBayesPlot_click",
                                 hover="contraBayesPlot_hover",
                                 width='100%')
               )
        ),
        column(6,
               h3("Anticipated results:"),
               div(class="well container-fluid",
                   h3(style="font-style:oblique; text-indent:50px", "Prospective study"),
                   fluidRow(
                     column(2, h4(style=
                                  "text-indent:60px; position: relative; top: 50%;
                                  transform: translateY(-50%)",
                                  "Data")),
                     column(10,
                            fluidRow(
                              column(2, ""),
                              column(3, numericInput("Npositives",
                                                     label = " #positive",
                                                     value=30)),
                              column(4,  numericInput("NtruePositives",
                                                      label = " #TRUE positive",
                                                      value=15))
                            ),
                            fluidRow(
                              column(2, ""),
                              column(3, numericInput("Nnegatives",
                                                     label = " #negative",
                                                     value=30)),
                              column(4, numericInput("NtrueNegatives",
                                                     label = " #TRUE negative",
                                                     value=29))
                            )
                     )
                   ),
                   h4(style="text-indent:60px", "Predictive intervals (95%)"),
                   fluidRow(column(9, offset=3,
                                   tableOutput("intervalsProspective")))
               ),
               br(),


               div(class="well container-fluid",
                   h3(style="font-style:oblique; text-indent:50px", "Retrospective study"),
                   fluidRow(
                     column(2, h4(style=
                                    "text-indent:60px; position: relative; top: 50%;
                              transform: translateY(-50%); vertical-align:center",
                                  "Data")),
                     column(10,
                            fluidRow(
                              column(2, ""),
                              column(3, numericInput("Ncases",
                                                     label = "#cases",
                                                     value=30)),
                              column(4, numericInput("NposCases",
                                                     label = "#positive cases",
                                                     value=15))
                            ),
                            fluidRow(
                              column(2, ""),
                              column(3, numericInput("Ncontrols",
                                                     label = "#controls",
                                                     value=30)),
                              column(4, numericInput("NnegControls",
                                                     label = "#negative controls",
                                                     value=28))
                            )
                     )
                   ),
                   h4(style="text-indent:60px", "Predictive intervals (95%)"),
                   fluidRow(column(9, offset=3,
                                   tableOutput("intervalsRetrospective")))
               )
        )
      )
  )
))
