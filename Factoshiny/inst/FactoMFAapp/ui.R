# ui2.R AFM

shinyUI(fluidPage(
  titlePanel(div(paste("MFA on the ",nomData," dataset"),style="color:#6E6E6E",align="center")),
  
  sidebarLayout(
      sidebarPanel(
        tags$head(
          tags$style("body {background-color: #E1D3FB; }"),
          tags$style(type='text/css', "#nameG1 { max-width: 100px; }"),
          tags$style(type='text/css', "#listvarG1 { max-width: 170px; }"),
          tags$style(type='text/css', "#nameG2 { max-width: 100px; }"),
          tags$style(type='text/css', "#listvarG2 { max-width: 170px; }"),
          tags$style(type='text/css', "#nameG3 { max-width: 100px; }"),
          tags$style(type='text/css', "#listvarG3 { max-width: 170px; }"),
          tags$style(type='text/css', "#nameG4 { max-width: 100px; }"),
          tags$style(type='text/css', "#listvarG4 { max-width: 170px; }"),
          tags$style(type='text/css', "#nameG5 { max-width: 100px; }"),
          tags$style(type='text/css', "#listvarG5 { max-width: 170px; }"),
          tags$style(type='text/css', "#nameG6 { max-width: 100px; }"),
          tags$style(type='text/css', "#listvarG6 { max-width: 170px; }"),
          tags$style(type='text/css', "#nameG7 { max-width: 100px; }"),
          tags$style(type='text/css', "#listvarG7 { max-width: 170px; }"),
          tags$style(type='text/css', "#nameG8 { max-width: 100px; }"),
          tags$style(type='text/css', "#listvarG8 { max-width: 170px; }"),
          tags$style(type='text/css', "#nameG9 { max-width: 100px; }"),
          tags$style(type='text/css', "#listvarG9 { max-width: 170px; }"),
          tags$style(type='text/css', "#nameG10 { max-width: 100px; }"),
          tags$style(type='text/css', "#nb1 { height:50px; }"),
          tags$style(type='text/css', "#nb2 { height:50px; }"),
          tags$style(type='text/css', "#listvarG10 { max-width: 170px; }"),
          tags$style(type='text/css', "#title1 { height: 25px; }"),
          tags$style(type='text/css', "#title2 { height: 25px; }"),
          tags$style(type='text/css', "#title3 { height: 25px; }"),
          tags$style(type='text/css', "#title4 { height: 25px; }"),
          tags$style(type='text/css', "#title5 { height: 25px; }")
        ),
        wellPanel(
        div(align="center",checkboxInput("graph","Show graphs options",FALSE)),
        conditionalPanel(
          condition="input.graph==true",
        div(align="center",selectInput("choixgraph",h6("Which graph would you like to modify ?"), choices=list("Groups"="group","Individuals"="ind","Quantitative variables"="quant","Frequences"="freq","Partial axes"="axes"),selected="group")),
        hr(),
        conditionalPanel(
          condition="input.choixgraph=='ind'",
          textInput("title2",h6("Title of the graph : "), title2),
          checkboxInput("meanind1","Draw points for the mean individuals",TRUE),
          checkboxInput("meanind","Draw labels for the mean individuals",TRUE),
          checkboxInput("qualind1","Draw points for the categories",TRUE),
          checkboxInput("qualind","Draw labels for the categories",TRUE),
          hr(),
          uiOutput("drawindiv"),
          conditionalPanel(
            condition="input.drawind=='c'",
            uiOutput("habillagequali")
            ),
          hr(),
          radioButtons("choixpartial","Partial points are drawn",choices=list("None"=1,"All"=2,"Choose"=3),selected=1,inline=TRUE),
          conditionalPanel(
            condition="input.choixpartial==3",
            selectInput("indivpartiel",label="Selectionnez les individus :",
                      choices=list(num=nom),multiple=TRUE)),
          hr(),
          conditionalPanel(
            condition="input.choixpartial!=1",
            checkboxInput("partind","Draw labels for the partial individuals",TRUE))
          ),
        conditionalPanel(
          condition="input.choixgraph=='quant'",
          textInput("title3",h6("Title of the graph : "), title3),
          radioButtons("selection","Select from",choices=list("No selection"="no","Contribution"="contrib","Cos2"="cos2"),selected="no"),
          uiOutput("slider1"),
          hr(),
          uiOutput("hide2"),
          checkboxInput("colorgroup","Color the variables by group",TRUE) 
        ),
        conditionalPanel(
          condition="input.choixgraph=='freq'",
          textInput("title5",h6("Title of the graph : "), title5),
          checkboxInput("affichind","Draw labels for the mean individuals",TRUE),
          checkboxInput("affichcol","Draw labels for the columns",TRUE)
        ),
        conditionalPanel(
          condition="input.choixgraph=='axes'",
          textInput("title4",h6("Title of the graph : "), title4),
          checkboxInput("coloraxe","Color the partial axe by group",TRUE)
        ),
        conditionalPanel(
          condition="input.choixgraph=='group'",
          textInput("title1",h6("Title of the graph : "), title1)
        ),
        fluidRow(
          column(5,selectInput("nb1", label = h6("  x axis"), 
                               choices = list("1" = 1, "2" = 2, "3" = 3,"4"= 4,"5" =5), selected = 1,width='80%')),
          column(5,selectInput("nb2", label =h6("   y axis"), 
                               choices = list("1" = 1, "2" = 2,"3" = 3,"4"= 4,"5" =5), selected = 2,width='80%')))
        )),
        wellPanel(
          h5("Save graphs as",align="center"),
          radioButtons("paramdown","",
                      choices=list("PNG"="png","JPG"="jpg","PDF"="pdf"),selected="png")
        ),
        div(align="center",actionButton("HCPCcode", "Get the MFA code")),
        br(),
        div(align="center",actionButton("Quit", "Quit the app"))
        ),
      
      mainPanel(
        tabsetPanel(id = "graph_sort",
                    tabPanel("Creation of groups",
                             br(),
                             checkboxInput("activemodif","Create the groups",FALSE),
                             conditionalPanel(
                               condition="input.activemodif==true",
                               h6("Group 1"),
                              fluidRow(
                              column(3,
                                     radioButtons("typeG1"," ",choices=list("Quantitative"="quant","Qualitative"="qual","Frequencies"="freq"),selected="quant")
                              ),
                              column(3,
                                     uiOutput("listvarG1")),
                              column(3,
                                     textInput("nameG1", label = " ", value = "Name 1"),
                                     conditionalPanel(
                                       condition="input.typeG1=='quant'",
                                       radioButtons("scale1","",choices=list("Scaled"="sc","Unscaled"="un"),selected="sc",inline=TRUE))),
                              column(2,
                                     radioButtons("typeG12","",choices=list("Active"="act","Supplementary"="sup"),selected="act"))),
                              br(),
                              checkboxInput("activeG2",h6("Create Group 2"),FALSE),
                              conditionalPanel(
                                condition="input.activeG2==true",
                                fluidRow(
                                column(3,
                                       radioButtons("typeG2"," ",choices=list("Quantitative"="quant","Qualitative"="qual","Frequencies"="freq"),selected="quant")
                                ),
                                column(3,
                                       uiOutput("listvarG2")),
                                column(3,
                                       textInput("nameG2", label = " ", value = "Name 2"),
                                       conditionalPanel(
                                         condition="input.typeG2=='quant'",
                                         radioButtons("scale2","",choices=list("Scaled"="sc","Unscaled"="un"),selected="sc",inline=TRUE))),
                                column(2,
                                       radioButtons("typeG22","",choices=list("Active"="act","Supplementary"="sup"),selected="act"))),
                                br(),
                                checkboxInput("activeG3",h6("Create Group 3"),FALSE),
                                conditionalPanel(
                                  condition="input.activeG3==true",
                                  fluidRow(
                                  column(3,
                                         radioButtons("typeG3"," ",choices=list("Quantitative"="quant","Qualitative"="qual","Frequencies"="freq"),selected="quant")
                                  ),
                                  column(3,
                                         uiOutput("listvarG3")),
                                  column(3,
                                         textInput("nameG3", label = " ", value = "Name 3"),
                                         conditionalPanel(
                                           condition="input.typeG3=='quant'",
                                           radioButtons("scale3","",choices=list("Scaled"="sc","Unscaled"="un"),selected="sc",inline=TRUE))),
                                  column(3,
                                         radioButtons("typeG32","",choices=list("Active"="act","Supplementary"="sup"),selected="act"))),
                                  br(),
                                  checkboxInput("activeG4",h6("Create Group 4"),FALSE),
                                  conditionalPanel(
                                    condition="input.activeG4==true",
                                    fluidRow(
                                    column(3,
                                           radioButtons("typeG4"," ",choices=list("Quantitative"="quant","Qualitative"="qual","Frequencies"="freq"),selected="quant")
                                    ),
                                    column(3,
                                           uiOutput("listvarG4")),
                                    column(3,
                                           textInput("nameG4", label = " ", value = "Name 4"),
                                           conditionalPanel(
                                             condition="input.typeG4=='quant'",
                                             br(),
                                             radioButtons("scale4","",choices=list("Scaled"="sc","Unscaled"="un"),selected="sc",inline=TRUE))),
                                    column(3,
                                           radioButtons("typeG42","",choices=list("Active"="act","Supplementary"="sup"),selected="act"))),
                                    br(),
                                    checkboxInput("activeG5",h6("Create Group 5"),FALSE),
                                    conditionalPanel(
                                      condition="input.activeG5==true",
                                      fluidRow(
                                      column(3,
                                             radioButtons("typeG5"," ",choices=list("Quantitative"="quant","Qualitative"="qual","Frequencies"="freq"),selected="quant")
                                      ),
                                      column(3,
                                             uiOutput("listvarG5")),
                                      column(3,
                                             textInput("nameG5", label = " ", value = "Name 5"),
                                             conditionalPanel(
                                               condition="input.typeG5=='quant'",
                                               br(),
                                               radioButtons("scale5","",choices=list("Scaled"="sc","Unscaled"="un"),selected="sc",inline=TRUE))),
                                      column(3,
                                             radioButtons("typeG52","",choices=list("Active"="act","Supplementary"="sup"),selected="act"))),
                                      br(),
                                      checkboxInput("activeG6",h6("Create Group 6"),FALSE),
                                      conditionalPanel(
                                        condition="input.activeG6==true",
                                        fluidRow(
                                        column(3,
                                               radioButtons("typeG6"," ",choices=list("Quantitative"="quant","Qualitative"="qual","Frequencies"="freq"),selected="quant")
                                        ),
                                        column(3,
                                               uiOutput("listvarG6")),
                                        column(3,
                                               textInput("nameG6", label = " ", value = "Name 6"),
                                               conditionalPanel(
                                                 condition="input.typeG6=='quant'",
                                                 br(),
                                                 radioButtons("scale6","",choices=list("Scaled"="sc","Unscaled"="un"),selected="sc",inline=TRUE))),
                                        column(3,
                                               radioButtons("typeG62","",choices=list("Active"="act","Supplementary"="sup"),selected="act"))),
                                        br(),
                                        checkboxInput("activeG7",h6("Create Group 7"),FALSE),
                                        conditionalPanel(
                                          condition="input.activeG7==true",
                                          fluidRow(
                                          column(3,
                                                 radioButtons("typeG7"," ",choices=list("Quantitative"="quant","Qualitative"="qual","Frequencies"="freq"),selected="quant")
                                          ),
                                          column(3,
                                                 uiOutput("listvarG7")),
                                          column(3,
                                                 textInput("nameG7", label = " ", value = "Name 7"),
                                                 conditionalPanel(
                                                   condition="input.typeG7=='quant'",
                                                   br(),
                                                   radioButtons("scale7","",choices=list("Scaled"="sc","Unscaled"="un"),selected="sc",inline=TRUE))),
                                          column(3,
                                                 radioButtons("typeG72","",choices=list("Active"="act","Supplementary"="sup"),selected="act"))),
                                          br(),
                                          checkboxInput("activeG8",h6("Create Group 8"),FALSE),
                                          conditionalPanel(
                                            condition="input.activeG8==true",
                                            fluidRow(
                                            column(3,
                                                   radioButtons("typeG8"," ",choices=list("Quantitative"="quant","Qualitative"="qual","Frequencies"="freq"),selected="quant")
                                            ),
                                            column(3,
                                                   uiOutput("listvarG8")),
                                            column(3,
                                                   textInput("nameG8", label = " ", value = "Name 8"),
                                                   conditionalPanel(
                                                     condition="input.typeG8=='quant'",
                                                     br(),
                                                     radioButtons("scale8","",choices=list("Scaled"="sc","Unscaled"="un"),selected="sc",inline=TRUE))),
                                            column(3,
                                                   radioButtons("typeG82","",choices=list("Active"="act","Supplementary"="sup"),selected="act"))),
                                            br(),
                                            checkboxInput("activeG9",h6("Create Group 9"),FALSE),
                                            conditionalPanel(
                                              condition="input.activeG9==true",
                                              fluidRow(
                                              column(3,
                                                     radioButtons("typeG9"," ",choices=list("Quantitative"="quant","Qualitative"="qual","Frequencies"="freq"),selected="quant")
                                              ),
                                              column(3,
                                                     uiOutput("listvarG9")),
                                              column(3,
                                                     textInput("nameG9", label = " ", value = "Name 9"),
                                                     conditionalPanel(
                                                       condition="input.typeG9=='quant'",
                                                       br(),
                                                       radioButtons("scale9","",choices=list("Scaled"="sc","Unscaled"="un"),selected="sc",inline=TRUE))),
                                              column(3,
                                                     radioButtons("typeG92","",choices=list("Active"="act","Supplementary"="sup"),selected="act"))),
                                              br(),
                                              checkboxInput("activeG10",h6("Create Group 10"),FALSE),
                                              conditionalPanel(
                                                condition="input.activeG10==true",
                                                fluidRow(
                                                column(3,
                                                       radioButtons("typeG10"," ",choices=list("Quantitative"="quant","Qualitative"="qual","Frequencies"="freq"),selected="quant")
                                                ),
                                                column(3,
                                                       uiOutput("listvarG10")),
                                                column(3,
                                                       textInput("nameG10", label = " ", value = "Name 10"),
                                                       conditionalPanel(
                                                         condition="input.typeG10=='quant'",
                                                         br(),
                                                         radioButtons("scale10","",choices=list("Scaled"="sc","Unscaled"="un"),selected="sc",inline=TRUE))),
                                                column(3,
                                                       radioButtons("typeG102","",choices=list("Active"="act","Supplementary"="sup"),selected="act")))
                                              )
                                            )
                                          )
                                        )
                                      )
                                    )
                                  )
                                  )
                                )
                             )),
                    tabPanel("Graphs",
                             br(),
                             div(align = "center",plotOutput("map5", width = 500, height=500)),
                             br(),
                             conditionalPanel(
                               condition="input.paramdown=='png'",
                               p(downloadButton("downloadData11","Download as png"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='jpg'",
                               p(downloadButton("downloadData12","Download as jpg"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='pdf'",
                               p(downloadButton("downloadData13","Download as pdf"),align="center")),
                             div(align = "center",plotOutput("map", width = 500, height=500)),
                             br(),
                             conditionalPanel(
                               condition="input.paramdown=='jpg'",
                               p(downloadButton("downloadData1","Download as jpg"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='png'",
                               p(downloadButton("downloadData","Download as png"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='pdf'",
                               p(downloadButton("downloadData2","Download as pdf"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='emf'",
                               p(downloadButton("downloadData7","Download as emf"),align="center")),
                             br(),
                             div(align = "center",uiOutput("map22")),
                             br(),
                             conditionalPanel(
                               condition="input.paramdown=='jpg'",
                               div(uiOutput("download4"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='png'",
                               div(uiOutput("download3"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='pdf'",
                               div(uiOutput("download5"),align="center")),
                             br(),
                             div(align = "center",plotOutput("map4", width = 500, height=500)),
                             br(),
                             conditionalPanel(
                               condition="input.paramdown=='png'",
                               p(downloadButton("downloadData15","Download as png"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='jpg'",
                               p(downloadButton("downloadData16","Download as jpg"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='pdf'",
                               p(downloadButton("downloadData17","Download as pdf"),align="center")),
                             div(align = "center",uiOutput("map66")),
                             br(),
                             conditionalPanel(
                               condition="input.paramdown=='png'",
                               div(uiOutput("download19"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='jpg'",
                               div(uiOutput("download20"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='pdf'",
                               div(uiOutput("download21"),align="center"))
                             ),

                    tabPanel("Values",
                             br(),
                             radioButtons("out","Which value would you like to see ?",
                                          choices=list("Summary of MFA"="MFA","Eigenvalues"="eig","Results for the individuals"="ind",
                                                       "Results for the quantitative variables"="quantvar","Results for the groups"="group","Results for the partial axes"="partaxe"),selected="eig",inline=TRUE),
                             conditionalPanel(
                               condition="input.out=='eig'",
                               tableOutput("sorties"),
                               plotOutput("map3")),
                             conditionalPanel(
                               condition="input.out=='MFA'",
                               verbatimTextOutput("summaryMFA")
                               ),
                             conditionalPanel(
                               condition="input.out=='ind'",
                               radioButtons("out2","What type of results ?",choices=list("Coordinates"="coord","Contributions"="contrib","CosÂ²"="cos2","Within inertia"="witi",
                                                                                         "Partial Coordinates"="partco","Within partial inertia"="wpi"),selected="coord",inline=TRUE),
                               conditionalPanel(
                                 condition="input.out2=='coord'",
                                 div(align="center",tableOutput("sorties1"))),
                               conditionalPanel(
                                 condition="input.out2=='contrib'",
                                 div(align="center",tableOutput("sorties2"))),
                               conditionalPanel(
                                 condition="input.out2=='cos2'",
                                 div(align="center",tableOutput("sorties3"))),
                               conditionalPanel(
                                 condition="input.out2=='witi'",
                                 div(align="center",tableOutput("sorties4"))),
                               conditionalPanel(
                                 condition="input.out2=='partco'",
                                 div(align="center",tableOutput("sorties5"))),
                               conditionalPanel(
                                 condition="input.out2=='wpi'",
                                 div(align="center",tableOutput("sorties6")))
                               ),
                             conditionalPanel(
                               condition="input.out=='quantvar'",
                               radioButtons("out3","What type of results ?",choices=list("Coordinates"="coord","Contributions"="contrib","CosÂ²"="cos2","Correlations"="cor"),selected="coord",inline=TRUE),
                               conditionalPanel(
                                 condition="input.out3=='coord'",
                                 div(align="center",tableOutput("sorties11"))),
                               conditionalPanel(
                                 condition="input.out3=='contrib'",
                                 div(align="center",tableOutput("sorties22"))),
                               conditionalPanel(
                                 condition="input.out3=='cos2'",
                                 div(align="center",tableOutput("sorties33"))),
                               conditionalPanel(
                                 condition="input.out3=='cor'",
                                 div(align="center",tableOutput("sorties44")))
                             ),
                             conditionalPanel(
                               condition="input.out=='group'",
                               div(align="center",tableOutput("sortiegroup"))
                               ),
                             conditionalPanel(
                               condition="input.out=='partaxe'",
                               radioButtons("out4","What type of results ?",choices=list("Coordinates"="coord","Correlations"="cor","Contribution"="contrib","Correlations between"="cor.btw"),selected="coord",inline=TRUE),
                               conditionalPanel(
                                 condition="input.out4=='coord'",
                                 div(align="center",tableOutput("sorties12"))),
                               conditionalPanel(
                                 condition="input.out4=='cor'",
                                 div(align="center",tableOutput("sorties23"))),
                               conditionalPanel(
                                 condition="input.out4=='contrib'",
                                 div(align="center",tableOutput("sorties34"))),
                               conditionalPanel(
                                 condition="input.out4=='cor.btw'",
                                 div(align="center",tableOutput("sorties45")))
                               )
                             
                             ),
                    tabPanel("Summary of dataset",
                             br(),
                             verbatimTextOutput("summary"),
                             selectInput("bam","Graphs for",choices=list(IdChoices=VariableChoices),multiple=FALSE),
                             plotOutput("histo")),
                    
                    tabPanel("Data",
                             br(),
                             dataTableOutput("JDD")
                             )
        )
      )
    )
))
