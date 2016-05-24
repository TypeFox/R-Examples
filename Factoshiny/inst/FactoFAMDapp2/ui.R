# ui script for FAMD2
shinyUI(fluidPage(
  titlePanel(div(paste("FAMD on the ",nomData," dataset"),style="color:#6E6E6E",align="center")),
  
  sidebarLayout(
    sidebarPanel(
      tags$head(
        tags$style("body {background-color: #FFF0C7; }"),
        tags$style(type='text/css', "#title1 { height: 25px; }"),
        tags$style(type='text/css', "#title2 { height: 25px; }"),
        tags$style(type='text/css', "#title3 { height: 25px; }")
      ),
      wellPanel(
      div(align="center",checkboxInput("pcaparam","Show FAMD parameters",FALSE)),
      conditionalPanel(
        condition="input.pcaparam==true",
        if(is.null(quantisup) && is.null(qualisup)){
        radioButtons("selecactive",label=h6("Choose the active variables"),
                          choices=list("All"="Toutes","Choose"="choix"),selected="Toutes")
        }
        else{
          radioButtons("selecactive",label=h6("Choose the active variables"),
                       choices=list("All"="Toutes","Choose"="choix"),selected="choix")
        },
        conditionalPanel(
          condition="input.selecactive=='choix'",
          h6("Select the supplementary quantitative variables"),
          if(length(VariableChoices)>1){
            selectInput("supvar",label="",
                    choices=list(IdChoices=VariableChoices),
                    selected=quantisup,multiple=TRUE)}
          else{
            selectInput("supvar",label="",
                        choices=VariableChoices,
                        selected=quantisup,multiple=TRUE)
          },
          h6("Select the supplementary qualitative variables"),
          if(length(QualiChoice)>1){
          selectInput("supvar1",label="",
                      choices=list(Idqualisup=QualiChoice),
                      selected=qualisup,
                      multiple=TRUE)}
          else{
          selectInput("supvar1",label="",
                        choices=quali,selected=qualisup,
                        multiple=TRUE)  
          }
          ),
        br(),
        h6("Select the supplementary individuals"),
        if(is.null(indsupl)){
          selectInput("indsup","",choices=list(num=nom), multiple=TRUE)
        }
        else{
          selectInput("indsup","",choices=list(num=nom), multiple=TRUE,selected=indsupl)
        }
      )
      ),
      wellPanel(
      div(align="center",checkboxInput("graph","Show graphs options",FALSE)),
      conditionalPanel(
        condition="input.graph==true",
        fluidRow(
          column(5,uiOutput("NB1")),
          column(5,uiOutput("NB2"))),
        hr(),
        div(align="center",selectInput("choixgraph",h6("Which graph would you like to modify ?"), choices=list("Individuals and categorical variables"="ind","Variables"="var","Quantitative variables"="quant"),selected="ind")),
        br(),
        conditionalPanel(
          condition="input.choixgraph=='ind'",
          textInput("title1",h6("Title of the graph : "), title1),
          sliderInput("cex",h6("Size of labels"),min=0.5,max=2.5,value=size,step=0.05,ticks=FALSE),
          br(),
          checkboxInput("labels2","Draw labels of individuals",labind),
          checkboxInput("labels","Draw labels of variables",labvar),
          selectInput("select",label=h6("Draw individuals according to :"),
                      choices=list("No selection"="NONE","cos2"="cos2","Contribution"="contrib","Manual"="Manuel"),selected=selection),
          conditionalPanel(
            condition="input.select=='cos2'",
            if(selection=="cos2"){
            div(align="center",sliderInput("slider1", label = "cos2",
                        min = 0, max = 1, value =as.numeric(selection2),step=0.05))}
            else{
              div(align="center",sliderInput("slider1", label = "cos2",
                                             min = 0, max = 1, value =0,step=0.05))
            }),
          conditionalPanel(
            condition="input.select=='contrib'",
            uiOutput("slider7")),
          conditionalPanel(
            condition="input.select=='Manuel'",
            if(selection=="Manuel"){
            selectInput("indiv",label="Select individuals :",
                        choices=list(num=nom),multiple=TRUE,selected=selection2)}
            else{
              selectInput("indiv",label="Selectionnez les individus :",
                          choices=list(num=nom),multiple=TRUE)
            }),
            if(is.null(habillageind)){
              checkboxInput("habi","Points colour depend on categorical variable",FALSE)
            }
            else{
              checkboxInput("habi","Points colour depend on categorical variable",TRUE)
            },
            conditionalPanel(
              condition="input.habi==true",
              uiOutput("habillage2")
            )
        ),
        conditionalPanel(
          condition="input.choixgraph=='var'",
          textInput("title2",h6("Title of the graph : "), title2),
          sliderInput("cex2",h6("Size of labels"),min=0.5,max=2.5,value=size2,step=0.05,ticks=FALSE),
          br(),
          selectInput("select0",label=h6("Draw variables according to :"),
                      choices=list("No selection"="NONE","cos2"="cos2","Contribution"="contrib"),selected=selection3),
          conditionalPanel(
            condition="input.select0=='contrib'",
            uiOutput("slider3")
            ),
          conditionalPanel(
            condition="input.select0=='cos2'",
            if(selection3=="cos2"){
              div(align="center",sliderInput("slider00", label = "cos2",
                                             min = 0, max = 1, value =as.numeric(selection4),step=0.05))  
            }
            else{
            div(align="center",sliderInput("slider00", label = "cos2",
                                           min = 0, max = 1, value =0,step=0.05))}
          )
        ),
        conditionalPanel(
          condition="input.choixgraph=='quant'",
          textInput("title3",h6("Title of the graph : "), title3),
          sliderInput("cex3",h6("Size of labels"),min=0.5,max=2.5,value=size3,step=0.05,ticks=FALSE),
          br(),
          selectInput("selecti",label=h6("Draw variables according to :"),
                      choices=list("No selection"="NONE","cos2"="cos2","Contribution"="contrib"),selected=selection5),
          conditionalPanel(
            condition="input.selecti=='contrib'",
            uiOutput("slider5")
          ),
          conditionalPanel(
            condition="input.selecti=='cos2'",
            if(selection3=="cos2"){
              div(align="center",sliderInput("slider000", label = "cos2",
                                             min = 0, max = 1, value =as.numeric(selection6),step=0.05))  
            }
            else{
              div(align="center",sliderInput("slider000", label = "cos2",
                                             min = 0, max = 1, value =0,step=0.05))}
          )
        )
      )
      ),
      wellPanel(
        h5("Save graphs as",align="center"),
        radioButtons("paramdown","",
                     choices=list("PNG"="png","JPG"="jpg","PDF"="pdf"),selected="png"),
        br(),
        div(align="center",actionButton("FAMDcode", "Get the FAMD code"))
        ),
      div(align="center",actionButton("Quit", "Quit the app"))
      ),
      
      mainPanel(
        tabsetPanel(id = "graph_sort",
                    tabPanel("Graphs",
                             br(),
                             div(align = "center",plotOutput("map2", width = 500, height=500)),
                             br(),
                             conditionalPanel(
                               condition="input.paramdown=='jpg'",
                               p(downloadButton("downloadData4","Download as jpg"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='png'",
                               p(downloadButton("downloadData3","Download as png"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='pdf'",
                               p(downloadButton("downloadData5","Download as pdf"),align="center")),
                             br(),
                             div(align="center",plotOutput("map", width = 500, height=500)),
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
                             br(),
                             div(align="center",plotOutput("map4", width = 500, height=500)),
                             br(),
                             conditionalPanel(
                               condition="input.paramdown=='jpg'",
                               p(downloadButton("downloadData6","Download as jpg"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='png'",
                               p(downloadButton("downloadData7","Download as png"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='pdf'",
                               p(downloadButton("downloadData8","Download as pdf"),align="center"))
                             ),

                    tabPanel("Values",
                             br(),
                             uiOutput("out22"),
                             br(),
                             conditionalPanel(
                               condition="input.out=='eig'",
                               div(align="center",tableOutput("sorties")),
                               plotOutput("map3")
                               ),
                             conditionalPanel(
                               condition="input.out=='resvar'",
                               h6("Coordinates"),
                               div(align="center",tableOutput("sorties2")),
                               br(),
                               h6("Contributions"),
                               div(align="center",tableOutput("sorties3")),
                               br(),
                               h6("Cos2"),
                               div(align="center",tableOutput("sorties4"))),
                             conditionalPanel(
                               condition="input.out=='resind'",
                               h6("Coordinates"),
                               div(align="center",tableOutput("sorties22")),
                               br(),
                               h6("Contributions"),
                               div(align="center",tableOutput("sorties33")),
                               br(),
                               h6("Cos2"),
                               div(align="center",tableOutput("sorties44"))),
                             conditionalPanel(
                               condition="input.out=='ACP'",
                               numericInput("nbele","Number of elements to print",value=10),
                               verbatimTextOutput("summaryFAMD"),
                               p(downloadButton("summary2","Download the summary"),align="center")
                               ),
                             conditionalPanel(
                               condition="input.out=='varsup'",
                               h6("Coordinates"),
                               div(align="center",tableOutput("sorties23")),
                               h6("Cos2"),
                               div(align="center",tableOutput("sorties32"))
                               ),
                  
                             conditionalPanel(
                               condition="input.out=='supind'",
                               h6("Coordinates"),
                               div(align="center",tableOutput("sorties36")),
                               h6("Cos2"),
                               div(align="center",tableOutput("sorties37"))
                             )
                             ),
                    tabPanel("Summary of dataset",
                             br(),
                             verbatimTextOutput("summary"),
                             br(),
                             selectInput("bam",h6("Graphs for "),choices=list(Idall=all),multiple=FALSE),
                             plotOutput("histo")),
                    
                    tabPanel("Data",
                             br(),
                             dataTableOutput("JDD")
                             )
        )
      )
    )
))
