# ui script for PCA2
shinyUI(fluidPage(
  titlePanel(div(paste("PCA on the ",nomData," dataset"),style="color:#6E6E6E",align="center")),
  
  sidebarLayout(
    sidebarPanel(
      tags$head(
        tags$style("body {background-color: #EFFBFB; }"),
        tags$style(type='text/css', "#title1 { height: 25px; }"),
        tags$style(type='text/css', "#title2 { height: 25px; }")
        #tags$style(type='text/css', "#nb1 { height: 20px; }")
      ),
      wellPanel(
      div(align="center",checkboxInput("pcaparam","Show PCA parameters",FALSE)),
      conditionalPanel(
        condition="input.pcaparam==true",
        if(is.null(quantisup)){
        radioButtons("selecactive",label=h6("Choose the active variables"),
                          choices=list("All"="Toutes","Choose"="choix"),selected="Toutes")
        }
        else{
          radioButtons("selecactive",label=h6("Choose the active variables"),
                       choices=list("All"="Toutes","Choose"="choix"),selected="choix")
        },
        conditionalPanel(
          condition="input.selecactive=='choix'",
          selectInput("supvar",label="Select the supplementary quantitative variables",
                    choices=list(IdChoices=VariableChoices),
                    selected=quantisup,multiple=TRUE)
          ),
        br(),      
        h6("Select the supplementary categorical variables"),
       
          if(length(QualiChoice)>1){
            if(is.null(qualisup)){
            selectInput("supquali",label="",choices=list(Idqualisup=as.vector(QualiChoice)),multiple=TRUE)
            }
            else{
            selectInput("supquali",label="",choices=list(Idqualisup=as.vector(QualiChoice)),multiple=TRUE,selected=qualisup)  
            }
          }
          else if (length(QualiChoice)==1){
            if(is.null(qualisup)){
              checkboxInput("supquali",QualiChoice,FALSE)
            }
            else{
              checkboxInput("supquali",QualiChoice,TRUE)
            }
          }
          else if(length(QualiChoice)==0){
            p("No categorical variable in your dataset")
          },
        
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
        div(align="center",radioButtons("ind_var","",
                     choices=list("Graph of Individuals  "="Ind","Graph of Variables"="Var"),selected="var",inline=TRUE)),
        br(),
        conditionalPanel(
          condition="input.ind_var=='Ind'",
          textInput("title1",h6("Title of the graph : "), titre1),
          sliderInput("cex",h6("Size of labels"),min=0.5,max=2.5,value=size,step=0.05,ticks=FALSE),
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
            if(selection=="contrib"){
            div(align="center",sliderInput("slider0", label = "Number of the most contributive individuals",
                                           min = 1, max = length(nom), value =as.numeric(selection2),step=1))}
            else{
              div(align="center",sliderInput("slider0", label = "Number of the most contributive individuals",
                                             min = 1, max = length(nom), value =length(nom),step=1)) 
            }),
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
          condition="input.ind_var=='Var'",
          textInput("title2",h6("Title of the graph : "), titre2),
          sliderInput("cex2",h6("Size of labels"),min=0.5,max=2.5,value=size2,step=0.05,ticks=FALSE),
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
        )
      )
      ),
      wellPanel(
        h5("Save graphs as",align="center"),
        radioButtons("paramdown","",
                     choices=list("PNG"="png","JPG"="jpg","PDF"="pdf"),selected="png"),
        br(),
        div(align="center",actionButton("PCAcode", "Get the PCA code"))
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
                               p(downloadButton("downloadData2","Download as pdf"),align="center"))
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
                               verbatimTextOutput("summaryPCA"),
                               p(downloadButton("summary2","Download the summary"),align="center")
                               ),
                             conditionalPanel(
                               condition="input.out=='varsup'",
                               h6("Coordinates"),
                               div(align="center",tableOutput("sorties23")),
                               h6("Correlations"),
                               div(align="center",tableOutput("sorties32"))
                               ),
                             conditionalPanel(
                               condition="input.out=='qualico'",
                               h6("Coordinates"),
                               div(align="center",tableOutput("sorties12")),
                               h6("V-test"),
                               div(align="center",tableOutput("sorties13"))
                               ),
                             conditionalPanel(
                               condition="input.out=='supind'",
                               h6("Coordinates"),
                               div(align="center",tableOutput("sorties36")),
                               h6("Cos2"),
                               div(align="center",tableOutput("sorties37"))
                             )
                             ),
                    tabPanel("Automatic description of axes",
                             br(),
                             radioButtons("Dim",label="Choose the dimensions",choices=list("Dimension 1"="Dim1","Dimension 2"="Dim2","Dimension 3"="Dim3"),selected="Dim1"),
                             conditionalPanel(
                               condition="input.Dim=='Dim1'",
                               p("Quantitative"),
                               div(align="center",tableOutput("sortieDimdesc3")),
                               p("Qualitative"),
                               div(align="center",tableOutput("sortieDimdesc4"))
                               
                             ),
                             br(),
                             conditionalPanel(
                               condition="input.Dim=='Dim2'",
                               p("Quantitative"),
                               div(align="center",tableOutput("sortieDimdesc33")),
                               p("Qualitative"),
                               div(align="center",tableOutput("sortieDimdesc44"))
                             ),
                             br(),
                             conditionalPanel(
                               condition="input.Dim=='Dim3'",
                               p("Quantitative"),
                               div(align="center",tableOutput("sortieDimdesc333")),
                               p("Qualitative"),
                               div(align="center",tableOutput("sortieDimdesc444"))
                             )
                    ),
                    tabPanel("Summary of dataset",
                             br(),
                             verbatimTextOutput("summary"),
                             br(),
                             selectInput("bam",h6("Graphs for "),choices=list(IdChoices=VariableChoices),multiple=FALSE),
                             plotOutput("histo")),
                    
                    tabPanel("Data",
                             br(),
                             dataTableOutput("JDD")
                             )
        )
      )
    )
))
