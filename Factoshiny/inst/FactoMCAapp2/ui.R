# ui script for MCA2

shinyUI(fluidPage(
  titlePanel(div(paste("MCA on the ",nomData," dataset"),style="color:#6E6E6E",align="center")),
  
  sidebarLayout(
    sidebarPanel(
      tags$head(
        tags$style("body {background-color: #FBEFEF; }"),
        tags$style(type='text/css', "#title1 { height: 25px; }"),
        tags$style(type='text/css', "#title2 { height: 25px; }"),
        tags$style(type='text/css', "#title3 { height: 25px; }")
      ),
      wellPanel(
        div(align="center",checkboxInput("mcaparam","Show MCA parameters",FALSE)),
        conditionalPanel(
          condition="input.mcaparam==true",
          if(is.null(supquali)){
            radioButtons("selecactive",label=h6("Choose the active qualitative variables"),
                         choices=list("All"="Toutes","Choose"="choix"),selected="Toutes")
          }
          else{
            radioButtons("selecactive",label=h6("Choose the active qualitative variables"),
                         choices=list("All"="Toutes","Choose"="choix"),selected="choix")
          },
          conditionalPanel(
            condition="input.selecactive=='choix'",
            selectInput("supvar",label=h6("Select the supplementary qualitative variables"),
                        choices=list(IdChoices=VariableChoices),
                        selected=supquali,multiple=TRUE)
          ),
          
          #Selection des variables quantitatives supplementaires
          h6("Select the supplementary quantitative variables"),
          if(length(QuantiChoice)>1){
            if(is.null(quantiS)){
            selectInput("supquanti",label="",choices=list(Idquantisup=as.vector(QuantiChoice)),multiple=TRUE)
          }
          else {
            selectInput("supquanti",label="",choices=list(Idquantisup=as.vector(QuantiChoice)),multiple=TRUE,selected=quantiS)
            
          }
          }
          else if (length(QuantiChoice)==1){
            if(is.null(quantiS)){
            checkboxInput("supquanti",QuantiChoice,FALSE)
          }
            else{
            checkboxInput("supquanti",QuantiChoice,TRUE)  
            }
          }
          else if(length(QuantiChoice)==0){
            p("No quantitative variable in your dataset")
          },
          
          
          h6("Supplementary individuals"),
          if(is.null(indsupl)){
            selectInput("indsup","",choices=list(num=nom), multiple=TRUE)
          }
          else{
            selectInput("indsup","",choices=list(num=nom), multiple=TRUE,selected=indsupl)
          }
        )
      ),
      
      #Prametres graphiques
      wellPanel(
        div(align="center",checkboxInput("graph","Show graphs options",FALSE)),
        conditionalPanel(
          condition="input.graph==true",
          
          #
          fluidRow(
            column(5,uiOutput("NB1")),
            column(5,uiOutput("NB2"))),
          hr(),
          uiOutput("choixchange"),
          hr(),
          conditionalPanel(
            condition="input.MCAgraph=='ind'",
            p("Graph of individuals and modalities",align="center"),
            uiOutput("choixindvar"),
            br(),
           # p("Draw labels for :",align="center"),
          #  uiOutput("pointlabel"),
            textInput("title1",h6("Title of the graph : "), title1),
            div(align="center",radioButtons("modind",h6("Select elements to modify",align="center"),choices=list("Individuals"="Ind","Modalities"="Mod"),selected="Mod",inline=TRUE)),
            br(),
            conditionalPanel(
            condition="input.modind=='Ind'",
            if(selection=="NONE"){
              selectInput("select",label=h6("Select individuals from "),
                          choices=list("No selection"="NONE","cos2"="cos2","Manual"="Manuel","Contribution"="Contrib"),selected="NONE")
            }
            else{
              selectInput("select",label=h6("Select individuals from "),
                          choices=list("No selection"="NONE","cos2"="cos2","Manual"="Manuel","Contribution"="Contrib"),selected=selection)
              },
             conditionalPanel(
              condition="input.select=='cos2'",
              if(selection=="cos2"){
                div(align="center",sliderInput("slider1", label = "cos2",
                                               min = 0, max = 1, value =selection2,step=0.05))
              }
              else{
                div(align="center",sliderInput("slider1", label = "cos2",
                                               min = 0, max = 1, value =0,step=0.05))
              }),
              
            conditionalPanel(
              condition="input.select=='Manuel'",
              if(selection=="Manuel"){
                selectInput("indiv",label="Selectionnez les individus :",
                            choices=list(num=nom),multiple=TRUE,selected=selection2) 
              }
              else{
                selectInput("indiv",label="Selectionnez les individus :",
                            choices=list(num=nom),multiple=TRUE)  
              }),
            conditionalPanel(
              condition="input.select=='Contrib'",
              if(selection=="Contrib"){
                sliderInput("sliderContrib",label="Nombre d'individus les plus contributifs",
                            min=1,max=length(nom),value=selection2,step=1)  
              }
              else{
                sliderInput("sliderContrib",label="Nombre d'individus les plus contributifs",
                            min=1,max=length(nom),value=length(nom),step=1)
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
            condition="input.modind=='Mod'",
            if(selection3=="NONE"){
              selectInput("selectMod",label=h6("Select modalities from "),
                          choices=list("No selection"="NONE","cos2"="cos2","Contribution"="Contrib"),selected="NONE") 
            }
            else{
              selectInput("selectMod",label=h6("Select modalities from "),
                          choices=list("No selection"="NONE","cos2"="cos2","Contribution"="Contrib"),selected=selection3)
              },
            conditionalPanel(
              condition="input.selectMod=='cos2'",
              if(selection3=="cos2"){
                div(align="center",sliderInput("sliderCosMod", label = "cos2",
                                               min = 0, max = 1, value =as.numeric(selection4),step=0.05))}  
              else{
                div(align="center",sliderInput("sliderCosMod", label = "cos2",
                                               min = 0, max = 1, value=0,step=0.05))  
              }),
            conditionalPanel(
              condition="input.selectMod=='Contrib'",
              uiOutput("slider3")
            )
          )),
          conditionalPanel(
            condition="input.MCAgraph=='var'",
            p("Graph of variables",align="center"),
            textInput("title2",h6("Title of the graph : "), title2),
          div(align="center",checkboxGroupInput("var_sup",h6(""),choices=list("Supplementary qualitative variables"="suplquali","Supplementary quantitative variables"="suplquanti","Active qualitative variables"="act"),selected=varsup))
        ),
        conditionalPanel(
          condition="input.MCAgraph=='quant'",
          p("Graph of the supplementary quantitative variables",align="center"),
          textInput("title3",h6("Title of the graph : "), title3)))
      ),
      
      
      wellPanel(
        h5("Save graphs as",align="center"),
        radioButtons("paramdown","",
                     choices=list("PNG"="png","JPG"="jpg","PDF"="pdf"),selected="png"),
      br(),
      div(align="center",actionButton("MCAcode", "Get the MCA code"))),

    br(),
    div(align="center",actionButton("Quit", "Quit the app"))
    ),
    mainPanel(
      tags$style(type = "text/css", "a{color: #B53977;}"),
      tabsetPanel(id = "graph_sort",
                  tabPanel("Graphs",
                           br(),
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
                           br(),
                           div(align="center",plotOutput("map4", width = 500, height=500)),
                           br(),
                           conditionalPanel(
                             condition="input.paramdown=='jpg'",
                             p(downloadButton("downloadData10","Download as jpg"),align="center")),
                           conditionalPanel(
                             condition="input.paramdown=='png'",
                             p(downloadButton("downloadData0","Download as png"),align="center")),
                           conditionalPanel(
                             condition="input.paramdown=='pdf'",
                             p(downloadButton("downloadData20","Download as pdf"),align="center")),
                           
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
                             div(uiOutput("download5"),align="center"))
                  ),
                  
                  ####
                  
                  tabPanel("Values",
                           br(),
                           uiOutput("out22"),
                           br(), 
                           conditionalPanel(
                             condition="input.out=='eig'",
                             div(align="center",tableOutput("sorties")),
                             div(align="center",plotOutput("map3",width = 500, height=300))),
                           
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
                           #Mettre un if pour le nombre d'individus
                           conditionalPanel(
                             condition="input.out=='resind'",
                             h6("Coordinates"),
                             div(align="center",dataTableOutput("sorties22")),
                             br(),
                             h6("Contributions"),
                             div(align="center",dataTableOutput("sorties33")),
                             br(),
                             h6("Cos2"),
                             div(align="center",dataTableOutput("sorties44"))),
                           
                           conditionalPanel(
                             condition="input.out=='MCA'",
                             numericInput("nbele",h6("Number of elements to print"),value=10),
                             br(),
                             verbatimTextOutput("summaryMCA"),
                             p(downloadButton("summary2","Download the summary"),align="center")
                           ),
                           conditionalPanel(
                             condition="input.out=='varsup'",
                             h6("Coordinates"),
                             div(align="center",tableOutput("sorties23")),
                             br(),
                             h6("cos2"),
                             div(align="center",tableOutput("sorties232")),
                             br(),
                             h6("v.test"),
                             div(align="center",tableOutput("sorties233"))),
                           
                           conditionalPanel(
                             condition="input.out=='quantico'",
                             h6("Coordinates"),
                             div(align="center",tableOutput("sorties43"))),
                           
                           conditionalPanel(
                             condition="input.out=='Isup'",
                             h6("Coordinates"),
                             div(align="center",tableOutput("sortiesIsupC")),
                             br(),
                             h6("Cos2"),
                             div(align="center",tableOutput("sortiesIsupCos")))
                  ),
                  
                  tabPanel("Automatic description of axes",
                           br(),
                           radioButtons("Dim",label="Choose the dimensions",choices=list("Dimension 1"="Dim1","Dimension 2"="Dim2","Dimension 3"="Dim3"),selected="Dim1"),
                           conditionalPanel(
                             condition="input.Dim=='Dim1'",
                             p("variables"),
                             div(align="center",tableOutput("sortieDimdesc2")),
                             br(),
                             p("categorical"),
                             div(align="center",tableOutput("sortieDimdesc")),
                             br(),
                             p("quantitative"),
                             div(align="center",tableOutput("sortieDimdesc3"))
                           ),
                           br(), 
                           conditionalPanel(
                             condition="input.Dim=='Dim2'",
                             p("variables"),
                             div(align="center",tableOutput("sortieDimdesc22")),
                             br(),
                             p("categorical"),
                             div(align="center",tableOutput("sortieDimdesc00")),
                             br(),
                             p("quantitative"),
                             div(align="center",tableOutput("sortieDimdesc33"))
                           ),
                           br(),
                           conditionalPanel(
                             condition="input.Dim=='Dim3'",
                             p("variables"),
                             div(align="center",tableOutput("sortieDimdesc222")),
                             br(),
                             p("categorical"),
                             div(align="center",tableOutput("sortieDimdesc000")),
                             br(),
                             p("quantitative"),
                             div(align="center",tableOutput("sortieDimdesc333"))
                           )
                  ),
                  
                  tabPanel("Summary of dataset",
                           br(),
                           verbatimTextOutput("summary"),
                           selectInput("bam","Graphs for ",choices=list(IdChoices=VariableChoices),multiple=FALSE),
                           
                           div(align = "center",plotOutput("histo", width = 500, height=500))),
                  
                  
                  tabPanel("Data",
                           br(),
                           dataTableOutput("JDD")
                  )                       
      )
    )
  )
))
