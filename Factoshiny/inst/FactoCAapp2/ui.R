# ui script for CA2

shinyUI(fluidPage(
  titlePanel(div(paste("CA on the ",nomData," dataset"),style="color:#2A0A29",align="center")),
  
  sidebarLayout(
    sidebarPanel(
      tags$head(
        tags$style("body {background-color: #F0E6E6; }"),
        tags$style(type='text/css', "#title1 { height: 25px; }")
      ),
      wellPanel(
      div(align="center",checkboxInput("caparam","Show CA parameters",FALSE)),
      conditionalPanel(
        condition="input.caparam==true",
        br(),
        if(is.null(colonnesup)){
      radioButtons("selecactive",label=h6("Choose the active columns"),
                         choices=list("All"="Toutes","Choose"="choix"),selected="Toutes")}
      else{
        radioButtons("selecactive",label=h6("Choose the active columns"),
                     choices=list("All"="Toutes","Choose"="choix"),selected="choix") 
      },
      conditionalPanel(
        condition="input.selecactive=='choix'",
        if(is.null(colonnesup)){
        selectInput("supvar",label=h6("Select the supplementary columns"),
                  choices=list(IdChoices=VariableChoices),
                  multiple=TRUE)}
        else{
          selectInput("supvar",label=h6("Select the supplementary columns"),
                      choices=list(IdChoices=VariableChoices),
                      multiple=TRUE,selected=colonnesup)  
        }
        ),
      if(is.null(catsup)){
      if(length(QualiChoice)>1){
        selectInput("supquali",label=h6("Select the supplementary categorical columns"),choices=list(Idqualisup=as.vector(QualiChoice)),multiple=TRUE)
      }
      if (length(QualiChoice)==1){
        h6("Select the supplementary categorical columns")
        checkboxInput("supquali",QualiChoice,FALSE)
      }}
      else{
        if(length(QualiChoice)>1){
          selectInput("supquali",label=h6("Select the supplementary categorical columns"),choices=list(Idqualisup=as.vector(QualiChoice)),multiple=TRUE,selected=catsup)
        }
        if (length(QualiChoice)==1){
          h6("Select the supplementary categorical columns")
          checkboxInput("supquali",QualiChoice,TRUE)
        } 
      },
      br(),
      if(is.null(lignesup)){
        selectInput("rowsupl",label=h6("Select the supplementary rows"),choices=list(num=nom),multiple=TRUE)}
      else{
        selectInput("rowsupl",label=h6("Select the supplementary rows"),choices=list(num=nom),multiple=TRUE,selected=lignesup)
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
        textInput("title1",h6("Title of the graph : "),title1),
        if(is.null(Invisible)){
          selectInput("invis",h6("Invisible points"),choices=list("Rows"="row","Columns"="col","Supplementary rows"="row.sup","Supplementary columns"="col.sup","Supplementary qualitative variable"="quali.sup"),multiple=TRUE)}
        else{
          selectInput("invis",h6("Invisible points"),choices=list("Rows"="row","Columns"="col","Supplementary rows"="row.sup","Supplementary columns"="col.sup","Supplementary qualitative variable"="quali.sup"),multiple=TRUE,selected=Invisible)
        },
        br(),
        sliderInput("cex",h6("Size of labels"),min=0.5,max=1.5,value=size,step=0.05),
        br(),
        radioButtons("seleccol",h6("Draw columns according to :"), choices=list("No selection"="no","Cos2"="cos2","Contribution"="contrib"),selected=selec1,inline=TRUE),
        conditionalPanel(
          condition="input.seleccol=='cos2'",
          if(selec1=="cos2"){
            sliderInput("slider3",label=h6("Select colums that have a cos2 greater than :"),
                        min=0,max=1,value=valueselec1,step=0.05)
          }
          else{
            sliderInput("slider3",label=h6("Select colums that have a cos2 greater than :"),
                        min=0,max=1,value=0,step=0.05) 
          }),
        conditionalPanel(
          condition="input.seleccol=='contrib'",
          uiOutput("contribcol")),
        br(),
        radioButtons("selecrow",h6("Draw rows according to :"), choices=list("No selection"="no","Cos2"="cos2","Contribution"="contrib"),selected=selec2,inline=TRUE),
        conditionalPanel(
          condition="input.selecrow=='cos2'",
          if(selec2=="cos2"){sliderInput("slider4",label=h6("Select rows that have a cos2 greater than :"),
                                         min=0,max=1,value=valueselec2,step=0.05)}
          else{sliderInput("slider4",label=h6("Select rows that have a cos2 greater than :"),
                           min=0,max=1,value=0,step=0.05)}),
        conditionalPanel(
          condition="input.selecrow=='contrib'",
          uiOutput("contribrow"))
      )
      ),
      wellPanel(
        h5("Save graphs as : ",align="center"),
        radioButtons("paramdown","",
                     choices=list("PNG"="png","JPG"="jpg","PDF"="pdf"),selected="png"),
      br(),
      div(align="center",actionButton("CAcode", "Get the CA code"))),
    div(align="center",actionButton("Quit", "Quit the app"))
    ),  
      mainPanel(
        tags$style(type = "text/css", "a{color: #2F0B3A;}"),
        tabsetPanel(id = "graph_sort",
                    tabPanel("Graphs",
                             br(),
                             div(verbatimTextOutput("warn")),
                             br(),
                             div(align="center",plotOutput("map",width=550,height=450)),
                             br(),
                             conditionalPanel(
                               condition="input.paramdown=='jpg'",
                               p(downloadButton("downloadData1","Download as jpg"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='png'",
                               p(downloadButton("downloadData","Download as png"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='pdf'",
                               p(downloadButton("downloadData2","Download as pdf"),align="center"))),
                    tabPanel("Values",
                             br(),
                             uiOutput("out22"),
                             br(),
                             conditionalPanel(
                               condition="input.out=='eig'",
                               div(align="center",tableOutput("sorties")),
                               div(align="center",plotOutput("map3"))),
                             conditionalPanel(
                               condition="input.out=='CA'",
                               numericInput("nbele",h6("Number of elements to print"),value=10),
                               verbatimTextOutput("summaryCA"),
                               p(downloadButton("summary2","Download the summary"),align="center")
                               ),
                             conditionalPanel(
                               condition="input.out=='var'",
                               h6("Coordinates"),
                               div(align="center",tableOutput("sorties1")),
                               h6("cos2"),
                               div(align="center",tableOutput("sorties2")),
                               h6("Contributions"),
                               div(align="center",tableOutput("sorties3"))
                               ),
                             conditionalPanel(
                               condition="input.out=='ind'",
                               h6("Coordinates"),
                               div(align="center",tableOutput("sorties4")),
                               h6("cos2"),
                               div(align="center",tableOutput("sorties5")),
                               h6("Contributions"),
                               div(align="center",tableOutput("sorties6"))
                             ),
                             conditionalPanel(
                               condition="input.out=='suprow'",
                               h6("Coordinates"),
                               div(align="center",tableOutput("sorties7")),
                               h6("cos2"),
                               div(align="center",tableOutput("sorties8"))
                             ),
                             conditionalPanel(
                               condition="input.out=='supcol'",
                               h6("Coordinates"),
                               div(align="center",tableOutput("sorties9")),
                               h6("cos2"),
                               div(align="center",tableOutput("sorties10"))
                             ),
                             conditionalPanel(
                               condition="input.out=='qualico'",
                               div(align="center",tableOutput("sorties11"))
                             )
                             ),
                    tabPanel("Summary of dataset",
                             br(),
                             verbatimTextOutput("summary")),
                    
                    tabPanel("Data",
                             br(),
                             dataTableOutput("JDD")
                             )
        )
      )
    )
))
