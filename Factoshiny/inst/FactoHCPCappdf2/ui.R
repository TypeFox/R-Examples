# ui script for HCPC for dataframe2

shinyUI(fluidPage(
  titlePanel(div(paste("HCPC on the ",nomData," dataset"),style="color:#0A2A12",align="center")),
#  titlePanel(div(paste("HCPC on the ",unlist(strsplit(nomData, split='[', fixed=TRUE))[1]," dataset"),style="color:#0A2A12",align="center")),
  
  sidebarLayout(
    sidebarPanel(
      tags$head(
        tags$style("body {background-color: #D2FAE5; }"),
        tags$style(type='text/css', "#title1 { height: 25px; }"),
        tags$style(type='text/css', "#title2 { height: 25px; }"),
        tags$style(type='text/css', "#title3 { height: 25px; }")
      ),
      wellPanel(
      div(align="center",checkboxInput("hcpcparam","Show HCPC parameters",FALSE)),
      conditionalPanel(
        condition="input.hcpcparam==true",
      uiOutput("clusters"),
      hr(),
      checkboxInput("consoli","Consolidation",consolidf),
      hr(),
      radioButtons("metric","Which metric would you like to use ?",choices=list("euclidean"="euc","manhattan"="manh"),inline=TRUE,select=metricdf)
      )),
      wellPanel(
      div(align="center",checkboxInput("graph","Show graphs options",FALSE)),
      conditionalPanel(
        condition="input.graph==true",
      fluidRow(
        column(5,selectInput("nb1", label = h6("x axis"), 
                             choices = list("1" = 1, "2" = 2, "3" = 3,"4"= 4,"5" =5), selected = nb1df,width='80%')),
        column(5,selectInput("nb2", label =h6("y axis"), 
                             choices = list("1" = 1, "2" = 2,"3" = 3,"4"= 4,"5" =5), selected = nb2df,width='80%'))),
      hr(),
      radioButtons("HCPCgraph",h6("Which graph do you want to modify ?"),
                   choices=list("Individual map"="ind","3D map"="3D","Tree map"="tree"),inline=TRUE),
      hr(),
      conditionalPanel(
        condition="input.HCPCgraph=='3D'",
        textInput("title1",h6("Title of the graph : "), title1),
        checkboxInput("nom3D","Names on 3D map",df),
        checkboxInput("center","Draw centers of clusters",centerdf),
        hr(),
        sliderInput("num","Angle (in degrees)",value=numdf,min=0,max=360,step=1)
        ),
      conditionalPanel(
        condition="input.HCPCgraph=='ind'",
        textInput("title2",h6("Title of the graph : "), title2),
        checkboxInput("drawtree","Draw tree",drawdf)
        ),
      conditionalPanel(
        condition="input.HCPCgraph=='tree'",
        textInput("title3",h6("Title of the graph : "), title3))
      )),
      wellPanel(
        h5("Save graphs as :",align="center"),
        radioButtons("paramdown","",
                     choices=list("PNG"="png","JPG"="jpg","PDF"="pdf"),selected="png"),
        br(),
        div(align="center",actionButton("HCPCcode", "Get the HCPC code"))),
      div(align="center",actionButton("Quit", "Quit the app"))
      ),
      
      mainPanel(
        tags$style(type = "text/css", "a{color: #0B6121;}"),
        tabsetPanel(id = "graph_sort",
                    tabPanel("Graphs",
                             br(),
                             div(align="center",plotOutput("map4",width = 650, height=500)),
                             br(),
                             conditionalPanel(
                               condition="input.paramdown=='jpg'",
                               p(downloadButton("downloadData6","Download as jpg"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='png'",
                               p(downloadButton("downloadData7","Download as png"),align="center")),
                             conditionalPanel(
                               condition="input.paramdown=='pdf'",
                               p(downloadButton("downloadData8","Download as pdf"),align="center")),
                             br(),
                             div(align="center",plotOutput("map",width = 500, height=500)),
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
                             div(align="center",plotOutput("map2",width=750,height=500)),
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
                             br()
                             ),

                    tabPanel("Values",
                             br(),
                             radioButtons("out"," Which outputs do you want?",
                                          choices=list("Description of classes by variables"="var","Description of classes by axes"="axe","Parangons"="para"),selected="var",inline=TRUE),
                             conditionalPanel(
                               condition="input.out=='var'",
                               div(align="center",tableOutput("descript"))
                               ),
                             conditionalPanel(
                               condition="input.out=='para'",
                               div(align="center",tableOutput("parangons"))),
                             conditionalPanel(
                               condition="input.out=='axe'",
                               div(align="center",tableOutput("axes")))
                             ),
                    
                    tabPanel("Data",
                             dataTableOutput("JDD")
                             )
        )
      )
    )
))
