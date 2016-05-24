# ui.R AFM2

shinyUI(fluidPage(
  titlePanel(div(paste("MFA on the ",nameJDD," dataset"),style="color:#6E6E6E",align="center")),
  sidebarLayout(
      sidebarPanel(
        tags$head(
          tags$style("body {background-color: #E1D3FB; }"),
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
          checkboxInput("meanind1","Draw points for the mean individuals",ind1),
          checkboxInput("meanind","Draw labels for the mean individuals",ind2),
          checkboxInput("qualind1","Draw points for the qualitative individuals",ind3),
          checkboxInput("qualind","Draw labels for the qualitative individuals",ind4),
          hr(),
          uiOutput("drawindiv"),
          conditionalPanel(
            condition="input.drawind=='c'",
            uiOutput("habillagequali")
            ),
          hr(),
          radioButtons("choixpartial",h6("Partial points are drawn"),choices=list("None"=1,"All"=2,"Choose"=3),selected=partial,inline=TRUE),
          conditionalPanel(
            condition="input.choixpartial==3",
            uiOutput("indivpartiel2")),
          hr(),
          conditionalPanel(
            condition="input.choixpartial!=1",
            checkboxInput("partind","Draw labels for the partial individuals",partial3))
          ),
        conditionalPanel(
          textInput("title3",h6("Title of the graph : "), title3),
          condition="input.choixgraph=='quant'",
          radioButtons("selection",h6("Select from"),choices=list("No selection"="no","Contribution"="contrib","Cos2"="cos2"),selected=selectvar),
          uiOutput("slider1"),
          hr(),
          uiOutput("hide2"),
          checkboxInput("colorgroup","Color the variables by group",colorvar) 
        ),
        conditionalPanel(
          condition="input.choixgraph=='freq'",
          textInput("title5",h6("Title of the graph : "), title5),
          checkboxInput("affichind","Draw labels for the mean individuals",freq1),
          checkboxInput("affichcol","Draw labels for the columns",freq2)
        ),
        conditionalPanel(
          condition="input.choixgraph=='axes'",
          textInput("title4",h6("Title of the graph : "), title4),
          checkboxInput("coloraxe","Color the partial axe by group",partaxe)
        ),
        conditionalPanel(
          condition="input.choixgraph=='group'",
          textInput("title1",h6("Title of the graph : "), title1)
        ),
        fluidRow(
          column(5,selectInput("nb1", label = h6("  x axis"), 
                               choices = list("1" = 1, "2" = 2, "3" = 3,"4"= 4,"5" =5), selected = axe1,width='80%')),
          column(5,selectInput("nb2", label =h6("   y axis"), 
                               choices = list("1" = 1, "2" = 2,"3" = 3,"4"= 4,"5" =5), selected = axe2,width='80%')))
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
                                                       "Results for the quantitative variables"="quantvar","Results for the groups"="group","Results for the partial axes"="partaxe"),selected="MFA",inline=TRUE),
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
                             verbatimTextOutput("summary")),
                    
                    tabPanel("Data",
                             br(),
                             dataTableOutput("JDD")
                             )
        )
      )
    )
))
