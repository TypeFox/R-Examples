fluidPage(#theme=shinytheme("united"),
  headerPanel(
    HTML('Distribuciones de variables aleatorias (test para UAA)
         <a href="http://snap.uaf.edu" target="_blank"><img align="right" alt="SNAP Logo" src="./img/SNAP_acronym_100px.png" /></a>'
    ), "Distributions of Random Variables"
    ),
  fluidRow(
    column(4,
           wellPanel( radioButtons("disttype","Tipo de distribución:",list("Discreta","Continua"),selected="Discreta") ),
           wellPanel(	uiOutput("distName") ),
           wellPanel(
             numericInput("n","Tamaño de la muestra:",1000),
             uiOutput("dist1"),
             uiOutput("dist2"),
             uiOutput("dist3")
           ),
           wellPanel(
             uiOutput("sampDens"),
             uiOutput("BW"),
             fluidRow(
               column(6, downloadButton("dlCurPlot", "Descargar Gráfico", class="btn-block btn-primary")),
               column(6, downloadButton("dldat", "Descargar Muestra", class="btn-block btn-warning"))
             )
           )
    ),
    column(8,
           tabsetPanel(
             tabPanel("Gráfico",plotOutput("plot", width="100%", height="auto"),verbatimTextOutput("summary")),
             #tabPanel("Summary",verbatimTextOutput("summary")),
             tabPanel("Muestra",tableOutput("table")),
             tabPanelAbout(),
             id="tsp"
           )
    )
  )
  )