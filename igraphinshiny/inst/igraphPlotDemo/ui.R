shiny::shinyUI(fluidPage(
  titlePanel("igraph plot demo"),
  fluidRow(
    column(4, wellPanel(
      fluidRow(
        column(7, wellPanel(
          selectInput(inputId="GraphType", label="Graph Type",
                      choices=c("Full Graph",
                                "Empty Graph",
                                "Star Graph",
                                "Lattice Graph",
                                "Ring Graph",
                                "Tree Graph",
                                "Erdos-Renyi",
                                "Watts-Strogatz",
                                "Barabasi-Albert",
                                "Adjacency Matrix"), selected="Full Graph")
        )),
        column(5, wellPanel(
          uiOutput("GraphTypeUI")
        ))
      )
    )),
    column(8, wellPanel(
      tabsetPanel(type="tabs",
                  tabPanel("Plot", wellPanel(
                    fluidRow(
                      column(12, wellPanel(
                        plotOutput("graphPlot")
                      ))
                    ),
                    fluidRow(
                      column(4, wellPanel(
                        radioButtons(inputId="PlotLayout", label="Plot Layout", choices=c("Auto","Random","Circle","Sphere","Fruchterman Reingold","Kamada Kawai","Drl","Spring","Reingold Tilford","Fruchterma Reingold Grid","Lgl","Graphout","SVD"), selected="Auto")
                      )),
                      column(4, wellPanel(
                        checkboxInput(inputId = "showNodeName", label = "Show Vertex Label",  value = TRUE),
                        sliderInput(inputId = "vertexSize", label = "Vertex Size",  value = 15, min=1, max=100)
                      )),
                      column(4, wellPanel(
                        downloadButton('downloadPlot', 'Download Plot in pdf')
                      ))
                    )
                  )),
                  tabPanel("Adjacency Matrix", tableOutput("AdjMatrix")),
                  tabPanel("Centralities", tableOutput("Centralities"))
      )
    ))
  )
))
