## I still keep load these two package for stand alone shiny apps
library(shiny)
library(igraph)

shinyServer(function(input, output) {
  output$GraphTypeUI <- renderUI({
    if (is.null(input$GraphType))
      return()
    switch(input$GraphType,
           "Full Graph" = wellPanel(numericInput("nNode", "Number of vertices", value = 10, max=100),
                            checkboxInput(inputId="isDirected", label="directed", value=FALSE),
                            checkboxInput(inputId="isLoops", label="loops", value=FALSE)),
           "Empty Graph" = wellPanel(numericInput("nNode", "Number of vertices", value = 10, max=100),
                             checkboxInput("isDirected", "directed", value=TRUE)),
           "Star Graph" = wellPanel(numericInput("nNode", "Number of vertices", value = 10, max=100),
                            selectInput(inputId="GraphMode", label="Graph Mode", choices=c("undirected","in","out","mutual"), selected="undirected")),
           "Lattice Graph" = wellPanel(numericInput("dimGraph", "Dimension of Lattice", value = 2),
                                       numericInput("lengthGraph", "Length of Lattice", value = 5),
                                       checkboxInput("isDirected", "directed", value=FALSE),
                                       checkboxInput("isMutual", "mutual", value=FALSE),
                                       checkboxInput("isCircular", "circular", value=TRUE)),
           "Ring Graph" = wellPanel(numericInput("nNode", "Number of vertices", value = 10, max=100),
                            checkboxInput("isDirected", "directed", value=FALSE),
                            checkboxInput("isMutual", "mutual", value=FALSE),
                            checkboxInput("isCircular", "circular", value=TRUE)),
           "Tree Graph" = wellPanel(numericInput("nNode", "Number of vertices", value = 10, max=100),
                                    numericInput("nChild", "Number of childern", value = 2, max=10),
                                    selectInput(inputId="GraphMode", label="Graph Mode", choices=c("undirected","in","out"), selected="undirected")),
           "Erdos-Renyi" = wellPanel(numericInput("nNode", "Number of vertices", value = 10, max=100),
                                      numericInput("pNode", "Link probability", value = 0.3, max=1),
                                      checkboxInput("isDirected", "directed", value=FALSE),
                                      checkboxInput(inputId="isLoops", label="loops", value=FALSE)),
           "Watts-Strogatz" =  wellPanel(numericInput("dimNode", "Dim of the lattice", value = 1, max=2),
                                          numericInput("sizeNode", "Size in dimension", value = 50),
                                          numericInput("neiNode", "Neighborhood", value = 3),
                                          numericInput("pNode", "Wire probability", value = 0.05, max=1),
                                          checkboxInput("isMultiple", "multiple", value=FALSE),
                                          checkboxInput(inputId="isLoops", label="loops", value=FALSE)),
           "Barabasi-Albert" = wellPanel(numericInput("nNode", "Number of vertices", value = 10, max=100),
                                          numericInput("powerGraph", "Preferential attachment power", value = 1),
                                          checkboxInput("isDirected", "directed", value=FALSE)),
           "Adjacency Matrix" = wellPanel(fileInput('file1', 'Choose CSV File', accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                          checkboxInput('header', 'Header', TRUE),
                                          radioButtons('sep', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
                                          radioButtons('quote', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"),'"')),
           "Adjacency Matrix Excel" = wellPanel(fileInput('file1', 'Choose xlsx File', accept=c('text/xlsx', 'text/comma-separated-values,text/plain', '.xlsx')),
                                          checkboxInput('header', 'Header', TRUE))
          )
  })

  realTimeGraph <- reactive({
    if (is.null(input$GraphType))
      return()
    g <- switch(input$GraphType,
                "Full Graph" = graph.full(n=input$nNode, directed=input$isDirected, loops=input$isLoops),
                "Empty Graph" = graph.empty(n=input$nNode, directed=input$isDirected),
                "Star Graph" = graph.star(n=input$nNode, mode=input$GraphMode),
                "Lattice Graph" = graph.lattice(length = input$lengthGraph, dim = input$dimGraph, directed=input$isDirected, mutual=input$isMutual, circular=input$isCircular),
                "Ring Graph" = graph.ring(n=input$nNode, directed=input$isDirected, mutual=input$isMutual, circular=input$isCircular),
                "Tree Graph" = graph.tree(n=input$nNode, children=input$nChild, mode=input$GraphMode),
                "Erdos-Renyi" = erdos.renyi.game(n=input$nNode, p=input$pNode, directed=input$isDirected, loops=input$isLoops),
                "Watts-Strogatz" = watts.strogatz.game(dim=input$dimNode, size=input$sizeNode, nei=input$nei$Node, p=input$pNode, multiple=input$isMultiple, loops=input$isLoops),
                "Barabasi-Albert" = barabasi.game(n=input$nNode, power=input$powerGraph, directed=input$isDirected),
                "Adjacency Matrix" = graph.adjacency(as.matrix(read.csv(input$file1$datapath, header=input$header, sep=input$sep, quote=input$quote)))
                )
  })

  plotGraph <- function(){
    g <- realTimeGraph()
    plotlayout <- switch(input$PlotLayout,
                      "Auto"=layout.auto(g),
                      "Random"=layout.random(g),
                      "Circle"=layout.circle(g),
                      "Sphere"=layout.sphere(g),
                      "Fruchterman Reingold"=layout.fruchterman.reingold(g),
                      "Kamada Kawai"=layout.kamada.kawai(g),
                      "Drl"=layout.drl(g),
                      "Spring"=layout.spring(g),
                      "Reingold Tilford"=layout.reingold.tilford(g),
                      "Fruchterma Reingold Grid"=layout.fruchterman.reingold.grid(g),
                      "Lgl"=layout.lgl(g),
                      "Graphopt"=layout.graphopt(g),
                      "SVD"=layout.svd(g)
                      )
    if(!input$showNodeName){
      V(g)$label = ""
    }
    V(g)$size = input$vertexSize
    plot(g, layout=plotlayout)
  }

  calculateCentrality <- function(){
    Centralities <- list()
    Centralities$Alpha <- as(alpha.centrality(realTimeGraph()),"matrix")
    Centralities$Bon <- as(bonpow(realTimeGraph()),"matrix")
    Centralities$Closeness <- as(closeness(realTimeGraph()),"matrix")
    Centralities$Evcent <- as(evcent(realTimeGraph())$vector,"matrix")
    Centralities$Kleinberg <- as(authority.score(realTimeGraph())$vector,"matrix")
    Centralities$PageRank <- as(page.rank(realTimeGraph())$vector,"matrix")
    Centralities$Betweenness <- as(betweenness(realTimeGraph()),"matrix")
    return(as.data.frame(Centralities))
  }

  output$graphPlot <- renderPlot({
    suppressWarnings(plotGraph())
  })

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste('NetworkGraph',format(Sys.time(),"%Y%m%d_%H%M%S"),'.pdf', sep='')
    },
    content = function(file) {
      pdf(file)
      print(suppressWarnings(plotGraph()))
      dev.off()
    }
  )

  output$AdjMatrix <- renderTable(as(get.adjacency(realTimeGraph()),"matrix"))

  output$Centralities <- renderTable(as(calculateCentrality(),"matrix"))
})
