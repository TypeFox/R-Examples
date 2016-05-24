

library(shiny)
library(geneSLOPE)

options(shiny.maxRequestSize=300*1024^2)

shinyServer(function(input, output) {

  phenotype <- reactive({
    name <- input$fileY
    validate(
      need(name != "", "Please upload phenotype data")
    )
    read_phenotype(name$datapath, sep=input$sep, header = input$header)
  })

  output$phenotypeOk <- reactive({
    length(phenotype()$y)
  })

  screening <- eventReactive(input$go, {
    DataSet <- input$file
    Map <- input$map.file
    validate(
      need(DataSet != "", "Please upload snp data")
    )
    screen_snps(DataSet$datapath, Map$datapath, phenotype = phenotype(),
                pValMax = input$pValCutoff, chunkSize = 1e2, verbose = FALSE)
  })

  output$summary <- renderText({
    paste("Phenotype data loaded.", length(phenotype()$y), "observations.")
  })

  clumping <-  eventReactive(input$go, {
    clump_snps(screening(), input$rho, verbose = FALSE)
  })

  slopeResult <- eventReactive(input$go, {
    select_snps(clumping(), input$fdr, verbose = FALSE)
  })

  output$clumpSummary <- renderPrint({
    summary(clumping())
  })

  output$slopePlot <- renderPlot({
    plot(slopeResult())
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("results_genSlope", format(Sys.time(),  "%y_%m_%y_%H_%M_%S"), ".csv")
    },
    content = function(file) {
      write.csv(slopeResult()$X, file, row.names = FALSE)
    }
  )
})
