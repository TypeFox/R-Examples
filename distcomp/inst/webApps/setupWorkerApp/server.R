shinyServer(function(input, output, session) {

  output$definitionSummary <- renderPrint({
    if (input$uploadDefinition == 0) return("")
    ## Create data frame from file
    isolate({
      inputFile <- input$definitionFile
      ## shiny::validate(
      ##   need(inputFile != "", "Please select a data set")
      ## )

      available <- availableComputations()

      ## Try to read the RDS fileReturn data frame or error as the case may be
      definition <- tryCatch(
        { readRDS(inputFile$datapath) },
        warning = function(x) x,
        error = function(x) x
      )

      if ( inherits(definition, "error") ||
          !is.data.frame(definition) || is.null(definition$compType) ||
          !(definition$compType %in% names(available)) ) {
        cat("Error! Bad definition file")
      } else {
        for (x in names(definition)) {
          setComputationInfo(x, definition[[x]])
        }
        str(definition)
      }
    })
  })

  observe({
    if (input$gotoDataInputs ==0 ) return()
    availableNames <- names(availableComputations())
    chosen <- match(getComputationInfo("compType"), availableNames)
    compType <- availableNames[chosen]
    stopApp(compType)
  })
})

