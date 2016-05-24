
shinyServer(function(input, output, session) {


  output$projectSummary <- renderText({
      if(input$printProjectSummary == 0) return()
          isolate({
              paste("Project Summary: ")
          })
    })


  output$name <- renderText({
      if(input$printProjectSummary == 0) return()
        isolate({
        paste("Project name: ", input$nameIn)
      })
    })


  output$desc <- renderText({
    if(input$printProjectSummary == 0) return()
      isolate({
        paste("Description: ", input$descIn)
    })
  })

  output$compType <- renderText({
    if (input$gotoDataInputs ==0 ) return()
    available <- availableComputations()
    availableDesc <- lapply(available, function(x) x$desc)
    chosen <- match(input$compType, availableDesc)
    compType <- names(available)[chosen]
    setComputationInfo("compType", compType)
    setComputationInfo("projectName", input$nameIn)
    setComputationInfo("projectDesc", input$descIn)
    paste(input$compType, "chosen.")
  })

  observe({
    if (input$gotoDataInputs ==0 ) return()
    available <- availableComputations()
    availableDesc <- lapply(available, function(x) x$desc)
    chosen <- match(input$compType, availableDesc)
    compType <- names(available)[chosen]
    setComputationInfo("compType", compType)
    setComputationInfo("projectName", input$nameIn)
    setComputationInfo("projectDesc", input$descIn)
    stopApp(compType)
  })
})


