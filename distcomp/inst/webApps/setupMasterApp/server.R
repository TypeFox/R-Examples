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

  output$definitionUploaded <- reactive({
    if (input$uploadDefinition == 0) return()
    ifelse(!is.null(getComputationInfo("id")), "Definition Uploaded; Proceed to specify sites", "")
  })

  output$writeRCode <- reactive({
    if (input$addSite == 0) return()
    ifelse(length(getComputationInfo("siteList")) > 1, "At least two sites added; so can write R code", "")
  })

  output$siteList <- renderPrint({
    if (input$addSite == 0) return()
    isolate({
      name <- stringr::str_trim(input$siteName)
      url <- stringr::str_to_lower(stringr::str_trim(input$ocpuURL))
      shiny::validate(
        need(name != "", "Please enter a non-empty name")
      )
      shiny::validate(
        need(url != "" && grepl("^http", url), "Please enter a URL")
      )

      siteList <- getComputationInfo("siteList")

      shiny::validate(
        need(!(name %in% names(siteList)), paste("Bad site: duplicate name", name))
      )

      n <- length(siteList)
      siteList[[n + 1]] <- list(name = name, url = url)
      setComputationInfo("siteList", siteList)
      str(siteList)
    })
  })

  output$codeSaved <- reactive({
    if (input$saveCode ==0 ) return()
    shiny::validate(
      need(stringr::str_trim(input$outputFile) != "", "Please enter a non-empty name")
    )

    sites <- getComputationInfo("siteList")
    defn <- makeDefinition(getComputationInfo("compType"))
    result <- tryCatch(writeCode(defn, sites, input$outputFile)
                     , error=function(x) x
                     , warning=function(x) x)

    if (inherits(result, "error")) {
      "Error!"
    } else {
      paste("Success: code saved to", input$outputFile)
    }

  })

  observe({
    if (input$exitApp ==0 ) return()
    stopApp(TRUE)
  })

})

