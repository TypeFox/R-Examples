shinyServer(function(input, output, session) {

  observe({
    if (input$exitApp > 0) stopApp(TRUE)
  })

  ## -- End: functions to make tabs active sequentially --

  ## Variables to detect various states the app can be in, depending
  ## on what actions the user has taken.

  ## output$datasetSpecified
  ## output$datasetChecked
  ## output$formulaEntered
  ## output$formulaChecked
  ## output$definitionSaved

  ## When the user chooses a file and clicks on the "Load Data" button
  ## this function is triggered
  output$dataFileContentSummary <- renderPrint({
    if (input$loadData == 0) return("")

    ## Create data frame from source
    isolate({
      if (input$input_type == 'CSV File') {
        inputFile <- input$dataFile
        shiny::validate(
          need(inputFile != "", "Please select a data set")
        )
        ## Parse missing value strings
        missingValueIndicators <- stringr::str_trim(scan(textConnection(input$missingIndicators),
                                                         what = character(0), sep=",", quiet=TRUE))
        ## Return data frame or error as the case may be
        dataResult <- tryCatch(
          { read.csv(file = inputFile$datapath, na.strings = missingValueIndicators) } ,
          warning = function(x) x,
          error = function(x) x
        )
      } else if (input$input_type == 'Redcap API') {
        shiny::validate(
          need(requireNamespace("redcapAPI", quietly = TRUE), "Please install the redcapAPI package"),
          need(input$redcapURL != "", "Please enter your Redcap URL"),
          need(input$redcapToken != "", "Please enter your Redcap Token")
        )

        dataResult <-tryCatch(
          { redcapAPI::exportRecords(redcapAPI::redcapConnection(url = input$redcapURL, token = input$redcapToken)) },
          warning = function(x) x,
          error = function(x) x
        )
      } else if (input$input_type == 'Postgres') {
        shiny::validate(
          need(requireNamespace("RPostgreSQL", quietly = TRUE), "Please install the RPostgreSQL package"),
          need(requireNamespace("dplyr", quietly = TRUE), "Please install the dplyr package"),
          need(input$dbName != "", "Please enter your Postgres database name"),
          need(input$dbHost != "", "Please enter your Postgres host"),
          need(!is.na(as.integer(input$dbPort)), "Please enter your Postgres port number"),
          need(input$dbUser != "", "Please enter your Postgres database user name"),
          need(input$dbPassword != "", "Please enter your Postgres database user password"),
          need(input$dbTable != "", "Please enter your Postgres database table name")
        )
        dataResult <-  tryCatch(
          {
            db <- dplyr::src_postgres(dbname = input$dbName, host = input$dbHost, port = input$dbPort,
                                      user = input$dbUser, password = input$dbPassword)
            ## CAUTION: Need to do better to prevent SQL injection...
            dplyr::tbl(db, paste("SELECT * from ", input$dbTable))
          },
          warning = function(x) x,
          error = function(x) x
        )
      } else {
        shiny::validate(
          need(FALSE, "Report bug to Package owner: unexpected Data source!")
        )
        dataResult <- NULL
      }

      if (is.data.frame(dataResult)){
        setComputationInfo("data", dataResult) ## Store data object
        updateTabsetPanel(session, inputId="navigationList", selected="Sanity Check")
        str(dataResult)
      } else {
        cat('Error!', dataResult$message)
      }
    })
  })

  output$dataLoaded <- reactive({
    if (input$loadData == 0) return()
    ifelse(is.data.frame(getComputationInfo("data")), "Data Loaded; Proceed to Sanity Check", "")
  })

  ## When the user clicks on the "Check Sanity" button
  ## this function is triggered
  output$sanityCheckResult <- reactive({
    if (input$checkSanity == 0) return()
    formula <- getComputationInfo("formula")
    isolate({
      result <- tryCatch(
        { CoxWorker$new(formula = as.formula(formula), data=getComputationInfo("data")) },
        warning = function(x) x,
        error = function(x) x)
      if ("CoxWorker" %in% class(result)) { ## Success
        "Success: data matches formula. Send to Opencpu Server."
      } else {
        paste("Error!", result$message)
      }
    })
  })

  output$populateResult <- renderPrint({
    if (input$populateServer == 0) return()
    isolate({
      shiny::validate(need(input$siteName != "", "Please enter a site name"))
      shiny::validate(need(input$ocpuURL != "", "Please enter an opencpu URL"))
      site <- list(name = input$siteName, url = input$ocpuURL)
      defn <- makeDefinition(getComputationInfo("compType"))
      data <- getComputationInfo("data")
      result <- tryCatch(uploadNewComputation(site=site, defn=defn, data=data),
                         error = function(x) x,
                         warning = function(x) x)
      if (inherits(result, "error") ) {
        "Error: Uploading the definition to server"
      } else {
        ifelse(result, "Success: definition uploaded to server",
               "Error while uploading definition to server")
      }
    })
  })

})

