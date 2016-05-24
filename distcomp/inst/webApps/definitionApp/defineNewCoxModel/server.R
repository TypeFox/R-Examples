shinyServer(function(input, output, session) {

  createProjectDefinition <- function() {
    data.frame(id = getComputationInfo("id"),
               compType = getComputationInfo("compType"),
               projectName = getComputationInfo("projectName"),
               projectDesc = getComputationInfo("projectDesc"),
               formula = getComputationInfo("formula"),
               stringsAsFactors=FALSE)
  }

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

  ## When the user chooses a file or Redcap source or Database source
  ## and clicks on the "Load Data" button this function is triggered
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
        updateTabsetPanel(session, inputId="navigationList", selected="Formula Check")
        str(dataResult)
      } else {
        cat('Error!', dataResult$message)
      }
    })
  })

  output$dataLoaded <- reactive({
    if (input$loadData == 0) return()
    ifelse(is.data.frame(getComputationInfo("data")), "Data loaded; Proceed to Formula Check", "")
  })

  ## When the user clicks on the "Check Formula" button
  ## this function is triggered
  output$checkFormulaResult <- renderPrint({
    if (input$checkFormula == 0) return()
    isolate({
      result <- tryCatch(
        { CoxWorker$new(formula = as.formula(input$formula), data=getComputationInfo("data")) },
        warning = function(x) x,
        error = function(x) x)
      if ("CoxWorker" %in% class(result)) { ## Success
        setComputationInfo("formula", input$formula)
        ## At this point, generate the definition id too.
        ## object <- list(compType = getComputationInfo("compType"),
        ##                projectName = getComputationInfo("projectName"),
        ##                projectDesc = getComputationInfo("projectDesc"),
        ##                formula = getComputationInfo("formula"),
        ##                random = runif(10))
        setComputationInfo("id", generateId(list(random=runif(10))))
        sprintf("Formula '%s' is OK!", input$formula)
      } else {
        cat("Error!", result$message)
      }
    })
  })

  output$formulaChecked <- reactive({
    if (input$checkFormula == 0) return()
    ifelse(is.null(getComputationInfo("formula")), "", "Formula Checked; Proceed to Output Result")
  })


  output$outputResult <- renderPrint({
    if (input$saveDefinition == 0) return()
    defn <- createProjectDefinition()
    str(defn)
  })

  output$definitionSaved <- reactive({
    if (input$saveDefinition == 0) return()
    defn <- createProjectDefinition()
    ##defnPath <- getConfig()$defnPath
    ##dirName <- paste(defnPath, defn$id, sep=.Platform$file.sep)
    dirName <- getComputationInfo("workingDir")
    fileName <- paste(dirName, input$outputFile, sep=.Platform$file.sep)
    result <- tryCatch(
      {
        ##dir.create(dirName)
        saveRDS(object = defn, file=fileName)
      },
      error = function(x) x)
    if (inherits(result, "error")) {
      paste("Error!", result$message)
    } else {
      paste0("Definition saved to ", input$outputFile)
    }
  })
})

