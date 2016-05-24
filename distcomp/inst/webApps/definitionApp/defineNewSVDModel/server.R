shinyServer(function(input, output, session) {

  createProjectDefinition <- function() {
    data.frame(id = getComputationInfo("id"),
               compType = getComputationInfo("compType"),
               projectName = getComputationInfo("projectName"),
               projectDesc = getComputationInfo("projectDesc"),
               rank = getComputationInfo("rank"),
               ncol = getComputationInfo("ncol"),
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
  ## output$rankEntered
  ## output$rankChecked
  ## output$definitionSaved

  ## When the user chooses a file and clicks on the "Load Data" button
  ## this function is triggered
  output$dataFileContentSummary <- renderPrint({
    if (input$loadData == 0) return("")
    ## Create data frame from file
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
          { read.csv(file = inputFile$datapath, na.strings = missingValueIndicators, header=FALSE) } ,
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
        dataResult <- as.matrix(dataResult)
        setComputationInfo("data", dataResult) ## Store data object
        ncol <- ncol(dataResult)
        if (ncol <= 1) {
          "Error: need a matrix"
        } else {
          setComputationInfo("ncol", ncol) ## Store number of columns
          updateTabsetPanel(session, inputId="navigationList", selected="Rank Check")
          str(dataResult)
        }
      } else {
        cat('Error!', dataResult$message)
      }
    })
  })

  output$dataLoaded <- reactive({
    if (input$loadData == 0) return()
    ifelse(is.matrix(getComputationInfo("data")), "Data Loaded; Proceed to Rank Check", "")
  })

  ## When the user clicks on the "Check Rank" button
  ## this function is triggered
  output$checkRankResult <- renderPrint({
    if (input$checkRank == 0) return()
    isolate({
      rank <- tryCatch(as.integer(input$rank),
                       warning = function(x) x,
                       error = function(x) x)
      if (inherits(rank, "error")) {
        cat("Error!", rank$message)
      } else {
        if (rank < 0 || rank > min(dim(getComputationInfo("data")))) {
          sprintf("Rank '%d' is invalid.", rank)
        } else {
          setComputationInfo("rank", as.integer(input$rank))
          ## At this point, generate the definition id too.
          ## object <- list(compType = getComputationInfo("compType"),
          ##                projectName = getComputationInfo("projectName"),
          ##                projectDesc = getComputationInfo("projectDesc"),
          ##                rank = getComputationInfo("rank"),
          ##                random = runif(10))
          setComputationInfo("id", generateId(list(random=runif(10))))
          sprintf("Rank '%s' is OK!", input$rank)
        }
      }
    })
  })

  output$rankChecked <- reactive({
    if (input$checkRank == 0) return()
    ifelse(is.null(getComputationInfo("rank")), "", "Rank Checked; Proceed to Output Result")
  })


  output$outputResult <- renderPrint({
    if (input$saveDefinition == 0) return()
    defn <- createProjectDefinition()
    str(defn)
  })

  output$definitionSaved <- reactive({
    if (input$saveDefinition == 0) return()
    defn <- createProjectDefinition()
    ## defnPath <- getConfig()$defnPath
    ## dirName <- paste(defnPath, defn$id, sep=.Platform$file.sep)
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

