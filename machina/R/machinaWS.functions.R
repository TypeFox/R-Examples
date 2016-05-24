# Set
Sys.setenv(TZ="UTC")
machina <- new.env(parent = emptyenv())

# Private functions

#' sendAndReceive()
#'
#' Sends message via web service to kafka and receives the corresponding response
sendAndReceive <-
  function(msg, needOpenModel=TRUE, verbose=FALSE)
  {
    if (verbose) cat("sendAndReceive()\n")

    #' check state of session
    if (class(machina$session) != "machina") { stop("Session not open") }
    if (is.null(machina$session$initialized) || !machina$session$initialized)
    {
      stop("Web session not started")
    }
    if (needOpenModel)
    {
      if (is.null(machina$session$modelName)) stop("Model not open")
    }

    url <- buildURL("/api/user/kafka/sendMessageToMW")
    if (verbose) cat("POST to ", url, "\n", sep="")
    rsp <-
      POST(
        url,
        add_headers(access_token = machina$session$accessToken),
        body = msg,
        encode = "json",
        if (verbose) { verbose() })

    if (rsp$status_code != 200)
    {
      if (verbose) cat("Error return\n")
      stop("[", rsp$errorType, "] ", rsp$errorMessage)
    }

    if (verbose)
    {
      cat("Content:\n")
      print(rawToChar(rsp$content))
    }

    content <- jsonlite::fromJSON(rawToChar(rsp$content))
    if (!is.null(content$error)) stop("Error: ", content$error)

    if (verbose) cat("Returning\n")
    return(content$result)
  }

#' buildXTS()
#'
#' Builds an XTS timeseries object from a Machina row
buildXTS <-
  function(row)
  {
    for (i in 1:length(row$days$data))
    {
      temp.times <- seq(row$days$startTime[i],
                        by=row$barDuration,
                        length.out=row$days$length[i])
      temp.data <- data.frame(temp.times, row$days$data[i])
      names(temp.data) <- c("Times", "Data")
      temp.ts <- xts(temp.data$Data, order.by=temp.data$Times)
      if (i == 1)
        out.ts <- temp.ts
      else
        out.ts <- rbind(out.ts, temp.ts)
    }
    return(out.ts)
  }

#' buildURL()
#'
#' Builds service URL from session's service host plus path
buildURL <-
  function(path)
  {
    return(paste(machina$session$protocol, "://", machina$session$serviceHost, path, sep=""))
  }

#' buildRow()
#'
#' Builds row from service response
buildRow <-
  function(rsp)
  {
    list(rowIndex=rsp$rowIndex,
         query=rsp$query,
         timeSeries=rsp$timeSeries,
         offset=rsp$offset,
         length=rsp$length,
         barDuration=rsp$barDuration,
         backtestSummary=rsp$backtestSummaryResults,
         backtest=rsp$backtestResults,
         days=rsp$days)
  }

# Public functions

#' Members of Machina session (machina$session)
#'
#' Informational:
#'   $serviceHost        name/port of web service
#'   $userName           name of user who started the session
#'   $modelName          name of current model
#'   $initialized        whether session was successfully initialized
#'
#' Functional:
#'   $modelCallback      function to call when model changes
#'   $protocol           protocol to use (http or https)
#'
#' Maintained by session:
#'   $accessToken        web session token
#'   $backtestConfig     last backtest configuration set to or received from model
#'   $backtestConfigList last list of backtest configurations received from model
#'   $model              last list of rows retrieved from the model
#'   $row                last row retrieve from the model

#' startSession()
#'
#' Starts Machina session by logging in to web service
#'
#' The service uses the following URLs:
#'   http://<host>/api/user/authenticateAccount    - login
#'   http://<host>/api/user/kafka/sendMessageToMW  - interact with ModelWorker
#'   http://<host>/api/user/logout                 - logout
#'
#' @param userName user name
#' @param password user password
#' @param serviceHost host of web service to use
#' @param protocol protocol to use (default = "https")
#' @param verbose whether to produce verbose output (default = FALSE)
#
#'  Side effects: Logs in to web service
startSession <-
  function(
    userName = NULL,
    password = NULL,
    serviceHost = "account.machi.na",
    protocol = "https",
    modelCallback = viewModel,
    verbose = FALSE)
  {
    if (interactive())
    {
      if (is.null(userName))
        userName <- readline("User name: ")
      if (is.null(password))
        password <- readline("Password: ")
    }
    if (is.null(userName) || is.null(password)) stop("User id and password must be supplied")
    if (!is.null(machina$session) && machina$session$initialized) stop("Session already started")

    machina$session <- list(serviceHost=serviceHost, initialized=FALSE)
    machina$session$userName <- userName
    machina$session$protocol <- protocol
    machina$session$modelCallback <- modelCallback
    class(machina$session) <- "machina"

    url <- buildURL("/api/user/authenticateAccount")
    if (verbose) cat("POST to ", url, "\n", sep="")
    rsp <- POST(url, authenticate(userName, password), if (verbose) { verbose() })

    machina$session$initialized = (rsp$status_code == 200)
    machina$session$accessToken <- rsp$headers$access_token

    if (!machina$session$initialized) stop("[", rsp$status_code, "]: ", rawToChar(rsp$content))
  }

#' endSession()
#'
#' Ends Machina session by logging out of web service interface
#'
#' @param verbose whether to produce verbose output (default = FALSE)
#'
#' @return Side effects: Logs out of web service
endSession <-
  function(verbose = FALSE)
  {
    if (class(machina$session) != "machina") stop("Session not open")
    if (!machina$session$initialized) stop("Session not started")

    url <- buildURL("/api/user/logout")
    if (verbose) cat("POST to ", url, "\n", sep="")
    rsp <-
      POST(
        url,
        add_headers(access_token = machina$session$accessToken),
        if (verbose) { verbose() })
    machina$session$initialized = (rsp$status_code != 200)

    if (machina$session$initialized) stop("[", rsp$status_code, "]: ", rawToChar(rsp$content))

    # remove session
    rm("session", envir = machina)
  }

#' listModels()
#'
#' Lists models available to user
#
#' @param verbose whether to produce verbose output (default = FALSE)
#
#' @return data.frame with models
listModels <-
  function(verbose = FALSE)
  {
    rsp <-
      sendAndReceive(
        msg = list(messageType = "ListModels"),
        needOpenModel = FALSE,
        verbose = verbose)
    return(rsp$models)
  }

#' openModel()
#'
#' Opens ModelWorker session
#'
#' @param modelName mode name of model to open
#' @param updateModel whether to populate the saved model (default = TRUE)
#' @param verbose whether to produce verbose output (default = FALSE)
#'
#' @return row list if updateModel is TRUE, else none
#'  Side effects: session members updated are $modelName and $model (if updateModel is TRUE)
openModel <-
  function(
    modelName = NULL,
    updateModel = TRUE,
    verbose = FALSE)
  {
    if (class(machina$session) != "machina") { stop("Session not open") }
    if (!is.null(machina$session$modelName)) stop("Model already open")

    if (is.null(modelName))
    {
      machina$session$modelName <- paste("model", floor(runif(1, 0, 10^12)), sep="")
      cat("Model name not supplied. Using", machina$session$modelName, "\n")
    }
    else
      machina$session$modelName <- modelName

    if (is.null(machina$session$modelName)) stop("Model name must be supplied")

    if (updateModel)
    {
      dummy <- getModel(verbose)
      return(machina$session$model)
    }
  }

#columns we want

#      "Row"     = character(0),
#      "Vector"  = character(0),
#      "Query"   = character(0),
#      "PnL"     = character(0),
#      "Sharpe"  = character(0),
#	      "TradesPerDay" = character(0),

viewModel <-
  function(verbose = FALSE)
  {
    if ((class(machina$session$model) != "data.frame") || (nrow(machina$session$model) == 0))
      cat("Empty model\n")
    else
    {
      if (!is.null(machina$session$model$backtest))
        rows <-
          data.frame(
            Query = machina$session$model$query,
            BacktestName = machina$session$model$backtest$backtestName,
            PandL = machina$session$model$backtest$pnl,
            NumTrades = machina$session$model$backtest$ntrades,
            SharpeRatio = machina$session$model$backtest$sharpe,
            SortinoRatio = machina$session$model$backtest$sortino)
      else
        rows <- data.frame(Query = machina$session$model$query)
      if (verbose) print(rows)

      View(
        rows,
        paste("Model Rows [", machina$session$userName, "/", machina$session$modelName, "]", sep=""))
    }
  }

#' closeModel()
#'
#' Closes ModelWorker session
#'
#' @param verbose whether to produce verbose output (default = FALSE)
#'
#' @return Side effects: model name is removed from session
closeModel <-
  function(verbose = FALSE)
  {
    machina$session$modelName <- NULL
  }

#' getModel()
#'
#' Gets model, returning date.frame with information on rows in model
#' If modelCallback is defined, it is called with the updated model
#' as its parameter.
#
#' @param verbose whether to produce verbose output (default = FALSE)
#'
#' @return data.frame listing rows in model
#'   Side effects: data.frame of rows is stored in session$model
getModel <-
  function(verbose = FALSE)
  {
    rsp <-
      sendAndReceive(
        msg =
          list(
            messageType = "GetModel",
            modelName = machina$session$modelName),
        verbose = verbose)
    machina$session$model <- rsp$rows
    if (is.function(machina$session$modelCallback))
    {
      if (verbose) cat("Calling model callback\n")
      machina$session$modelCallback(verbose = verbose)
    }
    return(machina$session$model)
  }

#' clearModel()
#'
#' Clears model by removing all rows
#'
#' @param updateModel whether to update the row list (default = TRUE)
#' @param verbose whether to produce verbose output (default = FALSE)
#'
#' @return updated (empty) row list if updateModel is TRUE, else none
#'   Side effects: (empty) data.frame of rows is stored in session$model if updateModel is TRUE
clearModel <-
  function(updateModel = TRUE, verbose = FALSE)
  {
    rsp <-
      sendAndReceive(
        msg =
          list(
            messageType = "ClearModel",
            modelName = machina$session$modelName),
        verbose = verbose)
    if (updateModel)
      return(getModel(verbose))
  }

#' undoModel()
#'
#' Undoes last operation on model (addRow or clear)
#'
#' @param updateModel whether to update the row list (default = TRUE)
#' @param verbose whether to produce verbose output (default = FALSE)
#'
#' @return updated (empty) row list if updateModel is TRUE, else none
#'   Side effects: data.frame of rows is stored in session$model if updateModel is TRUE
undoModel <-
  function(updateModel = TRUE, verbose = FALSE)
  {
    rsp <-
      sendAndReceive(
        msg =
          list(
            messageType = "Undo",
            modelName = machina$session$modelName),
        verbose = verbose)
    if (updateModel)
      return(getModel(verbose))
  }

#' addRow()
#'
#' Adds row to model
#'
#' If inclusion of data is requested, on return $ts has XTS timeseries of data for the row,
#' in both return value and machina$session$row.
#' startDate and endDate default to NULL (specify beginning and end, respectively,
#' of available data).  Either or both can be set using format "YYYYMMDD" respectively,
#' in order to reduce the amount of data returned.
#'
#' @param query Machina command line query
#' @param updateModel whether to update the row list (default = TRUE)
#' @param includeData whether to include data (default = FALSE)
#' @param startDate start date (default = NULL)
#' @param endDate end date (default = NULL)
#' @param verbose whether to produce verbose output (default = FALSE)
#'
#' @return list with row info
#'   Side effects: returned list is saved in session$row
addRow <-
  function(
    query,
    updateModel = TRUE,
    includeData = FALSE,
    startDate = NULL,
    endDate = NULL,
    verbose = FALSE)
  {
    rsp <-
      sendAndReceive(
        msg =
          list(
            messageType = "AddRow",
            modelName = machina$session$modelName,
            query = query),
        verbose = verbose)
    machina$session$row <- buildRow(rsp)
    if (updateModel) dummy <- getModel(verbose)
    if (includeData)
      dummy <-
        getRow(
          machina$session$row$rowIndex,
          includeData = includeData,
          startDate = startDate,
          endDate = endDate,
          verbose = verbose)
    return(machina$session$row)
  }

#' getRow()
#'
#' Gets row, optionally with data
#'
#' If inclusion of data is requested, on return $ts has XTS timeseries of data for the row,
#' in both return value and machina$session$row.
#' startDate and endDate default to NULL (specify beginning and end, respectively,
#' of available data).  Either or both can be set using format "YYYYMMDD" respectively,
#' in order to reduce the amount of data returned.
#'
#' @param rowIndex row number to get (1-based)
#' @param includeData whether to include data (default = FALSE)
#' @param startDate start date (default = NULL)
#' @param endDate end date (default = NULL)
#' @param verbose whether to produce verbose output (default = FALSE)
#'
#' @return list with info on row
#'   Side effects: returned list is saved in session$row
getRow <-
  function(rowIndex,
           includeData = FALSE,
           startDate = NULL,
           endDate = NULL,
           verbose = FALSE)
  {
    msg <-
      list(
        messageType = "GetRow",
        modelName = machina$session$modelName,
        rowIndex = rowIndex,
        includeData = as.character(includeData))
    if (!is.null(startDate)) msg$startTime <- paste(startDate, "T0931", sep="")
    if (!is.null(endDate)) msg$endTime <- paste(endDate, "T1601", sep="")

    rsp <- sendAndReceive(msg = msg, verbose = verbose)

    # parse duration; for now assume minute bars
    rsp$barDuration <- duration(1, "minutes")
    # convert values
    rsp$days$startTime <- ymd_hm(rsp$days$startTime)
    rsp$days$length <- as.numeric(rsp$days$length)
    rsp$days$offset <- as.numeric(rsp$days$offset)

    machina$session$row <- buildRow(rsp)

    ##' build timeseries from data, if requested
    if (includeData) machina$session$row$ts <- buildXTS(rsp)

    return(machina$session$row)
  }