#" @aliases get.safeObj,safeObj,character,list-method
#" @rdname get.safeObj-method
setMethod(f="get.safeObj", signature=c("safeObj", "character", "list"),
  definition=function(object, type, input) {
    if ( !type %in% c("dimInfo", "elapsedTime", "finalData", "nrNonDuplicatedCells",
        "nrPrimSupps", "nrSecondSupps", "nrPublishableCells", "suppMethod",
        "cellInfo", "cellID") ) {
      stop("get.safeObj:: argument 'type' is not valid!\n")
    }

    if ( type == "dimInfo" ) {
      return(g_dimInfo(object))
    }
    if ( type == "elapsedTime" ) {
      return(g_elapsedTime(object))
    }
    if ( type == "finalData" ) {
      return(g_finalData(object))
    }
    if ( type == "nrNonDuplicatedCells" ) {
      return(g_nrNonDuplicatedCells(object))
    }
    if ( type == "nrPrimSupps" ) {
      return(g_nrPrimSupps(object))
    }
    if ( type == "nrSecondSupps" ) {
      return(g_nrSecondSupps(object))
    }
    if ( type == "nrPublishableCells" ) {
      return(g_nrPublishableCells(object))
    }
    if ( type == "suppMethod" ) {
      return(g_suppMethod(object))
    }
    if ( type == "cellInfo" ) {
      return(g_getCellInfo(object, input))
    }
    if ( type == "cellID" ) {
      return(g_getCellID(object, input))
    }
  }
)

#' summarize \code{\link{safeObj-class}} objects
#'
#' extract and show information stored in \code{\link{safeObj-class}} objects
#' @param object an object of class \code{\link{safeObj-class}}
#' @param ... additional arguments, currently ignored
#' @export
#' @docType methods
#' @aliases summary,safeObj-method
#' @rdname summary.safeObj-method
setMethod(f="summary", signature="safeObj",
  definition=function (object, ...) {
    cat("\n######################################################\n")
    cat("### Summary of the result object of class 'safeObj' ###\n")
    cat("######################################################\n")
    cat(paste("--> The input data have been protected using algorithm ",g_suppMethod(object),".\n", sep=""))
    cat(paste("--> The algorithm ran for ",formatTime(g_elapsedTime(object))$time.str,".\n", sep=""))
    cat(paste("--> To protect ",g_nrPrimSupps(object)," primary sensitive cells, ", g_nrSecondSupps(object)," cells need to be additionally suppressed.\n", sep=""))
    cat(paste("--> A total of ",g_nrPublishableCells(object)," cells may be published.\n", sep=""))

    finalData <- g_finalData(object)
    nrCells <- nrow(finalData)
    nrNonDups <- g_nrNonDuplicatedCells(object)
    if ( nrCells > nrNonDups ) {
      cat(paste("--> Duplicated cells: Only ",nrNonDups," table cells are unique, the remaining ",nrCells-nrNonDups," cells are duplicates.\n", sep=""))
    }
    cat("\n###################################\n")
    cat("### Structure of protected Data ###\n")
    cat("###################################\n")
    print(str(finalData))
  }
)

#' show \code{\link{safeObj-class}} objects
#'
#' extract and show information stored in \code{\link{safeObj-class}} objects
#' @param object an object of class \code{\link{safeObj-class}}
#' @export
#' @docType methods
#' @aliases show,safeObj-method
#' @rdname show.safeObj-method
setMethod(f="show", signature="safeObj",
  definition=function (object) { print(str(object)) }
)

setMethod(f="g_dimInfo", signature="safeObj", definition=function(object) {
  object@dimInfo
})

setMethod(f="g_elapsedTime", signature="safeObj", definition=function(object) {
  object@elapsedTime
})

setMethod(f="g_finalData", signature="safeObj", definition=function(object) {
  object@finalData
})

setMethod(f="g_nrNonDuplicatedCells", signature="safeObj", definition=function(object) {
  object@nrNonDuplicatedCells
})

setMethod(f="g_nrPrimSupps", signature="safeObj", definition=function(object) {
  object@nrPrimSupps
})

setMethod(f="g_nrSecondSupps", signature="safeObj", definition=function(object) {
  object@nrSecondSupps
})

setMethod(f="g_nrPublishableCells", signature="safeObj", definition=function(object) {
  object@nrPublishableCells
})

setMethod(f="g_suppMethod", signature="safeObj", definition=function(object) {
  object@suppMethod
})

setMethod(f="g_getCellInfo", signature="safeObj", definition=function(object, input) {
  cellID <- g_getCellID(object, input=input)

  primSupp <- secondSupp <- NULL
  finalData <- g_finalData(object)
  cellStatus <- finalData[cellID, "sdcStatus"]

  if ( cellStatus == "u" ) {
    primSupp <- TRUE
    secondSupp <- FALSE
    if ( input[[3]] ) {
      cat ("The cell is a sensitive cell!\n")
    }
  }
  if ( cellStatus == "s" ) {
    primSupp <- FALSE
    secondSupp <- FALSE
    if ( input[[3]] ) {
      cat ("The cell can be published!\n")
    }
  }
  if ( cellStatus == "z" ) {
    primSupp <- FALSE
    secondSupp <- FALSE
    if ( input[[3]] ) {
      cat ("The cell will be enforced for publication!\n")
    }
  }
  if ( cellStatus == "x" ) {
    primSupp <- FALSE
    secondSupp <- TRUE
    if ( input[[3]] ) {
      cat ("The cell has been secondary suppressed!\n")
    }
  }
  return(list(cellID=cellID, data=finalData[cellID,], primSupp=primSupp, secondSupp=secondSupp))
})

setMethod(f="g_getCellID", signature="safeObj", definition=function(object, input) {
  para.names <- input[[1]]
  para.codes <- input[[2]]
  para.verbose <- input[[3]]

  finalData <- g_finalData(object)
  dimInfo <- g_dimInfo(object)

  vNames <- g_varname(dimInfo)
  vIndex <- g_pos_index(dimInfo)

  indexVar <- match(para.names, vNames)

  if ( length(input) != 3 ) {
    stop("cellID:: length of argument 'input' must equal 3!\n")
  }
  if ( length(para.names) != length(para.codes) ) {
    stop("cellID:: check argument 'input'!\n")
  }
  if ( !all(para.names %in% vNames) ) {
    stop("cellID:: check variable names in 'input[[1]]'!\n")
  }
  if ( !is.logical(para.verbose) ) {
    stop("cellID:: argument in 'input[[3]]' must be logical!\n")
  }

  cellID <- 1:nrow(finalData)
  for ( i in seq_along(para.names) ) {
    cellID <- intersect(cellID, which(!is.na(match(as.character(finalData[,indexVar[i]]), para.codes[i]))))
  }
  if ( length(cellID) != 1) {
    stop("cellID:: check argument 'input' --> 0 or > 1 cells identified!\n")
  }
  return(cellID)
})
