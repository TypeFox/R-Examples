########################################
### methods only for class 'dataObj' ###
########################################
#' @aliases get.dataObj,dataObj,character-method
#' @rdname get.dataObj-method
setMethod(f="get.dataObj", signature=c("dataObj", "character"),
  definition=function(object, type) {
    if ( !type %in% c("rawData", "dimVarInd", "freqVarInd",
        "numVarInd", "weightVarInd", "sampWeightInd",
        "isMicroData", "numVarNames", "freqVarName", "varName") ) {
      stop("get.dataObj:: argument 'type' is not valid!\n")
    }
    if ( type == "rawData" ) {
      return(g_raw_data(object))
    }
    if ( type == "dimVarInd" ) {
      return(g_dimvar_ind(object))
    }
    if ( type == "freqVarInd" ) {
      return(g_freqvar_ind(object))
    }
    if ( type == "numVarInd" ) {
      return(g_numvar_ind(object))
    }
    if ( type == "weightVarInd" ) {
      return(g_weightvar_ind(object))
    }
    if ( type == "sampWeightInd" ) {
      return(g_sampweight_ind(object))
    }
    if ( type == "isMicroData" ) {
      return(g_is_microdata(object))
    }
    if ( type == 'numVarNames' ) {
      return(names(g_raw_data(object))[g_numvar_ind(object)])
    }
    if ( type == 'freqVarName' ) {
      return(names(g_raw_data(object))[g_freqvar_ind(object)])
    }
    if ( type == 'varName' ) {
      return(names(g_raw_data(object))[g_dimvar_ind(object)])
    }
  }
)

#' @aliases set.dataObj,dataObj,character,listOrNULL-method
#' @rdname set.dataObj-method
setMethod(f='set.dataObj', signature=c("dataObj", "character", "listOrNULL"),
  definition=function(object, type, input) {
    if ( !type %in% c("rawData") ) {
      stop("set.dataObj:: check argument 'type'!\n")
    }
    if ( type == "rawData" ) {
      s_raw_data(object) <- input[[1]]
    }
    validObject(object)
    return(object)
  }
)

#' @aliases init.dataObj,list-method
#' @rdname init.dataObj-method
setMethod(f='init.dataObj', signature=c('list'),
  definition=function(input) {
    freq <- N <- .N <- NULL
    inputData <- input$inputData
    dimVarInd <- input$dimVarInd
    freqVarInd <- input$freqVarInd
    numVarInd <- input$numVarInd
    weightInd <- input$weightInd
    sampWeightInd <- input$sampWeightInd
    wExists <- FALSE
    if ( !is.null(weightInd)) {
      wExists <- ifelse(weightInd==freqVarInd | weightInd %in% numVarInd, TRUE, FALSE)
    }
    isMicroData <- FALSE

    ## aggregate data, use data.table
    datO <- data.table(inputData, key=colnames(inputData)[dimVarInd])
    rawData <- datO[, list(N=.N), by=key(datO)]
    N <- rawData$N
    rawData$N <- NULL

    if ( any(N != 1)  ) {
      isMicroData <- TRUE
      rawData <- copy(datO)
      rawData <- rawData[,key(rawData),with=FALSE]

      if ( is.null(freqVarInd) ) {
        set(rawData, NULL, "freq", 1)
        freqVarInd <- which(colnames(rawData)=="freq")
      } else {
        f <- datO[,get(colnames(datO)[freqVarInd]),]
        set(rawData, NULL, "freq", as.numeric(f))
      }
    } else {
      # data already aggregated
      if ( is.null(freqVarInd) ) {
        set(rawData, NULL, "freq", as.numeric(N))
        freqVarInd <- which(colnames(rawData)=="freq")
      } else {
        set(rawData, NULL, "freq", as.numeric(datO[[freqVarInd]]))
      }
    }
    freqVarInd <- which(colnames(rawData)=="freq")
    dimVarInd <- match(key(rawData), colnames(rawData))

    if ( !is.null(sampWeightInd) ) {
      if ( isMicroData == TRUE ) {
        sw <- datO[,get(colnames(datO)[sampWeightInd])]
      } else {
        sw <- datO[,list(sw=sum(get(colnames(datO)[sampWeightInd]))), by=key(datO)]$sw
      }
      set(rawData, NULL, colnames(datO)[sampWeightInd], sw)

      sampWeightInd <- ncol(rawData)
      set(rawData, NULL, "freq", as.numeric(sw*rawData$freq))
    }

    ## numvars
    nr_keys <- length(key(datO))
    if ( !is.null(numVarInd) ) {
      c.start <- ncol(rawData)
      cols <- colnames(datO)[numVarInd]

      if ( isMicroData == TRUE ) {
        for ( j in 1:length(cols)) {
          set(rawData, NULL, cols[j], as.numeric(datO[[cols[j]]]))
        }
      } else {
        xx <- datO[,lapply(.SD, sum), by=key(datO), .SDcols=cols]
        rawData <- merge(rawData, xx)
      }
      numVarInd <- (c.start+1):ncol(rawData)
    }

    ## weight var
    if ( !is.null(weightInd) ) {
      if ( isMicroData == TRUE ) {
        w <- datO[,get(colnames(datO)[weightInd])]
      } else {
        w <- datO[,list(w=sum(get(colnames(datO)[weightInd]))), by=key(datO)]$w
      }
      # we set this column only, if weightInd does not already exist in rawData
      # this could be the case, if a variable also used as numVar is specified as weightVar
      if ( !wExists ) {
        set(rawData, NULL, colnames(datO)[weightInd], as.numeric(w))
        weightInd <- ncol(rawData)
      } else {
        weightInd <- match(colnames(datO)[weightInd], names(rawData))
        if ( is.na(weightInd)) {
          weightInd <- match("freq", names(rawData))
        }
      }
    }

    ## do not use factors
    cols <- colnames(rawData)[dimVarInd]
    rawData[,cols] <- rawData[, lapply(.SD, as.character), .SDcols=cols]
    rm(datO)

    setkeyv(rawData, colnames(rawData)[dimVarInd])

    out <- new("dataObj",
      rawData=rawData,
      dimVarInd=dimVarInd,
      freqVarInd=freqVarInd,
      numVarInd=numVarInd,
      weightVarInd=weightInd,
      sampWeightInd=sampWeightInd,
      isMicroData=isMicroData
    )
    return(out)
  }
)

setMethod(f="g_raw_data", signature=c("dataObj"), definition=function(object) {
  return(object@rawData)
})

setMethod(f="g_dimvar_ind", signature=c("dataObj"), definition=function(object) {
  return(object@dimVarInd)
})

setMethod(f="g_freqvar_ind", signature=c("dataObj"), definition=function(object) {
  return(object@freqVarInd)
})

setMethod(f="g_numvar_ind", signature=c("dataObj"), definition=function(object) {
  return(object@numVarInd)
})

setMethod(f="g_weightvar_ind", signature=c("dataObj"), definition=function(object) {
  return(object@weightVarInd)
})

setMethod(f="g_sampweight_ind", signature=c("dataObj"), definition=function(object) {
  return(object@sampWeightInd)
})

setMethod(f="g_is_microdata", signature=c("dataObj"), definition=function(object) {
  return(object@isMicroData)
})

setMethod(f="g_numvar_names", signature=c("dataObj"), definition=function(object) {
  return(names(g_raw_data(object))[g_numvar_ind(object)])
})

setMethod(f="g_freqvar_name", signature=c("dataObj"), definition=function(object) {
  return(names(g_raw_data(object))[g_freqvar_ind(object)])
})

setMethod(f="g_var_name", signature=c("dataObj"), definition=function(object) {
  return(names(g_raw_data(object))[g_dimvar_ind(object)])
})

setReplaceMethod(f="s_raw_data", signature=c("dataObj", "list"), definition=function(object, value) {
  object@rawData <- value[[1]]
  return(object)
})










