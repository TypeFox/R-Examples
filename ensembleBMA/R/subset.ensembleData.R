`[.ensembleData` <-
function (x, i, j) 
{ 
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
    ncolx <- ncol(x)
    matchCall <- match.call()
    matchCall[[1]] <- as.name("[.data.frame")
    if (missing(i)) matchCall$i <- 1:nrow(x)
    nForcs <- nForecasts <- ensembleSize(x)
    exchangeable <- attr(x, "exchangeable")
    forecastHour <- attr(x, "forecastHour")
    initializationTime <- attr(x, "initializationTime")
    startupSpeed <- attr(x, "startupSpeed")
    if (!missing(j) && !is.null(j)) {
      if (is.logical(j)) {
        if (length(j) != nForecasts)
          stop("logical index must refer to the forecasts")
        j <- (1:nForecasts)[j]
        nForcs <- length(j)
        exchangeable <- exchangeable[j]
      }
      else if (is.character(j)) {
        m <- match(j, names(x)[1:nForecasts], nomatch = 0)
        if (any(!m))
          stop("character index must refer to the forecasts")
        if (any(duplicated(j))) 
          stop("repeated indexes not allowed")
        nForcs <- length(j)
        I <- 1:nForecasts
        names(I) <- ensembleMembers(x)
        j <- I[j]
        names(j) <- NULL
        exchangeable <- exchangeable[j]
      }
      else {
        if (any(abs(j) > nForcs)) 
           stop("column index must be confined to the forecasts")
        if (any(duplicated(j))) 
          stop("repeated indexes not allowed")
        j <- (1:nForecasts)[j]
        nForcs <- length(j)
        exchangeable <- exchangeable[j]
      }
      if (nForcs < ncolx) {
        matchCall$j <- c(j, (nForecasts+1):ncolx)
      }
      else matchCall$j <- j
    }
    else matchCall$j <- 1:ncolx
    if (!missing(i)) {
      v <- (1:nrow(x))
      names(v) <- dimnames(x)[[1]]
      i <- v[i]
      names(i) <- NULL
      if (any(duplicated(i))) 
        stop("repeated entries not allowed")
    }
    matchCall$drop <- FALSE
    listCall <- as.list(matchCall)
    nam <- names(listCall)
    listCall <- listCall[c(1,2,which(nam == "i"), which(nam == "j"), 
                  length(listCall))]
    names(listCall) <- NULL
    x <- eval(as.call(listCall), parent.frame())
    attr(x, "initializationTime") <- initializationTime
    attr(x, "forecastHour") <- forecastHour
    attr(x, "startupSpeed") <- startupSpeed
    attr(x, "exchangeable") <- exchangeable 
    attr(x, "ensembleSize") <- nForcs
    x
}

