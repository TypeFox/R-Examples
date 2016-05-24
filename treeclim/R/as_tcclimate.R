##' Create the correct data format for variable selection
##' 
##' Check the climate data that is either supplied as a matrix, a
##' data.frame, or a list of data.frames. Reformat the data for
##' further usage into the internal format.
##' @param x the supplied climate data
##' @param varnames optionally supplied variable names
##' @return a data.frame with years in rows and monthly observations
##' in columns
##' @keywords manip, internal
as_tcclimate <- function(x, varnames = NULL) {
  
  msg1 <- "Format of climate data was not recognized. It is absolutely necessary that only complete years (months 1-12) are provided."
  
  ## is it a list?
  if (any(class(x) == "list")) {       
    ## handle list case
    n <- length(x)
    minyrs <- maxyrs <- numeric(n)
    ## loop through list members, and get min and max years for later 
    ## reformatting.
    for (i in 1:n) {
      ## shortcut for current list member
      y <- x[[i]]                      
      ## explanation see non-list case
      if (dim(y)[2] == 13) {
        perf_seq <- seq(y[1,1], y[dim(y)[1],1], 1)
        if (length(y[,1]) != length(perf_seq)) {
          stop(msg1)
        }
        if (!any(y[,1] == perf_seq)) {
          stop(msg1)
        } else {
          minyrs[i] <- min(y[,1])
          maxyrs[i] <- max(y[,1])
        }
      }
    }
    
    yrs <- max(minyrs):min(maxyrs)
    nyrs <- length(yrs)
    output_matrix <- matrix(NA, ncol = n + 2, nrow = nyrs*12)
    output_matrix[,1] <- rep(yrs, each = 12)
    output_matrix[,2] <- rep(1:12, nyrs)
    ## loop through list again, and put everything in place in the new output
    ## matrix
    for (i in 1:n) {
      y <- x[[i]]
      for (j in 1:nyrs) {
        if (any(y[,1] == yrs[j])) {
          ## write elements of specific line into i+2 th row of output_matrix
          output_matrix[which(output_matrix[,1] == yrs[j]), 2+i] <-
            unlist(y[which(y[,1] == yrs[j]), 2:13]) 
        }                               
      }
    }

    ## handle non-list case
  } else {
    ## should have 12 months columns and one year column
    if (dim(x)[2] == 13) {              
      ## check if the first column is perfect sequence of integer years. if
      ## expression evaluates to FALSE, then this is the case, if
      ## TRUE: stop.
      perf_seq <- seq(x[1,1], x[dim(x)[1],1], 1)
      if (length(x[,1]) != length(perf_seq)) {
        stop(msg1)
      }
      if (!any(x[,1] == perf_seq)) { 
        stop(msg1)
      } else {                          
        ## this is most probably a dendroclim-formatted set of climate data
        yrs <- unique(x[,1])
        nyrs <- length(yrs)
        output_matrix <- matrix(NA, ncol = 3, nrow = nyrs*12)
        output_matrix[,1] <- rep(yrs, each = 12)
        output_matrix[,2] <- rep(1:12, nyrs)
        for (i in 1:nyrs) {
          ## loop through years and write respective rows in respective columns 
          ## in output_matrix
          output_matrix[which(output_matrix[,1] == yrs[i]), 3] <-
            unlist(x[which(x[,1] == yrs[i]), 2:13])
        }
      }
      ## could still be the originally intended format of data.
    } else {                            
      ## check if the first column is a perfect sequence of integer years, each
      ## repeated 12 times. if expression evaluates to FALSE, then this is the
      ## case, else stop.
      perf_seq <- rep(x[1,1]:x[dim(x)[1],1], each = 12)
      if (length(x[,1]) != length(perf_seq)) {
        stop(msg1)
      }
      if (!any(x[,1] == perf_seq)) {
        stop(msg1)
      } else {
        if (!(any(x[,2] == rep(1:12, length(unique(x[,1])))))) {
          ## check if the second column is perfect sequence of 1:12 as often as
          ## there are individual years in column 1. if expression evaluates to
          ## FALSE, then this is the case, else stop.
          stop(msg1)
        } else {
          ## pass data on directly
          output_matrix <- x
        }
      }
    } 
  }

  output <- data.frame(output_matrix)
  
  ## do we have names?
  if (!is.null(varnames)) {
    if (length(varnames) == dim(output[2])) {
      colnames(output)[-c(1,2)] <- varnames
    } else {
      stop("`varnames` has to be of the same length as the number of parameters.")
    }
  }

  class(output) <- c("tc_climate", "data.frame")
  output
}
