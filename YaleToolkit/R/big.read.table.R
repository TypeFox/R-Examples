#' @title Read in chunks from a large file with row/column filtering
#' to obtain a reasonable-sized data.frame. 
#' @param file the name of the file, obviously
#' @param nrows the chunk size; consider reducing this if there are
#' lots of columns
#' @param sep by default we expect a CSV file
#' @param header is \code{TRUE} by default
#' @param row.names I really dislike row names
#' @param cols for filtering column by name or number (supporting negative indexing)
#' @param rowfilter a function that is assumed to take a chunk as a
#' data frame and return a smaller data frame (with fewer rows), separately
#' from the column filtering.
#' @param as.is \code{TRUE} by default
#' @param estimate do a preliminary estimation of the work to be done,
#' and then have a chance to bail out if it looks like a bad idea
#' @note This is very much 'in development' and could be buggy.  I put it here
#' as I used some example in one of my courses, but then I needed to update
#' the package to keep CRAN happy.  So here it is.  Buyer Beware.  - Jay
#' @examples
#' data(CO2)
#' write.csv(CO2, "CO2.csv", row.names=FALSE)
#' x <- big.read.table("CO2.csv", nrows=10)
#' head(x)
#' @export
big.read.table <- function(file, nrows=100000, sep=",",
                           header=TRUE, row.names=NULL,
                           cols=NULL, rowfilter=NULL,
                           as.is=TRUE, estimate=FALSE)
{  
  if (estimate) {
    nlines <- getnrows(file)
    x <- read.table(file, sep=sep, row.names=row.names,
                    nrows=min(nlines, 1000), header=header)
    if (!is.null(cols)) x <- x[,cols,drop=FALSE]
    cat("Estimated read size without row filtering:",
        floor(object.size(x)*nlines/nrow(x)/1e6), "MB\n")
    if (interactive()) {
      ANSWER <- readline("Continue with read (Y/n)? ")
      if (substring(ANSWER, 1, 1) != "Y") {
        warning("Terminated read.")
        return(NULL)
      }
    }
  }
  
  iter <- iread.table(file, header=header,
                      row.names=row.names, sep=sep,
                      nrows=nrows, as.is=as.is)
  
  ans <- foreach(x=iter, .combine=rbind) %do% {
    if (!is.null(rowfilter)) x <- rowfilter(x)
    if (!is.null(cols)) x <- x[,cols,drop=FALSE]
    gc()
    return(x)
  }
  
  return(ans)
}