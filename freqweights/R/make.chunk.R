## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2014-04-25 07:55 emilio on emilio-Satellite-P100>
## ============================================================

## ============================================================
##
## ============================================================



##' Fast and friendly chunk file finagler
##' 
##' Read a file chunk by chunk 
##' 
##' It creates a function that reads sucesive chunks of
##' the data referenced by \code{input} usings the
##' \code{\link[data.table]{fread}} function. The \code{input} is characterized
##' in the help page of \code{\link[data.table]{fread}}. The data contained in the
##' \code{input} reference should not have any header.
##' 
##' This function is inspired by the \code{\link[biglm]{bigglm}} example.
##' 
##' @param input a length 1 character string. See Details. 
##' @param FUN any function applicated to each chunk 
##' @param chunksize number of lines for each chunk 
##' 
##' @return A function with an logical argument, \code{reset}. If this argument
##' is \code{TRUE}, it indicates that the data should be reread from the
##' beginning by subsequent calls. When it reads all the data, it automatically
##' resets the file. This function returns the value of \code{FUN} applied to
##' the chunk. By default, the chunk is returned as a
##' \code{\link[dplyr]{tbl_df}} object.
##' @seealso \code{\link[biglm]{bigglm}}, \code{\link[data.table]{fread}},  \code{\link[dplyr]{tbl_df}} 
##' @keywords IO manip
##' @export
##' @import dplyr
##' @importFrom data.table fread
##' @importFrom plyr aaply adply alply daply ddply dlply laply ldply llply a_ply d_ply l_ply
##' @examples
##' 
##' \dontrun{
##' library(hflights)
##' nrow(hflights) # Number of rows
##' 
##' ## We create a file with no header
##' input <- "hflights.csv"
##' write.table(hflights,file=input,sep=",",
##'             row.names=FALSE,col.names=FALSE)
##' 
##' ## Get the number of rows of each chunk
##' readchunk <- make.readchunk(input,FUN=function(x){NROW(x)})
##' 
##' a <- NULL
##' while(!is.null(b <- readchunk())) {
##'   if(is.null(a)) {
##'     a <- b
##'   } else {
##'     a <- a+b
##'   }
##' }
##' all.equal(a, nrow(hflights))
##' 
##' ## It resets automatically the file 
##' a <- NULL 
##' while(!is.null(b <- readchunk())) {
##'   if(is.null(a)) {
##'     a <- b
##'   } else {
##'     a <- a+b
##'   }
##' }
##' all.equal(a, nrow(hflights))
##' }
make.readchunk <- function(input, FUN=identity, chunksize=5000L){
  ## Adapted from
  ## biglm::bigglm
  ##  ## https://github.com/RevolutionAnalytics/rhdfs/blob/master/pkg/inst/examples/biglm.integration.r
  ## http://stackoverflow.com/questions/19894194/reading-in-chunks-at-a-time-using-fread-in-package-data-table
  chunk <- 1L
  done <- FALSE
  readedline <- 0L
  nrecords <- 0L
  if(chunksize < 2 && !verifylastrecord) {
    stop("make.readchunk: chunksize too small")
  }
  verifylastrecord <- TRUE
  if(verifylastrecord) {
    r <- data.table::fread(input, select=1, showProgress=FALSE)
    nrecords <- NROW(r)
  }
  f <- function(reset=FALSE){
    if(reset){
      ## When reset=TRUE it indicates that the data should be reread from the
      ## beginning by subsequent calls.
      chunk <<- 1L
      done <<-  FALSE
      readedline <<- 0L
    } else{
      ## it returns a data frame with the
      ## next chunk of data or NULL if no more data are available
      ## print(paste("file", input," nreadedlines", readedline))
      if(!done) {
        r <- tbl_df(data.table::fread(input,
                                             skip = readedline,
                                             nrows = chunksize,
                                             showProgress=FALSE))
      } else {
        r <- NULL
      }
      ##print(paste("dimension r",dim(r)))
      if( (nrow(r)< 2 &&  readedline%%chunksize) || # the last chunk was not fullfilled, but it reads a new record
         (nrow(r)< 2 &&  !readedline%%chunksize && # the last chunk was fulfilled, and now it reads just one record
          (!verifylastrecord || # If we do not verify the last record
           (verifylastrecord && readedline != nrecords -1)))){ # or this last record does not fit with the total records of the data set
        ## It is a degenerated fread (fread function always read at leat one row)
        ## we forget this record
        ## And we reset the chunk
        chunk <<- 1L
        done <<-  FALSE
        readedline <<- 0L
        return(NULL)
      } else {
        chunk <<- chunk + 1
        readedline <<- readedline + nrow(r)
        ##print(paste("chunks que van leidas ", readedline  ))
        return(FUN(r))
      }
    }
  }
  return(f)
}


## ============================================================
##
## ============================================================

## ============================================================
##
## ============================================================




## ============================================================
##
## ============================================================




