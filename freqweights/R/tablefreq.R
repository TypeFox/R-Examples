## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2015-06-21 09:20 emilio on emilio-despacho>
## ============================================================




##' Create a table of frequencies
##'
##' Based on the \code{\link[dplyr]{count}} function,
##' it can also  work with matrices or external data bases and the result may be updated.
##'
##' It creates a frequency table of the \code{data}, or just of the columns specified in \code{vars}.
##'
##' If you provide a \code{freq} formula, the cases are weighted by the result of the formula. Any variables in the  formula are removed from the data set. If the data set is a matrix, the \code{freq} formula is a classic R formula. Otherwise, the expresion of \code{freq} is treated as a mathematical expression.
##'
##' This function uses all the power of \code{\link[dplyr]{dplyr}} to create frequency tables. The main advantage of this function  is that it works with on-disk data stored in data bases, whereas \code{\link[plyr]{count}} only works with in-memory data sets.
##'
##' In general, in order to use the functions of this package, the frequency table obtained by this function should fit in memory. Otherwise you must use the 'chunk' versions (\code{link{clarachunk}}, \code{link{biglmfreq}}).
##'
##' The code of this function are adapted from a wish list of the devel page of  \code{\link[dplyr]{dplyr}} (See references). Prof. Wickham also provides a nice introduction about how to use it with databases.
##' @title Create a table of frequencies
##' @param tbl an object that can be coerced to a \code{\link[dplyr]{tbl}}.  It must contain all variables in \code{vars} and in \code{freq}
##' @param vars variables to count unique values of. It may be a character vector
##' @param freq a name of a variable of the tbl object specifying frequency weights. See Details
##' @return A  \code{\link[dplyr]{tbl}} object a with label and freq columns. When it is possible, the last column is named \code{freq} and it represents the frequency counts of the cases. This object of class \code{tablefreq}, has two attributes:
##' \item{freq}{the weighting variable used to create the frequency table}
##' \item{colweights}{Name of the column with the weighting counts}
##' @note The author would like to thank Prof. Hadley Wickham who allowed the reutilisation of part of his code.
##' When using the update function, be careful with non-integer weights: The precision of the final weights may be wrong due to the multiple sums.
##' @seealso \code{\link[plyr]{count}}, \code{\link[dplyr]{tbl}}
##' @keywords manip
##' @import dplyr
##' @export
##' @references
##' Hadley Wickham. Count function \url{https://github.com/hadley/dplyr/issues/358}
##' Hadley Wickham. Databases \url{http://cran.rstudio.com/web/packages/dplyr/vignettes/databases.html}
##' @examples
##' tablefreq(iris)
##' tablefreq(iris, c("Sepal.Length","Species"))
##' a <- tablefreq(iris,freq="Sepal.Length")
##' tablefreq(a, freq="Sepal.Width")
##'
##' library(dplyr)
##' iris %>% tablefreq("Species")
##' 
##' tfq <- tablefreq(iris[,c(1:2)])
##'
##' chunk1 <- iris[1:10,c(1:2)]
##' chunk2 <- iris[c(11:20),]
##' chunk3 <- iris[-c(1:20),]
##' a <- tablefreq(chunk1)
##' a <- update(a,chunk2)
##' a <- update(a,chunk3)
##' a
##'
##' \dontrun{
##' 
##' ## External databases
##' library(dplyr)
##' if(require(RSQLite)){
##'   hflights_sqlite <- tbl(hflights_sqlite(), "hflights")
##'   hflights_sqlite
##'   tbl_vars(hflights_sqlite)
##'   tablefreq(hflights_sqlite,vars=c("Year","Month"),freq="DayofMonth")
##' }
##'
##' ##
##' ## Graphs
##' ##
##' if(require(ggplot2) && require(hflights)){
##'   library(dplyr)
##'
##'   ## One variable
##'   ## Bar plot
##'   tt <- as.data.frame(tablefreq(hflights[,"ArrDelay"]))
##'   p <- ggplot() + geom_bar(aes(x=x, y=freq), data=tt, stat="identity")
##'   print(p)
##'
##'   ## Histogram
##'   p <- ggplot() + geom_histogram(aes(x=x, weight= freq), data = tt)
##'   print(p)
##'
##'   ## Density
##'   tt <- tt[complete.cases(tt),] ## remove missing values
##'   tt$w <- tt$freq / sum(tt$freq) ## weights must sum 1
##'   p <- ggplot() + geom_density(aes(x=x, weight= w), data = tt)
##'   print(p)
##'
##'   ##
##'   ## Two distributions
##'   ##
##'   ## A numeric and a factor variable
##'   td <- tablefreq(hflights[,c("TaxiIn","Origin")])
##'   td <- td[complete.cases(td),]
##'
##'   ## Bar plot
##'   p <- ggplot() + geom_bar(aes(x=TaxiIn, weight= freq, colour = Origin),
##'                            data = td, position ="dodge")
##'   print(p)
##'
##'   ## Density
##'   ## compute the relative frequencies for each group
##'   td <- td %>% group_by(Origin) %>%
##'                mutate( ngroup= sum(freq), wgroup= freq/ngroup)
##'   p <- ggplot() + geom_density(aes(x=TaxiIn, weight=wgroup, colour = Origin),
##'                                data = td)
##'   print(p)
##'
##'   ## For each group, plot its values
##'   p <- ggplot() + geom_point(aes(x=Origin, y=TaxiIn, size=freq),
##'                              data = td, alpha= 0.6)
##'   print(p)
##'
##'   ## Two numeric variables
##'   tc <- tablefreq(hflights[,c("TaxiIn","TaxiOut")])
##'   tc <- tc[complete.cases(tc),]
##'   p <- ggplot() + geom_point(aes(x=TaxiIn, y=TaxiOut, size=freq),
##'                              data = tc, color = "red", alpha=0.5)
##'   print(p)
##'
##'   ## Two factors
##'   tf <- tablefreq(hflights[,c("UniqueCarrier","Origin")])
##'   tf <- tf[complete.cases(tf),]
##'
##'   ## Bar plot
##'   p <- ggplot() + geom_bar(aes(x=Origin, fill=UniqueCarrier, weight= freq),
##'                            data = tf)
##'   print(p)
##' }
##' }
tablefreq <- function(tbl, vars=NULL, freq=NULL){
  if( !inherits(tbl,"tbl") ) {
    if(inherits(tbl,"data.frame")) {
      tbl <- tbl_df(tbl)
    } else {
      tbl <- tbl_df(as.data.frame(tbl))
    }
  }
  if(is.null(vars)) {
    vars <- tbl_vars(tbl)
  }
  if(!is.null(freq)) {
    if(!freq %in% tbl_vars(tbl)) {
      stop("tablefreq: freq variable not in tbl")
    }
    vars <- setdiff(vars,freq)
  }
  if(is.null(vars)) {
    stop("tablefreq: null vars")
  } else {
    ##tbl <- ungroup(tbl)
##    tbl <- tbl %>% ungroup() %>% evaldp(group_by, vars)
       tbl <- tbl  %>% evaldp(group_by, vars)
  }
  
  ## print(nmsfreq)
  ##print(grouped)
  if(is.null(freq)) {
    dots <- "n()"
  } else{
    tbl <- evaldp(tbl, filter, paste(freq, ">0"))
    ##dots <- paste( "sum(as.numeric(", freq,"))")
    dots <- paste( "sum(", freq,")")
   }
  ## Name of the last column
  nms <- c(vars,"freq")
  nmsfreq <- make.names(nms,unique=TRUE)[length(nms)]
  dots <- paste(nmsfreq, " = ", dots)
  ##  print("count: dots")
  ##  print(dots)
  x <- evaldp(tbl, summarise, dots) %>% ungroup() 
  vars <- tbl_vars(x)
  ##rownames(x) <- NULL
  attr(x,"freq") <- freq
  attr(x,"colweights") <- vars[length(vars)]
  class(x) <- c(class(x),"tablefreq")
  return(x)
}



## ============================================================
##
## ============================================================





## ============================================================
##
## ============================================================

##' @param object a \code{tablefreq} object
##' @param ... more data
##' @method update tablefreq
##' @import dplyr
##' @export
##' @rdname tablefreq
update.tablefreq <- function(object, ...) {
  if (!inherits(object, "tablefreq"))
    stop("update.table: error")
  dots <- list(...)
  moredata <- dots[[1]]
  if(is.null(moredata)) {
    vars <- tbl_vars(object) 
    lastcolumn <- vars[length(vars)]
    return(tablefreq(object,freq=lastcolumn))
  }
  varsobj <- colnames(object)[colnames(object) != attr(object,"colweights")]
  x  <- tablefreq(moredata, vars=varsobj, freq=attr(object,"freq"))
  ## Last column is the weighting variable. Assure that
  ## the name is the same in both data sets.
  if( any(!colnames(object) %in% colnames(x)) ) {
    stop("update. Different colnames")
  }
  ## Now, join both data sets and calculate a new table of frequencies
  x <- tablefreq(rbind_list(object, x[, colnames(object)]),
                 freq=attr(object,"colweights"))
  ## Preserve the original formula freq
  attr(x, "freq") <- attr(object,"freq")
  x
}
