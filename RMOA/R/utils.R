
#' Convert character strings to factors in a dataset
#'
#' Convert character strings to factors in a dataset
#'
#' @param x object of class data.frame
#' @param ... other parameters currently not used yet
#' @return a data.frame with the information in \code{x} where character columns are converted to factors
#' @export 
#' @examples
#' data(iris)
#' str(iris)
#' mydata <- factorise(iris)
#' str(mydata)
factorise <- function(x, ...){
  UseMethod(generic="factorise", object=x)
}

##' @S3method factorise data.frame
factorise.data.frame <- function(x, ...){  
  if(inherits(x, "data.frame")){
    tofactor <- sapply(x, FUN=function(x) inherits(x, "character"))
    tofactor <- tofactor[tofactor == TRUE]
    if(length(tofactor) == 0){
      return(x)
    }
  }
  x <- unclass(x)
  x <- as.data.frame(x, stringsAsFactors=TRUE)
  x
}



recoder <- function (x, from = c(), to = c()) {
  missing.levels <- unique(x)
  missing.levels <- missing.levels[!missing.levels %in% from]
  if (length(missing.levels) > 0) {
    from <- append(x = from, values = missing.levels)
    to <- append(x = to, values = missing.levels)
  }
  to[match(x, from)]
}
