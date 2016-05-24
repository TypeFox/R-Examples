##' Calculate the difference between two date in terms of unit of time.
##'
##' This function calculate the difference between two date in terms of unit of time, and age can be obtained when \sQuote{year} is specified as unit argument.
##' @param birthday a character or character vector specifying birthday or base date.
##' @param targetdate a character or character vector specifying target date.
##' @param unit a character specifying unit for calculating the difference between the two dates. Values of "year", "month" and "day" are supported.
##' @return a vector of age
##' @examples calcAge("1963-2-3")
##' @examples calcAge("1970-1-1", unit="day")
##' @export
calcAge <- function(birthday, targetdate=Sys.Date(), unit="year"){
  if(length(targetdate) == 1){
    sapply(birthday,
           function(x){
             tryCatch(
               {length(seq(as.Date(x), as.Date(targetdate), unit)) - 1},
               error=function(e){NA})
           },
           USE.NAMES=FALSE
           )
  } else {
    n <- length(birthday)
    # recycle rule
    targetdate <- rep(targetdate, length.out=n)[1:n]
    mapply(function(x, y){
      tryCatch(
               {length(seq(as.Date(x), as.Date(y), unit)) - 1},
               error=function(e){NA})
    },
           birthday, targetdate)
  }
}
