as.Date1970 <- function(x, ...){
#
  if(is.numeric(x)){
      dots <- list(...)
      dots$x <- x
      if(!('origin' %in% names(dots))){
          o1970 <- as.Date('1970-01-01', '%Y-%m-%d')
          dots$origin <- o1970
      }
      xD <- do.call(as.Date, dots)
      return(xD)
  }
  as.Date(x, ...)
}
