"dat2bernquaf" <- function(x, data, interval=NA, ...) {

   if(length(x) > 1) {
     warning("length of x is > 1, only first value will be used")
     x <- x[1]
   }
   if(length(data) < 2) {
     warning("length of data is < 2, returning NA")
     return(NA)
   }

   afunc <- function(f, ...) {
      xFbern <- dat2bernqua(f, data, ...)
      return(x - xFbern)
   }

   n <- length(data)
   if(length(interval) == 0 || is.na(interval)) {
      interval <- c( 1/(n+1), 1 - 1/(n+1) )
   } else if(length(interval) == 1) {
      interval <- c( interval, 1 - interval )
   } else if(length(interval) == 2) {
      # Do nothing
   } else {
      warning("Invalid logic leading to the development of the interval for the root solver, returning NA: interval=",interval)
      return(NA)
   }
   if(interval[2] < interval[1]) {
      warning("Backwards interval, sorting the two values and continuing on")
      interval <- sort(interval)
   }
   if(interval[1] <= 0) {
      warning("Left of interval <= 0, returning NA")
      return(NA)
   }
   if(interval[2] >= 1) {
      warning("Right of interval >= 1, returning NA")
      return(NA)
   }


   rt  <- NULL
   try(rt <- uniroot(afunc, interval, ...))
   if(is.null(rt)) {
     warning("Unable to root")
     return(NA)
   }
   zz <- list(x=x, f=rt$root, interval=interval,
              afunc.root=rt$f.root, iter=rt$iter, estim.prec=rt$estim.prec,
              source="dat2bernquaf")
   return(zz)
}


