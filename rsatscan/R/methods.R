#' Methods for satscan-class objects
#' 
#' These functions define the default methods for satscan-class objects, which
#' are the result objects from a call to satscan()
#' 
#' @param x is a satscan object
#' @param ... vestigial, for compatability with the default summary method
#' 
#' @return x, invisibly.  Side effect is to display ss$main, the SaTScan text report
#' 
#' @export

print.satscan = function(x, ...) {
   stopifnot(class(x)=="satscan")
   cat(x$main,fill=1)
   invisible(x)
}


#' Methods for satscan-class objects
#' 
#' These functions define the default methods for satscan-class objects, which
#' are the result objects from a call to satscan()
#' 
#' @param object is a satscan object
#' @param ... vestigial, for compatability with the default summary method
#' 
#' @return object, invisibly.  Side effect is to display minimal facts contained in ss
#' 
#' @export


summary.satscan = function(object, ...) {
  stopifnot(class(object)=="satscan")
  cat(object$main[c(9:11, 15:20,22)], fill=1)
  cat(paste0("There were ", dim(object$col)[1], " clusters identified."), fill=1)
  cat(paste0("There were ", sum(object$col$P_VALUE < .05), " clusters with p < .05."))
  invisible(object)
}



