#' @rdname crsm
#'
#' @export
#'
#' @method summary CRSM
#'
#' @param object object of class \code{CRSM}
#' @param \dots \dots{}


summary.CRSM <-
function(object, ...){

  cat("\n Call: ", deparse(object$call), "\n\n")

  cat("\n Number of iterations: \n")

  cat("\t Maximum number of iterations: \t", max(object$iterations), "\n", "\tMinimum number of iterations: \t", min(object$iterations), "\n\n\n")


    parall <- rbind(cbind("item estimates"=object$itempar, "SE"=object$itempar_se, "SE low"= object$itempar_se_low, "SE up"=object$itempar_se_up),"disppar"=c(object$disppar, object$disppar_se, object$disppar_se_low, object$disppar_se_up) )

  cat("Parameter estimates: \n")
  print(parall)

}
