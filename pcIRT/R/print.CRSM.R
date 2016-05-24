#' @rdname crsm
#'
#' @export
#'
#' @method print CRSM
#'
#' @param x object of class \code{CRSM}

print.CRSM <-
function(x, ...){

  cat("\n Call: ", deparse(x$call), "\n\n")

    parall <- rbind(cbind("item estimates"=x$itempar, "SE"=x$se.item.mean), "dispp"=c(x$disppar, x$se.distr.mean))
      parall <- rbind(cbind("item estimates"=x$itempar, "SE"=x$itempar_se),"dispp"=c(x$disppar, x$disppar_se) )

  cat("Parameter estimates: \n")
  print(parall)

}
