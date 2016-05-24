#' @rdname mprm
#'
#' @export
#'
#' @method print MPRM
#'
#' @param x object of class \code{MPRM}

print.MPRM <-
function(x, ...){
  
  cat("\n Call: ", deparse(x$call), "\n\n")
  
  cat("Function calls: \n")
    print(x$fun_calls)
   
  cat("Parameter estimates: \n")
  print(x$itempar)
  print(x$itempar_se)
  
}
