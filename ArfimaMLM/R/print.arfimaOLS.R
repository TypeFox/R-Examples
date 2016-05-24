#' @export
print.arfimaOLS <-
function(x, ...){
  if(!is.null(x$ecm)){
    cat("\n###################################\n")
    cat("Results Error Correction Mechanism: \n")
    print(x$ecm)
  }
  cat("\n###################################\n")
  cat("Fractional Differencing Parameters: \n\n")
  print(x$d)
  if(!is.null(x$arma)){
      cat("\n\n################\n")
      cat("AR/MA Estimates: \n\n")
      print(x$arma)
  }
  cat("\n##################\n")
  cat("Results OLS Model: \n")
  print(x$result)
}
