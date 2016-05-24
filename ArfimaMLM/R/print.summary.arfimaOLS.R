#' @export
print.summary.arfimaOLS <-
function(x, ...){
  if(!is.null(x$ecm)){
    cat("\n###################################\n")
    cat("Summary Error Correction Mechanism: \n")
    print(x$ecm)
  }
  
  cat("\n###################################\n")
  cat("Fractional Differencing Parameters: \n\n")
  print(x$d)
  
  if(!is.null(x$arma)){
      cat("\n\n########################\n")
      cat("Summary AR/MA Estimates: \n\n")
      print(x$arma,digits=3)
  }
  
  cat("\n\n##################\n")
  cat("Summary OLS Model: \n")
  print(x$result)
}
