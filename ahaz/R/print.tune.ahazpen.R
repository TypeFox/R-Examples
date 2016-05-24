"print.tune.ahazpen"<-function(x, digits = max(3, getOption("digits") - 3), ...){
    ## Purpose: Print 'tune.ahazpen' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   object: 'tune.ahazpen' object
    ##   digits: number of digits in output
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
  
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if(x$tune$type=="CV"){
    cat("Cross-validation: ",format(as.integer(x$nfolds), digits = digits)," folds \n")
  if(x$tune$rep>1)
    cat("Repetitions     : ",format(as.integer(x$nfolds), digits = digits),"\n\n")
  else
    cat("\n")
  } else {
    cat("BIC-type penalty parameter selection \n\n")
  }
  cat("Length of lambda sequence : ",format(length(x$lambda)),"\n")
  cat("Optimal lambda            : ",format(min(x$lambda.min), digits = digits),"\n")
  cat("d.f. at optimal lambda    : ",format(x$df[x$lambda==x$lambda.min][1], digits = digits),"\n\n")
}
