summary.PimpTest = function(object, pless=0.05,...){

  if (!inherits(object, "PimpTest")) stop(" is not of class PimpTest")

  if(object$para){
      cmat = cbind( object$meanPerVarImp,
                    object$sdPerVarImp,
                    object$p.ks.test)
      colnames(cmat) = c( "mean(PerVarImp)", "sd(PerVarImp)", "ks.p-value" )
      rownames(cmat) =  dimnames(object$VarImp)[[1]]
  } else{
    cmat = NULL
  }


  cmat2 = cbind(object$VarImp ,
               object$pvalue)
  colnames(cmat2) = c("VarImp","p-value")
  rownames(cmat2) =  dimnames(object$VarImp)[[1]]

  res = list(para = object$para,
             type = object$type,
             call = object$call,
             call.PIMP=object$call.PIMP,
             cmat = cmat,
             cmat2 = cmat2,
             pless = pless)

  class(res) <- "summary.PimpTest"
  return(res)

}


print.summary.PimpTest = function(x, ...){


  #######################################################################
  # Output of results
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$call.PIMP)
  cat("type:  ")
  print(x$type)
  cat("\n")
  if(x$para){
      cat("\n Kolmogorov-Smirnov test for the null importances:")
      cat("\n -------------------------------------------------")
      cat("\n")
      printCoefmat(x$cmat,has.Pvalue = T)

  }


  w05=which(x$cmat2[, 2] < x$pless)
  if(length(w05)>0){
    if(length(w05)>1){
      cat("\n")
      cat("\n  p-values less than ",x$pless,":")
      cat("\n ----------------------")
      cat("\n")
      printCoefmat(x$cmat2[w05,],has.Pvalue = TRUE)
    } else {
      cat("\n")
      cat("\n  p-values less than ",x$pless,":")
      cat("\n ----------------------")
      cat("\n")
      printCoefmat(matrix(x$cmat2[w05,],ncol = 2,
                          dimnames = list(rownames(x$cmat2[c(w05,w05+1),])[1],c("VarImp","p-value"))),
                   has.Pvalue = TRUE)

    }
  } else {
    cat("\n")
    cat("\n  p-values less than ",x$pless,":")
    cat("\n ----------------------")
    cat("\n")
    cat("None!!!")
  }
}
