summary.NTA  = function(object, pless=0.05,...){

  if (!inherits(object, "NTA")) stop(" is not of class NTA")

  cmat = cbind( object$PerVarImp,
                object$pvalue)
  colnames(cmat) = c( "CV-PerVarImp", "p-value" )
  rownames(cmat) =  dimnames(object$PerVarImp)[[1]]

  res = list(call = object$call,
             cmat = cmat,
             pless=pless)
  class(res) = "summary.NTA"

  return(res)
}

print.summary.NTA = function(x, ...){
  cat("Call:\n")
  print(x$call)

  w05=which(x$cmat[, 2] < x$pless)
  if(length(w05)>0){
    if(length(w05)>1){
      cat("\n")
      cat("\n  p-values less than ",x$pless,":")
      cat("\n ---------------------------")
      cat("\n")
      printCoefmat(x$cmat[w05,],has.Pvalue = TRUE)
    } else {
      cat("\n")
      cat("\n  p-values less than ",x$pless,":")
      cat("\n ---------------------------")
      cat("\n")
      printCoefmat(matrix(x$cmat[w05,],ncol = 2,
                          dimnames = list(rownames(x$cmat[c(w05,w05+1),])[1],c("VarImp","p-value"))),
                   has.Pvalue = TRUE)

    }
  } else {
    cat("\n")
    cat("\n  p-values less than ",x$pless,":")
    cat("\n ---------------------------")
    cat("\n")
    cat("None!!!")
  }
}
