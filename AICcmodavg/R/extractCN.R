##extract condition number
#Values of the condition number near 0 or negative indicate a problem, possibly 
#indicating fitting a model with too many parameters for the given data set.
extractCN <- function(mod, method = "svd", ...){
  UseMethod("extractCN", mod)
}


extractCN.default <- function(mod, method = "svd", ...) {
  stop("\nFunction not yet defined for this object class\n")
}


##unmarkedFit objects
extractCN.unmarkedFit <- function(mod, method = "svd", ...) {
  ##extract Hessian matrix
  hess <- mod@opt$hessian

  ##SVD
  if(identical(method, "svd")) {
    s <- svd(hess, nu = 0, nv = 0)$d
    CN <- max(s)/min(s[s > 0])
  }

  ##eigen
  if(identical(method, "eigen")) {
    eigenvals <- eigen(hess)$values
    CN <- max(eigenvals)/min(eigenvals)
  }

  ##compute log10
  logKappa <- log10(CN)
  
  ##arrange results
  out <- list("CN" = CN, "log10" = logKappa, "method" = method)

  class(out) <- "extractCN"
  return(out)
}



##print method
print.extractCN <- function(x, digits = 2, ...) {
  nice.vector <- c("Condition number" = x$CN,
                   "log10" = x$log10)
  print(round(nice.vector, digits = digits))
}
  
