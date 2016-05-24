# Method for accessing the R-squares for the submodels of a path model fitted by 'sempls'
rSquared <- function(object, ...){
  UseMethod("rSquared")
}

rSquared.sempls <- function(object, na.rm=FALSE, ...){
  #Y_hat <- object$factor_scores %*% object$path_coefficients
  Y_hat <- predict(object, what="LVs", scale="scaled", ...)
  if(sum(is.na(Y_hat)) > 0 & !na.rm) stop("Use argument 'na.rm=TRUE'!")
  R_squared <- apply(Y_hat, 2, var, na.rm=na.rm) / apply(object$factor_scores, 2, var, na.rm=na.rm)
  R_squared[R_squared==0] <- NA
  R_squared <- as.matrix(R_squared)
  colnames(R_squared) <- "R-squared"
  class(R_squared) <- "rSquared"
  return(R_squared)
}

print.rSquared <- function(x, na.print=".", digits=2, ...){
  print.table(x, na.print=na.print, digits=digits, ...)
  invisible(x)
}

