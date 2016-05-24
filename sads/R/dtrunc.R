dtrunc <- function(f, x, trunc, coef, log = FALSE){
  pf <- get(paste("p", f, sep=""), mode = "function")
  df <- get(paste("d", f, sep=""), mode = "function")
  tt <- rep(0, length(x))
  if (!missing(trunc)){
    tt[x > trunc] <- do.call(df, c(list(x = x[x>trunc]), coef))/(1 - do.call(pf, c(list(q = trunc), coef)))    
  } else{
    tt <- do.call(df, c(list(x = x), coef))
  }
  if (log) tt <- log(tt)
  return(tt)
}
