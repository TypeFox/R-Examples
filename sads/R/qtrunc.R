qtrunc <- function(f, p, trunc, coef, lower.tail = TRUE, log.p = FALSE){
  if(log.p) p <- exp(p)
  tt <- p
  pf <- get(paste("p", f, sep = ""), mode = "function")
  qf <- get(paste("q", f, sep = ""), mode = "function")
  if(!missing(trunc)){
    if(!lower.tail){
      aa <- do.call(pf, c(list(q = trunc), coef, lower.tail = T)) +
          (1-p)*(1 - do.call(pf, c(list(q = trunc), coef, lower.tail = T)))
      tt <- do.call(qf, c(list(p = aa), coef, lower.tail = T))
    }
    else{
      aa <- do.call(pf, c(list(q = trunc), coef, lower.tail = lower.tail)) +
          p*(1 - do.call(pf, c(list(q = trunc), coef, lower.tail = lower.tail)))
      tt <- do.call(qf, c(list(p = aa), coef, lower.tail = lower.tail))
    }
  }
  else{
    tt <- do.call(qf, c(list(p = p), coef, lower.tail = lower.tail))
  }
  return(tt)
}
