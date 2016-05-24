ptrunc <- function(f, q, trunc, coef, lower.tail=TRUE, log.p=FALSE){
  tt <- q
  pf <- get(paste("p", f, sep = ""), mode = "function")
  if(!missing(trunc)){
    aa <- rep(trunc, length(q))
    tt <- do.call(pf, c(list(q = apply(cbind(q, aa), 1, max)), coef,lower.tail=lower.tail))
    if(lower.tail){
      tt <- tt - do.call(pf, c(list(q = aa), coef))
    }
    tt <- tt/(1 - do.call(pf, c(list(q = aa), coef)))
  }
  else{
    tt <- do.call(pf, c(list(q = q), coef,lower.tail=lower.tail))
  }
  if(log.p)return(log(tt)) else return(tt)
}
