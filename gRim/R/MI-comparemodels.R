

compareModels.iModel <- function(object, object2, k=2, ...){
  comp <- compareGC(object$glist, object2$glist)
  logL1 <- logLik(object)
  logL2 <- logLik(object2)
  
  if (!comp["comparable"]){
    cat("Models are not nested and can hence not be compared\n")
    return()
  } else {
    smaller <- which(comp[1:2])
    if (smaller==1){
      larger <- 2
      d.logL <- c(logL2-logL1)
      d.df   <- c(attr(logL1,"df")-attr(logL2, "df"))
      large  <- object2$glist
      small  <- object$glist
    } else {
      larger <- 1
      d.logL <- c(logL1-logL2)
      d.df   <- c(attr(logL2,"df")-attr(logL1, "df"))
      large  <- object$glist
      small  <- object2$glist
    }  
  }
      
  ans <- list(small=small, large=large, d.logL=d.logL, d.df=d.df, 
              p.value=1-pchisq(d.logL, df=d.df),
              k=k, 
              aic=2*d.logL-k*d.df)      
    
  class(ans) <- "compareiModels"
  ans
}

print.compareiModels <- function(x,...){
  cat("Large:\n")
  str(x$large, give.head=FALSE,no.list=TRUE,comp.str=" ")
  cat("Small:\n")
  str(x$small, give.head=FALSE,no.list=TRUE,comp.str=" ")
  cat(sprintf("-2logL: %8.2f df: %i AIC(k=%4.1f): %8.2f p.value: %f\n",
              2*x$d.logL, x$d.df, x$k, x$aic, x$p.value))
}
