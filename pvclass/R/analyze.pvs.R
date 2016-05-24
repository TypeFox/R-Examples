analyze.pvs <-
function(pv, Y = NULL, alpha = 0.05, roc = TRUE, pvplot = TRUE, cex = 1){
  if(!is.null(Y)){
    Y <- factor(Y)
    T <- prtable(Y = unclass(Y), pv = pv, alpha = alpha)
    if(roc == TRUE & pvplot==TRUE){ par(ask = TRUE) }
    if(roc == TRUE){ rocs(Y = unclass(Y), pv = pv, cex = cex) }
    if(pvplot == TRUE){ pvplot(Y = Y, pv = pv, alpha = alpha) }
    invisible(T)
  } else{
    pvplot(pv = pv, alpha = alpha)
  }
}
  
