plot.Predict.Treat.T0T1.ContCont <- function(x, Xlab, Main, alpha=0.05, Cex.Legend=1, ...){ 
  
  Object <- x 
  
  if (missing(Xlab)) {Xlab <- expression(paste(Delta, "T"[j], "|S"[j]))}
  if (missing(Main)) {Main=" "} 
  
  min_T0T1 <- min(Object$T0T1)
  max_T0T1 <- max(Object$T0T1)
  
  crit_val <- qnorm(c(alpha/2), mean=0, sd=1, lower.tail=FALSE)
  
  if (class(Object)=="Predict.Treat.T0T1.ContCont") {
    user_req_SD_Delta.T_givenS <- sqrt(Object$Var_Delta.T_S)
    Pred_T <- Object$Pred_T 
    x <- ((seq(-4,4,length=1000)*user_req_SD_Delta.T_givenS))*1 + Pred_T
    hx <- dnorm(x, Pred_T, user_req_SD_Delta.T_givenS)
    plot(x=x, y=(hx*crit_val), type="l", xlab=Xlab, ylab="", main=Main, lwd=2, col=2, ...)
    abline(v = Pred_T, lty=2)
    
    SD_hier <- user_req_SD_Delta.T_givenS
    col_hier <- 2
    x1 <- Pred_T - (SD_hier * crit_val)
    x2 <- Pred_T + (SD_hier * crit_val)
    sam <- cbind(x, hx)
    y1 <- sam[,2][which.min(abs(sam[,1]-x1))]
    y2 <- sam[,2][which.min(abs(sam[,1]-x2))]
    
    segments(x0 = x1, y0 = 0, x1 = x1, y1 = y1, col=col_hier, lwd=2, lty=2)
    segments(x0 = x2, y0 = 0, x1 = x2, y1 = y2, col=col_hier, lwd=2, lty=2)
    
    
    
    legend("topright", inset=.05, legend=c(bquote(paste(rho[T0T1], "="~.(Object$T0T1))), expression()), 
           lwd=2, col=c(2), cex=Cex.Legend)
  }  
  
}
