plotTR <-function(object,...){
  #TODO : *Add CI around point estimates
#  require(lattice) # MJM20141101: lattice is imported
  #plot trend over time for all items
  itms <- object$itms
  tps <- object$mpoints
  pplgrps <- object$ngroups/itms
  trend <- object$etapar[((pplgrps-1)*itms*(tps-1)+1):((pplgrps-1)*itms*(tps-1)+(itms*(tps-1)))]
  tips <-rep(paste("t",1:tps,sep=""),each=itms)
  items <- rep(paste("Item",1:itms),tps)
  tr0 <- rep(0,itms)
  trend <- c(tr0,trend)
  plotdats <- data.frame(trend,items,tips)
  key.items <- list(space = "right", text = list(levels(plotdats$items)),
                   points = list(pch = 1:length(levels(plotdats$items)),
                   col = "black")
                   )
 plotout <- xyplot(trend~tips,data=plotdats,
                   aspect="fill", type="o",
                   groups=items, 
                   key=key.items,
                   lty=1,pch = 1:length(levels(plotdats$items)),
                   col.line = "darkgrey", col.symbol = "black",
                   xlab = "Time",
                   ylab = "Effect", 
                   main = "Trend effect plot for LLRA"
                   )
  print(plotout)
}
