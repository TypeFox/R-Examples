plotGR <- function(object,...){
  #TODO: *Add CI around point estimates
#  require(lattice) # MJM20141101: lattice is imported
  itms <- object$itms
  tps <- object$mpoints
  pplgrps <- object$ngroups/itms
  if(pplgrps<2) stop("There are no treatment effects in this analysis.")
  
  #treatment effects for all treatment groups at tps>1
  treat <- object$etapar[1:((pplgrps-1)*itms*(tps-1))]
  time <- factor(rep(paste("t",2:tps,sep=""),each=itms*(pplgrps-1)))
  item <- factor(rep(rep(paste("Item",1:itms),each=pplgrps-1),tps-1))
  names1 <- unique(names(object$groupvec))[1:(length(unique(names(object$groupvec))))-1]
  #labeling
  group <- factor(rep(names1,itms*(tps-1)))
  plotdats1 <- data.frame(treat,group,item,time)
  
  #effects (i.e. zeros) for all treatment groups at tp=1
  treat0 <- rep(0,itms*(pplgrps-1))
  time0 <- factor(rep("t1",each=itms*(pplgrps-1)))
  item0 <- factor(rep(paste("Item",1:itms),each=(pplgrps-1)))
  #labeling
  group0 <- factor(rep(names1,itms))
  plotdats0 <- data.frame(treat0,group0,item0,time0)
  names(plotdats0) <- c("treat","group","item","time")

  #effects (i.e. zeros) for control or baseline group for all tps 
  treat00 <- rep(0,itms*tps)
  time00 <- factor(rep(paste("t",1:tps,sep=""),each=itms))
  item00 <- factor(rep(paste("Item",1:itms),tps))
  group00 <- factor(rep(unique(names(object$groupvec))[length(unique(names(object$groupvec)))],itms*tps))
  plotdats00 <- data.frame(treat00,group00,item00,time00)
  names(plotdats00) <- c("treat","group","item","time")

  #all together
  plotdats <- rbind(plotdats00,plotdats0,plotdats1)

  #plot
  key.group <- list(space = "right", text = list(levels(plotdats$group)),
                    points = list(pch = 1:length(levels(plotdats$group)),
                    col = "black")
                    )
  plotout <- xyplot(treat ~ time | item, plotdats,
                    aspect = "xy", type = "o", 
                    groups = group, key = key.group,
                    lty = 1, pch = 1:length(levels(plotdats$group)),
                    col.line = "darkgrey", col.symbol = "black",
                    xlab = "Time",
                    ylab = "Effect", 
                    main = "Treatment effect plot for LLRA"
                    )
  print(plotout)
}
