plotLTT <- function(trees,col=rgb(.5,.5,.5,.5),xlab="Time",
                    ylab="Lineages",log="y",add=FALSE,
                    plot=TRUE,plot.tree=NA,type="s",...) 
{
  max.time <- sapply(trees,function(t) max(t[,1]))
  extant <- sapply(trees,function(t) sum(2*t[,2]-1))
  lineages <- lapply(trees,function(t) sum(2*t[,2]-1)+cumsum(1-2*t[,2]))
  max.lineages <- max(sapply(lineages,max))
  if (! plot) {
    return(lineages)
  } else {
    if (! add) {
      if (! is.na(plot.tree) & length(trees) == 1) {
        num.tips <- length(plot.tree[[1]]$tip.label)
        scale.factor <- 1.0*num.tips/log10(max.lineages)
        #ll <- lineages[[1]]*scale.factor
        ll <- log10(lineages[[1]])*scale.factor
        par(mar=c(4,4,1,1),oma=c(1,1,0,0))
        plot(plot.tree[[1]],show.tip.label=FALSE,
             edge.color=rgb(.3,.3,.3),root.edge=TRUE)
        lines(max.time-trees[[1]][,1],ll,type=type,lwd=5,col="white")
        lines(max.time-trees[[1]][,1],ll,type=type,lwd=3,col="black")
        xax <- pretty(-trees[[1]][,1])
        yax <- c(pretty(lineages[[1]]))
        if (0 %in% yax) { yax[yax==0] <- 1 }
        #axis(1,at=xax+max.time,labels=xax,pos=log(scale.factor)/log(num.tips)*num.tips)
        usr <- par("usr")
        dx <- (usr[2]-usr[1])
        dy <- (usr[4]-usr[3])
        inset <- 0.005
        lines(c(usr[1]+dx*inset,usr[2]),rep(usr[3]+dy*inset,2))
        lines(rep(usr[1]+dx*0.005,2),c(usr[3]+dy*inset,usr[4]))
        axis(1,at=xax+max.time,labels=xax,lwd=0)
        axis(2,at=log10(yax)*scale.factor,labels=yax,las=1,lwd=0)
        mtext("Time",1,outer=TRUE,cex=1.4,line=-1)
        mtext("Lineages",2,outer=TRUE,cex=1.4,line=-1)
      } else {
        plot(1,type="n",xlim=c(-max(max.time),0),ylim=c(1,max(unlist(lineages))),log=log,
             xlab=xlab,ylab=ylab,...)
        for (i in 1:length(lineages)) {
          if ("y" %in% strsplit(log,"")[[1]]) {
            lineages[[i]][lineages[[i]]==0] <- 0.5
          }
          lines(-rev(trees[[i]][,1]),rev(lineages[[i]]),col=col,type=type,...)
        }
        #lines(c(-trees[[i]][,1])rep(lineages[[i]][length(lineages[[i]])],2))
      }
    } else {
      for (i in 1:length(lineages)) {
        lines(-trees[[i]][,1],lineages[[i]],col=col,type=type,...)
      }
    }
  }
}

