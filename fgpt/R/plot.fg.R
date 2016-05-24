plot.fg <-
function(x, plane=0, plot.zero=TRUE, sign.level=0.05, ...){
  
  if (class(x)!="fg") {stop("Object needs to be of class \"fg\"")}
  if (plane!=0 & plane!=1 & plane!=2) {stop("plane needs to be either 0,1 or 2")}
  if (!is.logical(plot.zero)) {stop("plot.zero needs to be logical")}
  
  p.val <- function(perm,obsv){
    return(min((sum(perm>=obsv))/(sum(is.na(perm)==FALSE))*2,(sum(perm<=obsv))/(sum(is.na(perm)==FALSE))*2,1))   
  }
  
  par.fig <- par("fig")
  par.mar <- par("mar")
  par.ask <- par("ask")
 
  if(plot.zero==TRUE){
    xlims <- range(c(0,x$scales[-length(x$scales)]))
  }else{
    xlims <- range(x$scales)
  }
  
  if(x$pairwise==FALSE){
    
    pvals <- sapply(x$c.list,p.val,obsv=x$c.list[[1]][1])
    cols <- ifelse(pvals<=sign.level,"blue","black")

    par(mfrow=c(1,1))
    par(fig=c(0,0.9,0,1))
    par(mar=c(5,5,2,0))
    ylims <- range(x$c.list[[1]][1]-unlist(x$c.list))
    plot(x$scales,x$c.list[[1]][1]-sapply(x$c.list,mean,na.rm = TRUE), ylim=ylims, xlim=xlims, pch=16, col=cols[-length(x$scales)], xlab="grid cell size", ylab="unexplained Moran's I")
    segments(x$scales,x$c.list[[1]][1]-sapply(x$c.list,quantile,probs=0.025,na.rm = TRUE),x$scales,x$c.list[[1]][1]-sapply(x$c.list,quantile,probs=0.975,na.rm = TRUE), col=cols[-length(x$scales)])
    abline(h=0, lwd=2, col="red")
    
    par(fig=c(0.9,1,0,1),new=TRUE)
    par(mar=c(5,0.2,2,2))
    plot(0,x$c.list[[1]][1]-sapply(x$c.list,mean,na.rm = TRUE)[length(x$scales)], ylim=ylims, pch=16, col=cols[length(x$scales)], xlab="", ylab="",xaxt="n", yaxt="n")
    axis(1, at=0, labels=quote(infinity), cex.axis=1.5)
    segments(0,x$c.list[[1]][1]-sapply(x$c.list,quantile,probs=0.025,na.rm = TRUE)[length(x$scales)],0,x$c.list[[1]][1]-sapply(x$c.list,quantile,probs=0.975,na.rm = TRUE)[length(x$scales)],col=cols[length(x$scales)])
    abline(h=0, lwd=2, col="red")
    
  }else if (x$correlate==FALSE){
    if(plane %in% c(0,1)){
      
      pvals <- sapply(x$m.list,p.val,obsv=x$m.list[[1]][1])
      cols <- ifelse(pvals<=sign.level,"blue","black")
      
      par(mfrow=c(1,1))
      par(fig=c(0,0.9,0,1))
      par(mar=c(5,5,2,0))
      ylims <- range(unlist(x$m.list))
      plot(x$scales,sapply(x$m.list,mean,na.rm = TRUE), ylim=ylims, xlim=xlims, pch=16, col=cols[-length(x$scales)], xlab="grid cell size", ylab="mean")
      segments(x$scales,sapply(x$m.list,quantile,probs=0.025,na.rm = TRUE),x$scales,sapply(x$m.list,quantile,probs=0.975,na.rm = TRUE), col=cols[-length(x$scales)])
      abline(h=x$m.list[[1]][1], lwd=2, col="red")
    
      par(fig=c(0.9,1,0,1),new=TRUE)
      par(mar=c(5,0.2,2,2))
      plot(0,sapply(x$m.list,mean,na.rm = TRUE)[length(x$scales)], ylim=ylims, pch=16, xlab="", ylab="",xaxt="n", yaxt="n",col=cols[length(x$scales)])
      axis(1, at=0, labels=quote(infinity), cex.axis=1.5)
      segments(0,sapply(x$m.list,quantile,probs=0.025,na.rm = TRUE)[length(x$scales)],0,sapply(x$m.list,quantile,probs=0.975,na.rm = TRUE)[length(x$scales)],col=cols[length(x$scales)])
      abline(h=x$m.list[[1]][1], lwd=2, col="red")
      }
  
    if(plane==0){par(ask=TRUE)}
  
    if(plane %in% c(0,2)){
      
      pvals <- sapply(x$v.list,p.val,obsv=x$v.list[[1]][1])
      cols <- ifelse(pvals<=sign.level,"blue","black")
      
      par(mfrow=c(1,1))
      par(fig=c(0,0.9,0,1))
      par(mar=c(5,5,2,0))
      ylims <- range(unlist(x$v.list))
      plot(x$scales,sapply(x$v.list,mean,na.rm = TRUE), ylim=ylims, xlim=xlims, pch=16, col=cols[-length(x$scales)], xlab="grid cell size", ylab="variance")
      segments(x$scales,sapply(x$v.list,quantile,probs=0.025,na.rm = TRUE),x$scales,sapply(x$v.list,quantile,probs=0.975,na.rm = TRUE), col=cols[-length(x$scales)])
      abline(h=x$v.list[[1]][1], lwd=2, col="red")
    
      par(fig=c(0.9,1,0,1),new=TRUE)
      par(mar=c(5,0.2,2,2))
      plot(0,sapply(x$v.list,mean,na.rm = TRUE)[length(x$scales)], ylim=ylims, pch=16, xlab="", ylab="",xaxt="n", yaxt="n", col=cols[length(x$scales)])
      axis(1, at=0, labels=quote(infinity), cex.axis=1.5)
      segments(0,sapply(x$v.list,quantile,probs=0.025,na.rm = TRUE)[length(x$scales)],0,sapply(x$v.list,quantile,probs=0.975,na.rm = TRUE)[length(x$scales)], col=cols[length(x$scales)])
      abline(h=x$v.list[[1]][1], lwd=2, col="red")
    }

  }else{
    
    pvals <- sapply(x$c.list,p.val,obsv=x$c.list[[1]][1])
    cols <- ifelse(pvals<=sign.level,"blue","black")
    
    par(mfrow=c(1,1))
    par(fig=c(0,0.9,0,1))
    par(mar=c(5,5,2,0))
    
    ylims <- range(x$c.list[[1]][1]-unlist(x$c.list))
    plot(x$scales,sapply(x$c.list,mean,na.rm = TRUE), ylim=ylims, xlim=xlims, pch=16, xlab="grid cell size", ylab="correlation", col=cols[-length(x$scales)])
    segments(x$scales,sapply(x$c.list,quantile,probs=0.025,na.rm = TRUE),x$scales,sapply(x$c.list,quantile,probs=0.975,na.rm = TRUE), col=cols[-length(x$scales)])
    abline(h=x$c.list[[1]][1], lwd=2, col="red")
    
    par(fig=c(0.9,1,0,1),new=TRUE)
    par(mar=c(5,0.2,2,2))
    plot(0,sapply(x$c.list,mean,na.rm = TRUE)[length(x$scales)], ylim=ylims, pch=16, xlab="", ylab="",xaxt="n", yaxt="n", col=cols[length(x$scales)])
    axis(1, at=0, labels=quote(infinity), cex.axis=1.5)
    segments(0,sapply(x$c.list,quantile,probs=0.025,na.rm = TRUE)[length(x$scales)],0,sapply(x$c.list,quantile,probs=0.975,na.rm = TRUE)[length(x$scales)], col=cols[length(x$scales)])
    abline(h=x$c.list[[1]][1], lwd=2, col="red")
    }
  par(fig=par.fig)
  par(mar=par.mar)
  par(ask=par.ask)
  }
