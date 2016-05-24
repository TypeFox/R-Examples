`accumplot` <-
function(xr,addit=F,labels="",col=1,ci=2,pch=1,type="p",cex=1,xlim=c(1,xmax),ylim=c(1,rich),xlab="sites",ylab="species richness",...) {
    x <- xr
    xmax <- max(x$sites)
    rich <- max(x$richness)    
    if(addit==F) {graphics::plot(x$sites,x$richness,xlab=xlab,ylab=ylab,bty="l",type=type,col=col,pch=pch,cex=cex,xlim=xlim,ylim=ylim)} 
    if(addit==T) {graphics::points(x$sites,x$richness,type=type,col=col,pch=pch,cex=cex)}
    graphics::plot(x,add=T,ci=ci,col=col,...)
    if(labels!="") {
        l <- length(x$sites)
        graphics::text(x$sites[1],x$richness[1],labels=labels,col=col,pos=2,cex=cex)
        graphics::text(x$sites[l],x$richness[l],labels=labels,col=col,pos=4,cex=cex)
    }
}

