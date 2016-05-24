plot.csMCA <-
 function(x,type='v',axes=1:2,points='all',col='dodgerblue4', app=0, ...) {
  tit1 <- paste('Dim ',axes[1],' (',round(x$eig$mrate[axes[1]],1),'%)',sep='')
  tit2 <- paste('Dim ',axes[2],' (',round(x$eig$mrate[axes[2]],1),'%)',sep='')
  if (type=='v') {
    cmin <- apply(x$var$coord[,axes],2,min)*1.1
    cmax <- apply(x$var$coord[,axes],2,max)*1.1
    clim <- cbind(cmin,cmax)
    nv <- nrow(x$var$coord)
    if(points=='all') condi <- 1:nv
    if (points=='besth') condi <- x$var$contrib[,axes[1]]>=100/nv
    if (points=='bestv') condi <- x$var$contrib[,axes[2]]>=100/nv
    if (points=='best') condi <- x$var$contrib[,axes[1]]>=100/nv | x$var$contrib[,axes[2]]>=100/nv
    coord <- x$var$coord[condi,axes]
    prop <- round(x$var$weight[-x$call$excl]/nrow(x$ind$coord)*2+0.5,1)[condi]
    plot(coord,col='white',xlim=clim[1,],ylim=clim[2,],xlab=tit1,ylab=tit2,...)
    if(app==0) text(coord,rownames(coord),col=col,cex=1)
    if(app==1) text(coord,rownames(coord),col=col,cex=prop)
    if(app==2) {
	points(coord,pch=17,col=col,cex=prop)
	text(coord,rownames(coord),pos=3,col=col,cex=1)
	}
    }
  if (type %in% c('i','inames')) { # new
    cmin <- apply(x$ind$coord[,axes],2,min)*1.1
    cmax <- apply(x$ind$coord[,axes],2,max)*1.1
    clim <- cbind(cmin,cmax)
    ni <- nrow(x$ind$coord) # new
    if(points=='all') condi <- 1:ni # new
    if (points=='besth') condi <- x$ind$contrib[,axes[1]]>=100/ni # new
    if (points=='bestv') condi <- x$ind$contrib[,axes[2]]>=100/ni # new
    if (points=='best') condi <- x$ind$contrib[,axes[1]]>=100/ni | x$ind$contrib[,axes[2]]>=100/ni # new
    coord <- x$ind$coord[condi,axes] # new
    if(type=='i') pcol <- col #new
    if(type=='inames') pcol <- 'white' #new
    plot(coord,col=pcol,xlim=clim[1,],ylim=clim[2,],xlab=tit1,ylab=tit2,pch=19,cex=0.2,...) #new
    if(type=='inames') text(coord,rownames(coord),col=col) #new
    }
  abline(h=0,v=0,col='grey')
}
