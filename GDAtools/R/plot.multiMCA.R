plot.multiMCA <- function(x,type='v',axes=c(1,2),points='all',groups=1:x$call$ngroups,col=rainbow(x$call$ngroups),app=0, ...) {
  tit1 <- paste('Dim ',axes[1],' (',round(x$eig[[2]][axes[1]],1),'%)',sep='')
  tit2 <- paste('Dim ',axes[2],' (',round(x$eig[[2]][axes[2]],1),'%)',sep='')
  if (type=='v') {
    cmin <- apply(do.call('rbind',lapply(x$VAR,function(x) apply(x$coord[,axes],2,min))),2,min)*1.1
    cmax <- apply(do.call('rbind',lapply(x$VAR,function(x) apply(x$coord[,axes],2,max))),2,max)*1.1
    clim <- cbind(cmin,cmax)
    plot(x$ind$coord,col='white',xlim=clim[1,],ylim=clim[2,],xlab=tit1,ylab=tit2,...)
    abline(h=0,v=0,col='grey')
    for(i in 1:length(groups)) {
      var <- x$VAR[[groups[i]]]$coord
      if(points=='all') condi <- 1:nrow(var)
      if(points=='besth') condi <- abs(x$VAR[[groups[i]]]$v.test[,axes[1]])>=2.58
      if(points=='bestv') condi <- abs(x$VAR[[groups[i]]]$v.test[,axes[2]])>=2.58
      if(points=='best') condi <- abs(x$VAR[[groups[i]]]$v.test[,axes[1]])>=2.58 | abs(x$VAR[[groups[i]]]$v.test[,axes[2]])>=2.58
      coord <- var[condi,axes]
      prop <- round(x$VAR[[groups[i]]]$weight[-x$my.mca[[groups[i]]]$call$excl]/nrow(x$ind$coord)*2+0.5,1)[condi]
      if(app==0) text(coord,rownames(coord),col=col[i],cex=1)
      if(app==1) text(coord,rownames(coord),col=col[i],cex=prop)
      if(app==2) {
	points(coord,pch=17,col=col[i],cex=prop)
	text(coord,rownames(coord),pos=3,col=col[i],cex=1)
	}
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
    abline(h=0,v=0,col='grey')
    }
}
