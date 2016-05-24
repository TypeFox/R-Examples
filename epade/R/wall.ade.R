wall.ade <-
function(vnames=NULL, main=NULL, xlab=NULL, ylab=NULL, glab=NULL, legendon='topright', xlim=NULL, ylim=NULL, lwd=1, pch=16,  lty=1, xticks=NULL, yticks=NULL, col=NULL, tcol=NULL,  bgcol=NULL, lcol=NULL, wall=0, v=NULL, h=NULL){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr', 'plt',  'pin', 'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))

if(is.null(xlim))  xlim<-c(-1, 1)
if(is.null(ylim))  ylim<-c(-1, 1)

vnames2<-NULL
if(!is.null(vnames) & is.list(vnames)){
vnames1 <-vnames[[1]]
vnames2 <-vnames[[2]]
vnames  <-vnames1
}


#####################################
# Colors
if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'

if(is.null(col) & wall!=0 & is.null(vnames)) col <- rgb(0.3,0.3,0.45)
if(is.null(col) & wall==0 & is.null(vnames)) col <- 'gray30'

if(!is.null(vnames)){
a.Ng<-length(vnames)
if(is.null(col)) col<-a.getcol.ade(a.Ng)
}

if(is.null(lcol))  lcol<- tcol
#####################################





#####  Style 0 #################################################################
if(wall==0){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=3){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes=FALSE)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1, at=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks), labels=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks))
if(!is.null(yticks) & length(yticks)==1) a2<-axis(2, col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1, at=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks), labels=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks))
if(!is.null(xticks) & length(xticks)>1) a1<-axis(1, col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1, at=xticks, labels=xticks)
if(!is.null(yticks) & length(yticks)>1) a2<-axis(2, col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1, at=yticks, labels=yticks)
if(is.null(xticks))  a1<-axis(1, col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1)
if(is.null(yticks))  a2<-axis(2, col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)
box( col=bgcol)
title(main)
}

#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
legend(legendon, legend=g,  pt.cex=2, title=glab, pch=c(rep(pch, length(g)))  ,lwd=c(rep(0, length(g)))  ,col = col, box.lwd=1, lty = c(rep(0, length(g))),  merge = TRUE, box.col=bgcol, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(c(g, glab),font = 2)))
}


}
################################################################################



#####  Style 1 #################################################################
if(wall==1){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=3){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim ,  axes=FALSE)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks), labels=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks))
if(!is.null(yticks) & length(yticks)==1) a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks), labels=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks))
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=xticks, labels=xticks)
if(!is.null(yticks) & length(yticks)>1)  a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=yticks, labels=yticks)
if(is.null(xticks))  a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
if(is.null(yticks))  a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)
box(lwd=1, col=rgb(1,1,1))
title(main)
}

#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
n<- length(g)
legend(legendon, legend=g,  pt.cex=2,   title=glab,pch=c(rep(pch, n)) , col = a.alpha.ade(col, 1), lty = c(rep(0, n)), box.col=rgb(1,1,1) , box.lwd=2, text.col=tcol, bg=bgcol,  text.width=max(strwidth(c(g, glab),font = 2)))
}
}
################################################################################


#####  Style 2 #################################################################
if(wall==2){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=3){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim ,  axes=FALSE)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1, at=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks), labels=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks))
if(!is.null(yticks) & length(yticks)==1) a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1, at=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks), labels=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks))
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1, at=xticks, labels=xticks)
if(!is.null(yticks) & length(yticks)>1)  a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1, at=yticks, labels=yticks)
if(is.null(xticks)) a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
if(is.null(yticks)) a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=bgcol, lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)
box(col=a.coladd.ade(bgcol, -75))
title(main)
}


#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
n<- length(g)
legend(legendon, legend=g,  pt.cex=2,  title=glab, pch=c(rep(pch, n)) ,col = a.alpha.ade(a.coladd.ade(col, 10), 1), lty = c(rep(0, n)), box.col=a.coladd.ade(bgcol, -75) , box.lwd=1, text.col=tcol, bg=rgb(1,1,1),  text.width=max(strwidth(c(g, glab),font = 2)))
legend(legendon, legend=g,  pt.cex=2,  title=glab, bg=rgb(1,1,1,0),    pch=c(rep(pch-15, n)) , col = a.alpha.ade(a.coladd.ade(col, -50), 1), lty = c(rep(0, n)), box.col=a.coladd.ade(bgcol, -75) , box.lwd=0, text.col=tcol,  text.width=max(strwidth(c(g, glab),font = 2)))
}
}
################################################################################


#####  Style 3 #################################################################
if(wall==3){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)



#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=3){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=rgb(1,1,1), axes=FALSE)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -50), lwd.ticks=1, at=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks), labels=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks))
if(!is.null(yticks) & length(yticks)==1) a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -50), lwd.ticks=1, at=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks), labels=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks))
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -50), lwd.ticks=1, at=xticks, labels=xticks)
if(!is.null(yticks) & length(yticks)>1)  a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -50), lwd.ticks=1, at=yticks, labels=yticks)
if(is.null(xticks))  a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -50), lwd.ticks=1)
if(is.null(yticks))  a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -50), lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=a.coladd.ade(bgcol, -50), lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)
box(lwd=1, col=a.coladd.ade(bgcol, -50))
title(main)
}


#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
n<- length(g)
legend(legendon, legend=g,  pt.cex=2, title=glab,  pch=c(rep(pch, n)) ,col = a.alpha.ade(a.coladd.ade(col, 10), 1), lty = c(rep(0, n)), box.col=a.coladd.ade(bgcol, -50) , box.lwd=1, text.col=tcol, bg=rgb(1,1,1),  text.width=max(strwidth(c(g, glab),font = 2)))
legend(legendon, legend=g,  pt.cex=2, title=glab,  bg=rgb(1,1,1,0),  pch=c(rep(pch-15, n)) , col = a.alpha.ade(a.coladd.ade(col, -50), 1), lty = c(rep(0, n)), box.col=a.coladd.ade(bgcol, -50) , box.lwd=0, text.col=tcol,  text.width=max(strwidth(c(g, glab),font = 2)))
}
}
################################################################################




#####  Style 4 #################################################################
if(wall==4){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
par(col.main=rgb(1,1,1))
par(font=2)


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=3){
plot(0, 0, type='s', xlab='', ylab='', xlim=xlim, ylim=ylim ,  axes=FALSE)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks), labels=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks))
if(!is.null(yticks) & length(yticks)==1) a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks), labels=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks))
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=xticks, labels=xticks)
if(!is.null(yticks) & length(yticks)>1)  a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=yticks, labels=yticks)
if(is.null(xticks)) a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
if(is.null(yticks)) a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)

# Outer
par(xpd=TRUE)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0, 2.75,  2.75, 0)), col=tcol, border=rgb(1,1,1))
if(!is.null(ylab)) if(ylab!='' & ylab!=' ') polygon( a.glc(side=2, line=c(3.5, 3.5, 2, 2)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
if(!is.null(ylab)) if(xlab!='' & xlab!=' ') polygon( a.glc(side=c(2, 2, 4, 4), line=0),     a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))
text(a.glc(side=0), a.glc(side=3, line=1),    labels=main, cex = 1.25, font=2, col=rgb(1,1,1), adj=c(0.5,0))
text(a.glc(side=0), a.glc(side=1, line=3.5),  labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=2.5), a.glc(side=5),  labels=ylab, cex = 1.1,  font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)
box(lwd=1, col=rgb(1,1,1))
}


#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
n<- length(g)
legend(legendon, legend=g,  pt.cex=2, title=glab, pch=c(rep(pch, n)) ,      col = a.alpha.ade(col, 1), lty = c(rep(0, n)), box.col=rgb(1,1,1) , box.lwd=1, text.col=rgb(1,1,1), bg=tcol,  text.width=max(strwidth(c(g, glab),font = 2)))
legend(legendon, legend=g,  pt.cex=2, title=glab, pch=c(rep(pch-15, n)) ,   col = rgb(1,1,1), lty = c(rep(0, n)), box.col=rgb(1,1,1, 0) , box.lwd=1, text.col=rgb(1,1,1, 0), bg=rgb(1,1,1,0),  text.width=max(strwidth(c(g, glab),font = 2)))

}
}
################################################################################



#####  Style 5 #################################################################
if(wall==5){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
par(col.main=tcol)
par(font=2)
newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[2]<1) newmai[2]<- 1 - oldmai[2]
if(oldmai[3]>0.75 & oldmai[3]<=0.82) newmai[3]<- 0.75-oldmai[3]
if(oldmai[4]>0.25 & oldmai[4]<=0.42) newmai[4]<- 0.25-oldmai[4]
par(mai=(oldmai+newmai))


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=3){
plot(0, 0, type='s', xlab='', ylab='', xlim=xlim, ylim=ylim ,  axes=FALSE)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks), labels=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks))
if(!is.null(yticks) & length(yticks)==1) a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks), labels=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks))
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=xticks, labels=xticks)
if(!is.null(yticks) & length(yticks)>1)  a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=yticks, labels=yticks)
if(is.null(xticks)) a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
if(is.null(yticks)) a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)
# Outer
par(xpd=TRUE)
polygon(a.glc(side=2, line=c(4.25, 4.25, 0, 0)), a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)),   a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
polygon(a.glc(side=2, line=c(4.25, 4.25 ,3.65, 3.65)),  a.glc(side=c(1,3,3,1), line=c(2.6, 0.6, 0.6, 2.6)), col=bgcol,  border=tcol)
polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=tcol)
polygon(a.glc(side=2, line=c(4.25, 4.25, 0, 0)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
polygon(a.glc(side=c(2, 2, 4, 4), line=0), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0, 0.6, 0.6)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
text(a.glc(side=0), a.glc(side=3, line=1.5),  labels=main, cex = 1.25, font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=0), a.glc(side=1, line=3.75), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=2.5), a.glc(side=5), labels=ylab, cex = 1.1,   font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)
box(lwd=1, col=tcol)

}


#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
n<- length(g)
xr<-diff(par('usr')[1:2])/10
yr<-diff(par('usr')[3:4])/10
legend(legendon, legend=g,  pt.cex=2, title=glab, pch=c(rep(pch, n)) ,  col = a.alpha.ade(col, 1), lty = c(rep(0, n)), box.col=tcol , box.lwd=1, text.col=tcol, bg=rgb(1,1,1),  text.width=max(strwidth(c(g, glab),font = 2)))
legend(legendon, legend=g,  pt.cex=2, title=glab,pch=c(rep(pch-15, n)) ,col = tcol, lty = c(rep(0, n)), box.col=tcol , box.lwd=1, text.col=tcol, bg=rgb(1,1,1,0),  text.width=max(strwidth(c(g, glab),font = 2)))
}
}
################################################################################



#####  Style 6 #################################################################
if(wall==6){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=3){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim ,  axes=FALSE)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3, at=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks), labels=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks))
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, col=rgb(1,1,1), col.ticks=rgb(1,1,1), lwd.ticks=1, at=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks), labels=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks))
if(!is.null(yticks) & length(yticks)==1) a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3, at=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks), labels=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks))
if(!is.null(yticks) & length(yticks)==1) a2<-axis(2, col=rgb(1,1,1), col.ticks=rgb(1,1,1), lwd.ticks=1, at=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks), labels=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks))
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3, at=xticks, labels=xticks)
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, col=rgb(1,1,1), col.ticks=rgb(1,1,1), lwd.ticks=1, at=xticks, labels=xticks)
if(!is.null(yticks) & length(yticks)>1)  a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3, at=yticks, labels=yticks)
if(!is.null(yticks) & length(yticks)>1)  a2<-axis(2, col=rgb(1,1,1), col.ticks=rgb(1,1,1), lwd.ticks=1, at=yticks, labels=yticks)
if(is.null(xticks))  a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
if(is.null(xticks))  a1<-axis(1, col=rgb(1,1,1), col.ticks=rgb(1,1,1), lwd.ticks=1)
if(is.null(yticks))  a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
if(is.null(yticks))  a2<-axis(2, col=rgb(1,1,1), col.ticks=rgb(1,1,1), lwd.ticks=1)

abline(v=a1, h=a2, lty=1, col=a.coladd.ade(bgcol, -35), lwd=3)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)
box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=a.coladd.ade(bgcol, -35))
title(main)
}

#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
n<- length(g)
legend(legendon, legend=g,  pt.cex=2,   title=glab,pch=c(rep(pch, n)) , col = a.alpha.ade(col, 1), lty = c(rep(0, n)), box.col=a.coladd.ade(bgcol, -35) ,  box.lwd=3, text.col=tcol, bg=bgcol,  text.width=max(strwidth(c(g, glab),font = 2)))
legend(legendon, legend=g,  pt.cex=2,   title=glab,pch=c(rep(pch, n)) , col = a.alpha.ade(col, 1), lty = c(rep(0, n)), box.col=rgb(1,1,1) ,                box.lwd=1, text.col=tcol, bg=bgcol,  text.width=max(strwidth(c(g, glab),font = 2)))

box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=a.coladd.ade(bgcol, -35))
}
}
################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


################################################################################
# Plotaufruf
plot.box.ade( xlab=xlab, ylab=ylab, main=main, xlim=xlim, ylim=ylim, lwd=lwd)

if(!is.null(vnames)){
legens.ade(ylim[1], vnames, xlim=xlim, ylim=ylim, lwd=lwd, pch=pch, col=col, legendon=legendon)
}
################################################################################         



}
