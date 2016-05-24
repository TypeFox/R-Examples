KM.plot.ade <-
function( time, event, group=NULL, data=NULL, vnames=NULL, main='Kaplan-Meier Plot', xlab='Follow-Up Time', ylab='Cumulative Survival',xlim=NULL, ylim=NULL, xticks=NULL, legendon='bottomleft', lwd=2, lty=1, col=NULL, tcol=NULL, bgcol=NULL, pdigs=4, CI=FALSE, ycut=TRUE, zenspoints=FALSE, test=FALSE, wall=0){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr',  'plt',  'pin',  'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))


library(survival)
library(plotrix)


################################################################################
# Ohne Dataframe
ismitdata=TRUE
if(is.numeric(time) & is.numeric(event)){
ismitdata=FALSE
if(length(time)!=length(event))  stop("time and event length different.")
data<-time
data<-as.data.frame(data)
data$zmy<-event
data$tmy<-time
if(!is.null(group)){
eval(parse(text=paste( 'data$', gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(group))), '<-group')))
group<-gsub('.*[$]', '' ,deparse(substitute(group)))
}
time<-'tmy'
event<-'zmy'
}
################################################################################



if (is.null(vnames)) vnames <- group
if(ismitdata){if(!is.data.frame(data))  stop("(data) must be a data.frame!") }
if(!is.null(group)) g<-eval(parse(text=paste("data$",group, sep='')))
if(!is.null(group)) g <- as.factor(g)
if(is.null(group))  xrange<-range(eval(parse(text=paste("data$",time, sep=''))), na.rm=TRUE)
if(!is.null(group)) xrange<-range(eval(parse(text=paste("data$",time, sep='')))[!is.na(g)], na.rm=TRUE)
if(is.null(xlim))   xlim<-c(0, xrange[2])



################################################################################
# Colors
a.N<-1
if(!is.null(group)) a.N<-nlevels(g)

if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'

if(is.null(col) & a.N==1) col <- tcol
if(is.null(col) & a.N> 1) col <- a.getcol.ade(a.N)


################################################################################



#####  Function  Zenst  ########################################################  
zenst<-function(s,t,t2){   
k=1;
s2<-NULL
t<-sort(t)
for(i in 1:length(t)){
if(k<length(t2)){ while(t[i]>=t2[k+1] & k<length(t2)) k<-k+1;  }
if(t[i]>=t2[k]) s2[i]<-s[k]
if(t[i]< t2[k] & k==1) s2[i]<-1
}
return(s2)
}
################################################################################


#########################
a.line.ade<-function(x,y, x2,y2, u,l, xlim, ylim, lwd=2, col=1, zenspoints=FALSE, CI=FALSE, lty=1){


# CI intervals
if(CI) {
u[is.na(u)] <-u[!is.na(u)][length(u[!is.na(u)])]
l[is.na(l)] <-l[!is.na(l)][length(l[!is.na(l)])]
alpha=0.25
if(wall==0) alpha=1
points(x[-length(x)], u,  type='s', lty=2, lwd=lwd, col=a.alpha.ade(col, alpha))
points(x[-length(x)], l,  type='s', lty=2, lwd=lwd, col=a.alpha.ade(col, alpha))
rect(x[-length(x)], l[-length(l)], x[-1], u[-length(u)], density = NULL, angle = 45, col = a.alpha.ade(col, 0.1), border = a.alpha.ade(col, 0))
}

xwd<-(lwd/3)
if(xwd<1) xwd<-1
if(zenspoints) points(x2 , y2 , cex=xwd   ,pch = "|", col=col)
points(x, y, type='s', lwd=lwd, col=col, lty=lty)
}
#########################


################################################################################
################################################################################
################################################################################

#####  Style 0 #################################################################
if(wall==0){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

#########################
#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=1){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=rgb(1,1,1,0), axes=F)
if(is.null(xticks)) axis(1, col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1)
if(!is.null(xticks) & length(xticks)==1) axis(1, at=pretty(xlim, n=xticks),  col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1)
if(!is.null(xticks) & length(xticks)>1)  axis(1, at=xticks,  col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1)
axis(2, col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1)
title(main, col=tcol)
box(col=bgcol)
if(ylim[1]>0.05) axis.break(axis=2,breakpos=ylim[1],bgcol=rgb(1,1,1),breakcol=bgcol, style="slash",brw=0.02)
}
#########################

#########################
#  Legend  #
legens.ade<-function(x, y, vnames, g, p, xlim, ylim, lwd, lty=1,  col, plot=TRUE, test){
n<- nlevels(g)
if(is.null(lty)) lty <- c(rep(1,n) ,0)
if(test)  letext<-c(paste(vnames,': ', levels(g)), paste('p Value: ', p))
if(!test) letext<-c(paste(vnames,': ', levels(g)))
lrect<-legend(legendon, inset = 0.0375,  legend=letext, plot=plot, lwd=c(rep(lwd+1, n),0), bg=rgb(1,1,1,0), col = c(col, rgb(1,1,1,0)), box.lwd=1, text.col=tcol, lty = c(lty, 0),  merge = TRUE, yjust=0, box.col=bgcol, text.width=max(strwidth(letext,font = 2)))
return(lrect)
}
#########################
}
################################################################################



#####  Style 1 #################################################################
if(wall==1){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

#########################
#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=1){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=rgb(1,1,1,0), axes=F)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(is.null(xticks)) a1<-axis(1, col=rgb(1,1,1, 0), col.ticks=tcol, lwd.ticks=1)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(xlim, n=xticks),  col=rgb(1,1,1, 0), col.ticks=tcol, lwd.ticks=1)
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, at=xticks,  col=rgb(1,1,1, 0), col.ticks=tcol, lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1, 0), col.ticks=tcol, lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)

title(main, col=tcol)
box(col=rgb(1,1,1))
if(ylim[1]>0.05 & ycut){
ktx<-diff(par('usr')[1:2])/50
kty<-diff(par('usr')[3:4])/150
polygon( c(par('usr')[1]+ktx, par('usr')[1], par('usr')[1], par('usr')[1]+ktx), c(ylim[1]+kty*4, ylim[1]+kty, ylim[1]+kty*3, ylim[1]+kty*6), col=rgb(1,1,1), border=FALSE)
polygon( c(par('usr')[1]+ktx, par('usr')[1], par('usr')[1], par('usr')[1]+ktx), c(ylim[1]+kty, ylim[1]-kty*2, ylim[1], ylim[1]+kty*3), col=rgb(1,1,1), border=FALSE)

}
}
#########################

#########################
#  Legend  #
legens.ade<-function(x, y, vnames, g, p, xlim, ylim, lwd, lty=1, col, plot=TRUE, test){
n<- nlevels(g)
if(is.null(lty)) lty <- c(rep(1,n) ,0)
if(test)  letext<-c(paste(vnames,': ', levels(g)), paste('p Value: ', p))
if(!test) letext<-c(paste(vnames,': ', levels(g)))
lrect<-legend(legendon, inset = 0.0375,  legend=letext, plot=plot, lwd=c(rep(lwd+1, n),0), col = c(col, rgb(1,1,1,0)), bg=bgcol, box.lwd=2, text.col=tcol, lty = c(lty ,0),  merge = TRUE, yjust=0, box.col=rgb(1,1,1), text.width=max(strwidth(letext,font = 2)))
return(lrect)
}
#########################
}
################################################################################



#####  Style 2 #################################################################
if(wall==2){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

#########################
#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=1){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=rgb(1,1,1,0), axes=F)
if(is.null(xticks)) a1<-axis(1, col=rgb(1,1,1, 0), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(xlim, n=xticks),  col=rgb(1,1,1, 0), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, at=xticks,  col=rgb(1,1,1, 0), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1, 0), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=bgcol, lwd=1)

title(main, col=tcol)
box(col=a.coladd.ade(bgcol, -75))
if(ylim[1]>0.05 & ycut){
ktx<-diff(par('usr')[1:2])/75
kty<-diff(par('usr')[3:4])/100
par(xpd=TRUE)
segments( par('usr')[1], ylim[1]+kty*0.5, par('usr')[1], ylim[1]+kty*2.5 , col=rgb(1,1,1))
segments( par('usr')[1]-ktx, ylim[1]+kty, par('usr')[1]+ktx, ylim[1]+kty*4 , col=a.coladd.ade(bgcol, -50), lwd=1)
segments( par('usr')[1]-ktx, ylim[1]-kty, par('usr')[1]+ktx, ylim[1]+kty*2 , col=a.coladd.ade(bgcol, -50), lwd=1)
par(xpd=FALSE)
}
}
#########################

#########################
#  Legend  #
legens.ade<-function(x, y, vnames, g, p, xlim, ylim, lwd, lty=1, col, plot=TRUE, test){
n<- nlevels(g)
if(is.null(lty)) lty <- c(rep(1,n) ,0)
if(test)  letext<-c(paste(vnames,': ', levels(g)), paste('p Value: ', p))
if(!test) letext<-c(paste(vnames,': ', levels(g)))
lrect<-legend(legendon, inset = 0.0375,  legend=letext, plot=plot, lwd=c(rep(lwd+1, n),0), col = c(col, rgb(1,1,1,0)), bg=rgb(1,1,1), box.lwd=1, text.col=tcol, lty = c(lty ,0),  merge = TRUE, yjust=0, box.col=a.coladd.ade(bgcol, -75), text.width=max(strwidth(letext,font = 2)))
return(lrect)
}
#########################
}
################################################################################




#####  Style 3 #################################################################
if(wall==3){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

#########################
#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=1){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=rgb(1,1,1,0), axes=F)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(is.null(xticks)) a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -50), lwd.ticks=1)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(xlim, n=xticks),  col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -50), lwd.ticks=1)
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, at=xticks,  col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -50), lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -50), lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=a.coladd.ade(bgcol, -50), lwd=1)
title(main, col=tcol)
box(col=a.coladd.ade(bgcol, -75))
if(ylim[1]>0.05 & ycut){
ktx<-diff(par('usr')[1:2])/75
kty<-diff(par('usr')[3:4])/100
par(xpd=TRUE)
segments( par('usr')[1], ylim[1]+kty*0.5, par('usr')[1], ylim[1]+kty*2.5 , col=bgcol)
segments( par('usr')[1]-ktx, ylim[1]+kty, par('usr')[1]+ktx, ylim[1]+kty*4 , col=a.coladd.ade(bgcol, -50), lwd=1)
segments( par('usr')[1]-ktx, ylim[1]-kty, par('usr')[1]+ktx, ylim[1]+kty*2 , col=a.coladd.ade(bgcol, -50), lwd=1)
par(xpd=FALSE)
}
}
#########################

#########################
#  Legend  #
legens.ade<-function(x, y, vnames, g, p, xlim, ylim, lwd, lty=1, col, plot=TRUE, test){
n<- nlevels(g)
if(is.null(lty)) lty <- c(rep(1,n) ,0)
if(test)  letext<-c(paste(vnames,': ', levels(g)), paste('p Value: ', p))
if(!test) letext<-c(paste(vnames,': ', levels(g)))
lrect<-legend(legendon, inset = 0.0375,  legend=letext, plot=plot, lwd=c(rep(lwd+1, n),0), col = c(col, rgb(1,1,1,0)), bg=rgb(1,1,1), box.lwd=1, text.col=tcol, lty = c(lty ,0),  merge = TRUE, yjust=0, box.col=a.coladd.ade(bgcol, -75), text.width=max(strwidth(letext,font = 2)))
return(lrect)
}
#########################
}
################################################################################


#####  Style 4 #################################################################
if(wall==4){
par(col.axis=tcol)
par(col.lab=rgb(1,1,1))
par(col.main=rgb(1,1,1))
par(font=2)

#########################
#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=1){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=rgb(1,1,1,0), axes=F)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(is.null(xticks)) a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(xlim, n=xticks),  col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, at=xticks,  col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
par(xpd=TRUE)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0, 2.75,  2.75, 0)), col=tcol, border=rgb(1,1,1))
if(ylab!='' & ylab!=' ') polygon( a.glc(side=2, line=c(3.5, 3.5, 2, 2)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
if(xlab!='' & xlab!=' ') polygon( a.glc(side=c(2, 2, 4, 4), line=0),     a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))
text(a.glc(side=0), a.glc(side=3, line=1),    labels=main, cex = 1.25, font=2, col=rgb(1,1,1), adj=c(0.5,0))
text(a.glc(side=0), a.glc(side=1, line=3.5), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=2.5), a.glc(side=5),  labels=ylab, cex = 1.1,  font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)

box(col=rgb(1,1,1))
if(ylim[1]>0.05 & ycut){
ktx<-diff(par('usr')[1:2])/50
kty<-diff(par('usr')[3:4])/150
polygon( c(par('usr')[1]+ktx, par('usr')[1], par('usr')[1], par('usr')[1]+ktx), c(ylim[1]+kty*4, ylim[1]+kty, ylim[1]+kty*3, ylim[1]+kty*6), col=rgb(1,1,1), border=FALSE)
polygon( c(par('usr')[1]+ktx, par('usr')[1], par('usr')[1], par('usr')[1]+ktx), c(ylim[1]+kty, ylim[1]-kty*2, ylim[1], ylim[1]+kty*3), col=rgb(1,1,1), border=FALSE)

}
}
#########################

#########################
#  Legend  #
legens.ade<-function(x, y, vnames, g, p, xlim, ylim, lwd, lty=1, col, plot=TRUE,test){
n<- nlevels(g)
if(is.null(lty)) lty <- c(rep(1,n) ,0)
if(test)  letext<-c(paste(vnames,': ', levels(g)), paste('p Value: ', p))
if(!test) letext<-c(paste(vnames,': ', levels(g)))
rrect<- legend(legendon, inset = 0.0375,  legend=letext, lwd=c(rep(lwd+4, n),0)  ,col = c( rgb(1,1,1)) , bg=tcol, box.lwd=1,   lty = c(rep(1, n) ,0),  merge = TRUE, yjust=0, box.col=rgb(1,1,1), text.col=rgb(1,1,1), text.width=max(strwidth(letext,font = 2)))
lrect<- legend(legendon, inset = 0.0375,  legend=letext, lwd=c(rep(lwd+2, n),0)  ,col = c(col, rgb(1,1,1,0))         , bg=a.alpha.ade(tcol, 0), box.lwd=1, lty = c(lty ,0),  merge = TRUE, yjust=0, box.col=rgb(1,1,1), text.col=rgb(1,1,1), text.width=max(strwidth(letext,font = 2)))
return(lrect)
}
#########################
}
################################################################################



#####  Style 5 #################################################################
if(wall==5){
newmai<-rep(0, 4)
oldmai<-par('mai')

if(oldmai[2]<1) newmai[2]<- 1 - oldmai[2]
if(oldmai[3]>0.75 & oldmai[3]<=0.82) newmai[3]<- 0.75-oldmai[3]
if(oldmai[4]>0.25 & oldmai[4]<=0.42) newmai[4]<- 0.25-oldmai[4]
par(mai=(oldmai+newmai))


par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
par(font=2)

#########################
#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=1){
plot(0, 0, type='s', xlab='', ylab='', xlim=xlim, ylim=ylim, col=rgb(1,1,1,0), axes=F)
if(is.null(xticks)) a1<-axis(1, col=rgb(1,1,1, 0), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(xlim, n=xticks),  col=rgb(1,1,1, 0), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, at=xticks,  col=rgb(1,1,1, 0), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1, 0), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=bgcol, lwd=1)

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

box(col=tcol)
if(ylim[1]>0.05 & ycut) axis.break(axis=2,breakpos=ylim[1],bgcol="white",breakcol=tcol, style="slash", brw=0.02)
}
#########################

#########################
#  Legend  #
legens.ade<-function(x, y, vnames, g, p, xlim, ylim, lwd,  lty=1, col, plot=TRUE, test){
n<- nlevels(g)
if(is.null(lty)) lty <- c(rep(1,n) ,0)
if(test)  letext<-c(paste(vnames,': ', levels(g)), paste('p Value: ', p))
if(!test) letext<-c(paste(vnames,': ', levels(g)))
lrect<- legend(legendon, inset = 0.0375,  legend=letext, lwd=c(rep(lwd+1, n),0)  ,col = c(col, rgb(1,1,1,0)), bg=rgb(1,1,1),  box.lwd=1, lty = c(lty ,0),  merge = TRUE, yjust=0, box.col=tcol, text.col=tcol, text.width=max(strwidth(letext,font = 2)))
return(lrect)
}
#########################
}
################################################################################


#####  Style 6 #################################################################
if(wall==6){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

#########################
#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=1){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=rgb(1,1,1,0), axes=F)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)

if(is.null(xticks)) a1<-axis(1, col=rgb(1,1,1, 0), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
if(is.null(xticks)) a1<-axis(1, col=rgb(1,1,1, 0), col.ticks=rgb(1,1,1), lwd.ticks=1)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(xlim, n=xticks),  col=rgb(1,1,1, 0), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, at=pretty(xlim, n=xticks),  col=rgb(1,1,1, 0), col.ticks=rgb(1,1,1), lwd.ticks=1)
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, at=xticks,  col=rgb(1,1,1, 0), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, at=xticks,  col=rgb(1,1,1, 0), col.ticks=rgb(1,1,1), lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1, 0), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
a2<-axis(2, col=rgb(1,1,1, 0), col.ticks=rgb(1,1,1), lwd.ticks=1)

abline(v=a1, h=a2, lty=1, col=a.coladd.ade(bgcol, -35), lwd=3)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)


title(main, col=tcol)
box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=a.coladd.ade(bgcol, -35))


if(ylim[1]>0.05 & ycut){
ktx<-diff(par('usr')[1:2])/50
kty<-diff(par('usr')[3:4])/150
polygon( c(par('usr')[1]+ktx, par('usr')[1], par('usr')[1], par('usr')[1]+ktx), c(ylim[1]+kty*4, ylim[1]+kty, ylim[1]+kty*3, ylim[1]+kty*6), col=rgb(1,1,1,0), border=a.coladd.ade(bgcol, -35), lwd=2)
polygon( c(par('usr')[1]+ktx, par('usr')[1], par('usr')[1], par('usr')[1]+ktx), c(ylim[1]+kty, ylim[1]-kty*2, ylim[1], ylim[1]+kty*3), col=rgb(1,1,1,0), border=a.coladd.ade(bgcol, -35), lwd=2)

polygon( c(par('usr')[1]+ktx, par('usr')[1], par('usr')[1], par('usr')[1]+ktx), c(ylim[1]+kty*4, ylim[1]+kty, ylim[1]+kty*3, ylim[1]+kty*6), col=rgb(1,1,1), border=rgb(1,1,1,0))
polygon( c(par('usr')[1]+ktx, par('usr')[1], par('usr')[1], par('usr')[1]+ktx), c(ylim[1]+kty, ylim[1]-kty*2, ylim[1], ylim[1]+kty*3), col=rgb(1,1,1), border=rgb(1,1,1,0))


}
}
#########################

#########################
#  Legend  #
legens.ade<-function(x, y, vnames, g, p, xlim, ylim, lwd, lty=1, col, plot=TRUE,test){
n<- nlevels(g)
if(is.null(lty)) lty <- c(rep(1,n) ,0)
if(test)  letext<-c(paste(vnames,': ', levels(g)), paste('p Value: ', p))
if(!test) letext<-c(paste(vnames,': ', levels(g)))
lrect<-legend(legendon, inset = 0.0375, legend=letext, plot=plot, lwd=c(rep(lwd+1, n),0), col = c(col, rgb(1,1,1,0)), bg=bgcol, text.col=tcol,        lty = c(lty ,0),  merge = TRUE, yjust=0, box.col=rgb(1,1,1), box.lwd=3, text.width=max(strwidth(letext,font = 2)))
lrect<-legend(legendon, inset = 0.0375, legend=letext, plot=plot, lwd=c(rep(lwd+1, n),0), col = c(col, rgb(1,1,1,0)), bg=rgb(0,0,0,0), text.col=tcol, lty = c(lty ,0),  merge = TRUE, yjust=0, box.col=a.coladd.ade(bgcol, -35) , box.lwd=1, text.width=max(strwidth(letext,font = 2)))
return(lrect)
}
#########################
}
################################################################################


################################################################################
################################################################################
################################################################################

################################################################################
# 1 Gruppe
if(is.null(group)){
a.time<-eval(parse(text=paste('data$', time, sep='')))
a.event<-eval(parse(text=paste('data$', event, sep='')))
Sn <- survfit(Surv(a.time[!is.na(a.time) & !is.na(a.event)], a.event[!is.na(a.time) & !is.na(a.event)]) ~ 1, data=data, na.action=na.exclude)

if(ycut &  is.null(ylim))   ylim<- c(min(Sn$surv, na.rm=TRUE), 1)
if(!ycut)  ylim<- c(0, 1)

x.k<-sort(eval(parse(text=paste("data$",time, sep='')))[eval(parse(text=paste("data$",event, sep='')))==0])
y.k<-zenst(Sn$surv, sort(eval(parse(text=paste("data$",time, sep='')))[eval(parse(text=paste("data$",event, sep='')))==0]) ,Sn$time)

plot.box.ade( xlab, ylab, main, xlim=xlim, ylim=ylim , lwd=lwd)

a.line.ade(c(0, Sn$time, xrange[2]), c(1, Sn$surv, Sn$surv[length(Sn$surv)]), x.k, y.k,   c(1, Sn$upper),  c(1, Sn$lower), xlim=xlim, ylim=ylim, lwd=lwd, col=col ,CI=CI, zenspoints=zenspoints, lty=lty)
}
################################################################################


################################################################################
# Mehrere Gruppen
if(!is.null(group)){
if(nlevels(g)>1){
ylims<-NULL
Sn<-NULL
x.k<-NULL
y.k<-NULL
for(i in 1:nlevels(g)){
subdata <-  subset(data, eval(parse(text=group))==levels(g)[i])
Sn[[i]] <- eval(parse(text=paste("survfit(Surv(",time, "," ,event, ") ~ 1, data=subdata)", sep='')))

x.k[[i]]<-sort(eval(parse(text=paste("subdata$",time, sep='')))[eval(parse(text=paste("subdata$",event, sep='')))==0])
y.k[[i]]<-zenst(Sn[[i]]$surv, sort(eval(parse(text=paste("subdata$",time, sep='')))[eval(parse(text=paste("subdata$",event, sep='')))==0]) ,Sn[[i]]$time)
ylims<-min(c(Sn[[i]]$surv, ylims), na.rm=TRUE)
}
move=F
if(!ycut)  ylims<-0
if(is.null(ylim)){
ylim<-c(ylims,1)
move=T
}
plot.box.ade( xlab, ylab, main, xlim=xlim, ylim=ylim, lwd=lwd)

#################
#   p Value
testwerte<-survdiff(Surv(eval(parse(text=paste("data$",time, sep=''))) , eval(parse(text=paste("data$",event, sep='')))) ~ eval(parse(text=paste("data$",group, sep=''))), rho=0)
p<- (1 - pchisq(testwerte$chisq, (nlevels(g)-1)))
p <- format_p.ade(p, pdigs)
#################
lty<-rep(lty, nlevels(g))

for(i in 1:nlevels(g)){
survmy <-Sn[[i]]$surv
if(length(Sn[[i]]$time)-length(Sn[[i]]$surv)==1)   survmy<-c(1,survmy)
a.line.ade(c(0, Sn[[i]]$time, xrange[2]), c(1, survmy , Sn[[i]]$surv[length(Sn[[i]]$surv)]), x.k[[i]], y.k[[i]], c(1, Sn[[i]]$upper), c(1,Sn[[i]]$lower), xlim=xlim, ylim=c(0, 1), lwd=lwd, col=col[i], CI=CI, zenspoints=zenspoints, lty=lty[i])
}
zz<-legens.ade(0, ylim[1], vnames, g, p, xlim=xlim, ylim=ylim, lwd, lty=lty, col, plot=TRUE, test=test)

}
}
################################################################################
################################################################################



}
