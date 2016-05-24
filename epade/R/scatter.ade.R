scatter.ade <-
function(x, y=NULL, group=NULL, z=NULL, data=NULL, vnames=NULL, main=NULL, xlab=NULL, ylab=NULL, glab=NULL, zlab=NULL, legendon='topright', xlim=NULL, ylim=NULL, zlim=NULL, lwd=1, cex=1,  pch=16,  lty=1, xticks=NULL, yticks=NULL, zticks=NULL,  col=NULL, tcol=NULL,  bgcol=NULL, lcol=NULL, alpha=NULL, fitline=0, wall=0, v=NULL, h=NULL, diag=FALSE, span = 0.75){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr', 'plt',  'pin', 'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))



##############################
if(!is.character(x)){
xt<-deparse(substitute(x))
if(regexpr('~', xt)>=0){
y<-gsub('[~].*$', '', xt)
xpart<-gsub('^.*[~]', '', xt)
x<-gsub('[+].*$', '', xpart)
if(nchar(gsub('[^+]', '', xpart))==1) group<-gsub('^.*[+]', '', xpart)
if(nchar(gsub('[^+]', '', xpart))==2){
group<-gsub('^.*[+]', '', xpart)
z<-gsub('[+].*$', '', gsub(paste(x,'[+]', sep=''), '', xpart))
}
}}
##############################



vnames2<-NULL
if(!is.null(vnames) & is.list(vnames)){
vnames1 <-vnames[[1]]
vnames2 <-vnames[[2]]
vnames  <-vnames1
}



#####################################
# Nur X gegen ID
if(is.null(y)){
if(is.numeric(x)){
y<-1:length(x)
if(is.null(ylab))  ylab<-gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(x)))
}
if(is.character(x)){
data$my.id.fun<-1:dim(data)[1]
y<-'my.id.fun'
if(is.null(ylab)) ylab<-x
}
yret<-y
y<-x
x<-yret
if(is.null(xlab)) xlab<-'ID'
}

#####################################



#####################################
# without Data.Frame or with
ismitdata<-FALSE
if(is.numeric(x) & is.numeric(y) ){
data<-NULL
data<-as.data.frame(x)
data$ymy<-y
data$xmy<-x
data$zmy<-z
if(length(group)>3){
data$groupmy<-group
group<-'groupmy'
}
if(is.null(xlab))  xlab<- gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(x)))
if(is.null(ylab))  ylab<- gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(y)))
if(is.null(zlab))  zlab<- gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(z)))


x<-'xmy'
y<-'ymy'
if(!is.null(z)) z<-'zmy'
ismitdata=TRUE

}

if(!ismitdata){
if(!is.null(data)){ if(!is.data.frame(data))  stop("(data) must be a data.frame!") }
}

if(!is.null(group)){
g<-eval(parse(text=paste("data$",group, sep='')))
if(!is.factor(g)) g <- as.factor(g)
if(is.null(vnames)) vnames <- levels(g)
}
if(is.null(xlab))  xlab<-x
if(is.null(ylab))  ylab<-y
if(is.null(zlab))  zlab<-z
if(!is.null(zlab)) if(zlab=='' | zlab==' ' | zlab=='  ')  zlab<-NULL
#####################################

#####################################
# Errors

if(is.null(eval(parse(text=paste("data$",x)))))  stop('Variable x not found')
if(!is.null(y)) if(is.null(eval(parse(text=paste("data$",y)))))  stop('Variable y not found')
if(!is.null(z)) if(is.null(eval(parse(text=paste("data$",z)))))  stop('Variable z not found')
if(!is.null(group)) if(is.null(eval(parse(text=paste("data$",group)))))  stop('Variable gruop not found')



#####################################
# Z size
legendon2<-NULL

if(length(legendon)>1){
legendon2<-legendon[2]
legendon <-legendon[1]
}


if(!is.null(z)){
a.z<-eval(parse(text=paste("data$",z)))
a.z<-as.numeric(a.z)

cex<-a.z
zrange<-range(a.z, na.rm=TRUE)
cex<- (cex-min(cex, na.rm=TRUE))/max((cex-min(cex, na.rm=TRUE)), na.rm=TRUE)
if(is.null(zlim))   zlim<-c(0.75, 4)
cex<-(cex*zlim[2])+zlim[1]


nv2 <- length(unique(a.z))
if(nv2>3) nv2<-3

if(!is.null(zticks) & length(unique(a.z))>3) if(length(zticks)==1) nv2<-zticks

if(nv2==2){
if(is.null(vnames2))  vnames2<- c(format_n.ade(unique(a.z)))
z.ptcex  <- c(min(cex, na.rm=TRUE)+((max(cex, na.rm=TRUE)-min(cex, na.rm=TRUE))/5),  max(cex, na.rm=TRUE)-((max(cex, na.rm=TRUE)-min(cex, na.rm=TRUE))/5))
cex<- z.ptcex
}

if(nv2<=6 & nv2>2){
if(is.null(vnames2))  vnames2<- c(format_n.ade(quantile(a.z, seq(0.05, 0.95, length.out=nv2), na.rm=TRUE)))
z.ptcex  <- quantile(cex, seq(0.05, 0.95, length.out=nv2), na.rm=TRUE)
}

if(nv2>6){
if(is.null(vnames2))  vnames2<- c(format_n.ade(quantile(a.z, seq(0.01, 0.99, length.out=nv2), na.rm=TRUE)))
z.ptcex  <- quantile(cex, seq(0.01, 0.99, length.out=nv2), na.rm=TRUE)
}


if(length(zticks)>1){
vnames2<- zticks
z.ptcex  <- (zticks-min(a.z, na.rm=TRUE))/max((a.z-min(a.z, na.rm=TRUE)), na.rm=TRUE)
z.ptcex<-(z.ptcex*zlim[2])+zlim[1]
nv2<-length(zticks)
}



if(!is.null(group) & is.null(legendon2)){
if(legendon=='topleft')     legendon2<-'bottomleft'
if(legendon=='topright')    legendon2<-'bottomright'
if(legendon=='bottomleft')  legendon2<-'bottomright'
if(legendon=='bottomright') legendon2<-'bottomleft'
if(legendon=='bottom')      legendon2<-'top'
if(legendon=='top')         legendon2<-'bottom'
if(legendon=='left')        legendon2<-'right'
if(legendon=='right')       legendon2<-'left'
if(legendon=='center')      legendon2<-'top'
}
if(is.null(group))  legendon2<- legendon

}
#####################################


#####################################
# Colors
if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'

if(is.null(col) & wall!=0 & is.null(group)) col <- rgb(0.3,0.3,0.45)
if(is.null(col) & wall==0 & is.null(group)) col <- 'gray30'

if(!is.null(group)){
a.Ng<-nlevels(g)
if(is.null(col)) col<-a.getcol.ade(a.Ng)
}

if(is.null(lcol))  lcol<- tcol
#####################################


#####################################
# Alpha
if(is.null(alpha)){
a.N<-length(eval(parse(text=paste("data$",x))))
if(a.N<=1000) alpha<-1
if(a.N >1000){
alpha<- 1/(a.N/1500)
}
if(alpha<0.025) alpha<- 0.025
if(wall==0)    alpha<- 1
}
if(alpha>1)    alpha<- 1
#####################################





#####  Style 0 #################################################################
if(wall==0){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


######################
# Points
plot.ade<-function(x, y, xlim, ylim, cex=1, pch=16, col=1 , lwd=2, lty=1, alpha=1){

rgbc<-col2rgb(col, alpha = FALSE)
points(x, y, type='p', pch=pch, cex=cex, col=rgb(rgbc[1], rgbc[2], rgbc[3], alpha*255, maxColorValue=255))

# Lines
if(fitline==1 | fitline=='lm' | fitline=='lin' | fitline=='mean')   abline(lm(y~x), col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lwd=lwd, lty=lty)
if(fitline==2 | fitline=='lowess' | fitline=='lowes' | fitline=='low'){
loew<-loess(y[!is.na(x) & !is.na(y)]~x[!is.na(x) & !is.na(y)], span = span)
lines(loew$x[order(loew$x)],predict(loew)[order(loew$x)], col = rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lwd=lwd, lty=lty)
}
if(fitline==3 | fitline=='poly' | fitline=='polynom'){
myx<-seq(range(x, na.rm=TRUE)[1], range(x, na.rm=TRUE)[2], length.out=100)
newdat<-as.data.frame(cbind(myx, myx^2, myx^3, myx^4, myx^5, myx^6) )
names(newdat)<-c('x', 'x2', 'x3', 'x4', 'x5', 'x6')
pred<-predict(lm(y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6), na.action=na.exclude),newdat,interval="confidence")
points(newdat$x, pred[,1] , lwd=lwd, lty=lty, type='l', col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255))
}
}
######################


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
if(diag) abline(a=0, b=1, col=bgcol, lty=lty, lwd=lwd)
box( col=bgcol)
title(main)
}

#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
if(!is.null(group))  legend(legendon, legend=g,  pt.cex=2, title=glab, pch=c(rep(pch, length(g)))  ,lwd=c(rep(0, length(g)))  ,col = col, box.lwd=1, lty = c(rep(0, length(g))),  merge = TRUE, box.col=bgcol, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(c(g, glab),font = 2)))
if(!is.null(z))      legend(legendon2, y.intersp=1.75, x.intersp=1.25, legend=vnames2,  pt.cex=z.ptcex, pch=c(rep(pch, nv2))  ,lwd=c(rep(0, nv2))  ,col = tcol, box.lwd=1, lty = c(rep(0, nv2)),  merge = TRUE, box.col=bgcol, text.col=tcol, bg=rgb(1,1,1, 0), title=zlab, text.width=max(strwidth(c(vnames2, zlab),font = 2)))
}


}
################################################################################



#####  Style 1 #################################################################
if(wall==1){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

# Points
plot.ade<-function(x,y,xlim, ylim, cex=1, pch=16, col=1 , lwd=2, lty=1, alpha=1){

colb<- a.coladd.ade(col,  10)
cold<- a.coladd.ade(col, -50)
rgbc1<-col2rgb(colb, alpha = FALSE)
rgbc2<-col2rgb(cold, alpha = FALSE)
points(x, y, type='p', pch=pch, cex=cex, col=rgb(rgbc1[1], rgbc1[2], rgbc1[3], alpha*255, maxColorValue=255))
if(!is.null(z)) points(x, y, type='p', pch=pch-15  , cex=cex, col=a.coladd.ade(rgb(rgbc1[1], rgbc1[2], rgbc1[3], alpha*255, maxColorValue=255), -50))

rgbc<-rgbc2

# Lines
if(fitline==1 | fitline=='lm' | fitline=='lin' | fitline=='mean')   abline(lm(y~x), col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lty=lty, lwd=lwd)
if(fitline==2 | fitline=='lowess' | fitline=='lowes' | fitline=='low'){
loew<-loess(y[!is.na(x) & !is.na(y)]~x[!is.na(x) & !is.na(y)], span = span)
lines(loew$x[order(loew$x)],predict(loew)[order(loew$x)], col = rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lty=lty, lwd=lwd)
}
if(fitline==3 | fitline=='poly' | fitline=='polynom'){
myx<-seq(range(x, na.rm=TRUE)[1], range(x, na.rm=TRUE)[2], length.out=100)
newdat<-as.data.frame(cbind(myx, myx^2, myx^3, myx^4, myx^5, myx^6) )
names(newdat)<-c('x', 'x2', 'x3', 'x4', 'x5', 'x6')
pred<-predict(lm(y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6), na.action=na.exclude),newdat,interval="confidence")
points(newdat$x, pred[,1] , lwd=lwd, type='l', lty=lty, col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255))
}
}

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
if(diag) abline(a=0, b=1, col=tcol, lty=lty, lwd=lwd)
box(lwd=1, col=rgb(1,1,1))
title(main)
}

#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
n<- length(g)
if(!is.null(group)) legend(legendon, legend=g,  pt.cex=2,   title=glab,pch=c(rep(pch, n)) , col = a.alpha.ade(col, 1), lty = c(rep(0, n)), box.col=rgb(1,1,1) , box.lwd=2, text.col=tcol, bg=bgcol,  text.width=max(strwidth(c(g, glab),font = 2)))
if(!is.null(z))  legend(legendon2, y.intersp=1.75, x.intersp=1.25, legend=vnames2,  pt.cex=z.ptcex, pch=c(rep(pch, nv2))  ,lwd=c(rep(0, nv2))  ,col = tcol, box.lwd=2, lty = c(rep(0, nv2)),  merge = TRUE, box.col=rgb(1,1,1), text.col=tcol, bg=bgcol,title=zlab, text.width=max(strwidth(c(vnames2, zlab),font = 2)))
}
}
################################################################################


#####  Style 2 #################################################################
if(wall==2){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

# Points
plot.ade<-function(x,y,xlim, ylim, cex=1, pch=16, col=1 ,  lwd=2, lty=1, alpha=1){
colb<- a.coladd.ade(col, 10)
cold<- a.coladd.ade(col, -50)
rgbc1<-col2rgb(colb, alpha = FALSE)
rgbc2<-col2rgb(cold, alpha = FALSE)
points(x, y, type='p', pch=pch, cex=cex, col=rgb(rgbc1[1], rgbc1[2], rgbc1[3], alpha*255, maxColorValue=255))
points(x, y, type='p', pch=pch-15  , cex=cex, col=a.coladd.ade(rgb(rgbc1[1], rgbc1[2], rgbc1[3], alpha*255, maxColorValue=255), -25))
rgbc<-rgbc2
# Lines
if(fitline==1 | fitline=='lm' | fitline=='lin' | fitline=='mean')   abline(lm(y~x), col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lty=lty, lwd=lwd)
if(fitline==2 | fitline=='lowess' | fitline=='lowes' | fitline=='low'){
loew<-loess(y[!is.na(x) & !is.na(y)]~x[!is.na(x) & !is.na(y)], span = span)
lines(loew$x[order(loew$x)],predict(loew)[order(loew$x)], col = rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lty=lty, lwd=lwd)
}
if(fitline==3 | fitline=='poly' | fitline=='polynom'){
myx<-seq(range(x, na.rm=TRUE)[1], range(x, na.rm=TRUE)[2], length.out=100)
newdat<-as.data.frame(cbind(myx, myx^2, myx^3, myx^4, myx^5, myx^6) )
names(newdat)<-c('x', 'x2', 'x3', 'x4', 'x5', 'x6')
pred<-predict(lm(y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6), na.action=na.exclude),newdat,interval="confidence")
points(newdat$x, pred[,1] , lwd=lwd, type='l', lty=lty, col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255))
}
}

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
if(diag) abline(a=0, b=1, col=a.coladd.ade(bgcol, -75), lty=lty, lwd=lwd)
box(col=a.coladd.ade(bgcol, -75))
title(main)
}


#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
n<- length(g)
if(!is.null(group)) legend(legendon, legend=g,  pt.cex=2,  title=glab, pch=c(rep(pch, n)) ,col = a.alpha.ade(a.coladd.ade(col, 10), 1), lty = c(rep(0, n)), box.col=a.coladd.ade(bgcol, -75) , box.lwd=1, text.col=tcol, bg=rgb(1,1,1),  text.width=max(strwidth(c(g, glab),font = 2)))
if(!is.null(group)) legend(legendon, legend=g,  pt.cex=2,  title=glab, bg=rgb(1,1,1,0),    pch=c(rep(pch-15, n)) , col = a.alpha.ade(a.coladd.ade(col, -50), 1), lty = c(rep(0, n)), box.col=a.coladd.ade(bgcol, -75) , box.lwd=0, text.col=tcol,  text.width=max(strwidth(c(g, glab),font = 2)))
if(!is.null(z))  legend(legendon2, y.intersp=1.75, x.intersp=1.25, legend=vnames2,  pt.cex=z.ptcex, pch=c(rep(pch, nv2))  ,lwd=c(rep(0, nv2))  ,col = tcol, box.lwd=1, lty = c(rep(0, nv2)),  merge = TRUE, box.col=a.coladd.ade(bgcol, -75), text.col=tcol, title=zlab, bg=rgb(1,1,1), text.width=max(strwidth(c(vnames2, zlab),font = 2)))
}
}
################################################################################


#####  Style 3 #################################################################
if(wall==3){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

# Points
plot.ade<-function(x,y,xlim, ylim, cex=1, pch=16, col=1 , lwd=2, lty=1, alpha=1){
colb<- a.coladd.ade(col, 10)
cold<- a.coladd.ade(col, -50)
rgbc1<-col2rgb(colb, alpha = FALSE)
rgbc2<-col2rgb(cold, alpha = FALSE)
points(x, y, type='p', pch=pch, cex=cex, col=rgb(rgbc1[1], rgbc1[2], rgbc1[3], alpha*255, maxColorValue=255))
points(x, y, type='p', pch=pch-15  , cex=cex, col=a.coladd.ade(rgb(rgbc1[1], rgbc1[2], rgbc1[3], alpha*255, maxColorValue=255), -25))
rgbc<-rgbc2
# Lines
if(fitline==1 | fitline=='lm' | fitline=='lin' | fitline=='mean')   abline(lm(y~x), col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lwd=lwd, lty=lty)
if(fitline==2 | fitline=='lowess' | fitline=='lowes' | fitline=='low'){
loew<-loess(y[!is.na(x) & !is.na(y)]~x[!is.na(x) & !is.na(y)], span = span)
lines(loew$x[order(loew$x)],predict(loew)[order(loew$x)], col = rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lwd=lwd, lty=lty)
}
if(fitline==3 | fitline=='poly' | fitline=='polynom'){
myx<-seq(range(x, na.rm=TRUE)[1], range(x, na.rm=TRUE)[2], length.out=100)
newdat<-as.data.frame(cbind(myx, myx^2, myx^3, myx^4, myx^5, myx^6) )
names(newdat)<-c('x', 'x2', 'x3', 'x4', 'x5', 'x6')
pred<-predict(lm(y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6), na.action=na.exclude),newdat,interval="confidence")
points(newdat$x, pred[,1] , lwd=lwd, lty=lty,  type='l', col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255))
}
}

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
if(diag) abline(a=0, b=1, col=tcol, lty=lty, lwd=lwd)
box(lwd=1, col=a.coladd.ade(bgcol, -50))
title(main)
}


#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
n<- length(g)
if(!is.null(group)) legend(legendon, legend=g,  pt.cex=2, title=glab,  pch=c(rep(pch, n)) ,col = a.alpha.ade(a.coladd.ade(col, 10), 1), lty = c(rep(0, n)), box.col=a.coladd.ade(bgcol, -50) , box.lwd=1, text.col=tcol, bg=rgb(1,1,1),  text.width=max(strwidth(c(g, glab),font = 2)))
if(!is.null(group)) legend(legendon, legend=g,  pt.cex=2, title=glab,  bg=rgb(1,1,1,0),  pch=c(rep(pch-15, n)) , col = a.alpha.ade(a.coladd.ade(col, -50), 1), lty = c(rep(0, n)), box.col=a.coladd.ade(bgcol, -50) , box.lwd=0, text.col=tcol,  text.width=max(strwidth(c(g, glab),font = 2)))
if(!is.null(z))  legend(legendon2, y.intersp=1.75, x.intersp=1.25, legend=vnames2,  pt.cex=z.ptcex, pch=c(rep(pch, nv2))  ,lwd=c(rep(0, nv2))  ,col = tcol, box.lwd=1, lty = c(rep(0, nv2)),  merge = TRUE, box.col=a.coladd.ade(bgcol, -50), text.col=tcol, title=zlab, bg=rgb(1,1,1), text.width=max(strwidth(c(vnames2, zlab),font = 2)))
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

# Points
plot.ade<-function(x,y,xlim, ylim, cex=1, pch=16, col=1 , lwd=2, lty=1, alpha=1){
colb<- a.coladd.ade(col, 10)
cold<- a.coladd.ade(col, -50)
rgbc1<-col2rgb(colb, alpha = FALSE)
rgbc2<-col2rgb(cold, alpha = FALSE)
points(x, y, type='p', pch=pch, cex=cex, col=rgb(rgbc1[1], rgbc1[2], rgbc1[3], alpha*255, maxColorValue=255))
if(!is.null(z)) points(x, y, type='p', pch=pch-15  , cex=cex, col=a.coladd.ade(rgb(rgbc1[1], rgbc1[2], rgbc1[3], alpha*255, maxColorValue=255), -50))
rgbc<-rgbc2
# Lines
if(fitline==1 | fitline=='lm' | fitline=='lin' | fitline=='mean')   abline(lm(y~x), col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lwd=lwd, lty=lty)
if(fitline==2 | fitline=='lowess' | fitline=='lowes' | fitline=='low'){
loew<-loess(y[!is.na(x) & !is.na(y)]~x[!is.na(x) & !is.na(y)], span = span)
lines(loew$x[order(loew$x)],predict(loew)[order(loew$x)], col = rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lwd=lwd, lty=lty)
}
if(fitline==3 | fitline=='poly' | fitline=='polynom'){
myx<-seq(range(x, na.rm=TRUE)[1], range(x, na.rm=TRUE)[2], length.out=100)
newdat<-as.data.frame(cbind(myx, myx^2, myx^3, myx^4, myx^5, myx^6) )
names(newdat)<-c('x', 'x2', 'x3', 'x4', 'x5', 'x6')
pred<-predict(lm(y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6), na.action=na.exclude),newdat,interval="confidence")
points(newdat$x, pred[,1] , lwd=lwd, lty=lty, type='l', col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255))
}
}

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
if(diag) abline(a=0, b=1, col=tcol, lty=lty, lwd=lwd)
# Outer
par(xpd=TRUE)
polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0, 2.75,  2.75, 0)), col=tcol, border=rgb(1,1,1))
if(ylab!='' & ylab!=' ') polygon( a.glc(side=2, line=c(3.5, 3.5, 2, 2)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
if(xlab!='' & xlab!=' ') polygon( a.glc(side=c(2, 2, 4, 4), line=0),     a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))
text(a.glc(side=0), a.glc(side=3, line=1),    labels=main, cex = 1.25, font=2, col=rgb(1,1,1), adj=c(0.5,0))
text(a.glc(side=0), a.glc(side=1, line=3.5), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=2.5), a.glc(side=5),  labels=ylab, cex = 1.1,  font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)
box(lwd=1, col=rgb(1,1,1))

}



#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
n<- length(g)
if(!is.null(group)) legend(legendon, legend=g,  pt.cex=2, title=glab, pch=c(rep(pch, n)) ,   col = a.alpha.ade(col, 1), lty = c(rep(0, n)), box.col=rgb(1,1,1) , box.lwd=1, text.col=rgb(1,1,1), bg=tcol,  text.width=max(strwidth(c(g, glab),font = 2)))
if(!is.null(group)) legend(legendon, legend=g,  pt.cex=2, title=glab, pch=c(rep(pch-15, n)) ,col = rgb(1,1,1), lty = c(rep(0, n)), box.col=rgb(1,1,1) , box.lwd=1, text.col=rgb(1,1,1), bg=a.alpha.ade(col, 0),  text.width=max(strwidth(c(g, glab),font = 2)))
if(!is.null(z))  legend(legendon2, y.intersp=1.75, x.intersp=1.25, legend=vnames2,  pt.cex=z.ptcex, pch=c(rep(pch, nv2))  ,lwd=c(rep(0, nv2))  ,col = rgb(1,1,1), box.lwd=1, lty = c(rep(0, nv2)),  merge = TRUE, box.col=rgb(1,1,1), text.col=rgb(1,1,1), title=zlab, bg=tcol, text.width=max(strwidth(c(vnames2, zlab),font = 2)))

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
#if(oldmai[1]<1.2) newmai[1]<- 1.2 - oldmai[1]
if(oldmai[2]<1) newmai[2]<- 1 - oldmai[2]
if(oldmai[3]>0.75 & oldmai[3]<=0.82) newmai[3]<- 0.75-oldmai[3]
if(oldmai[4]>0.25 & oldmai[4]<=0.42) newmai[4]<- 0.25-oldmai[4]
par(mai=(oldmai+newmai))



# Points
plot.ade<-function(x,y,xlim, ylim, cex=1, pch=16, col=1 , lwd=2, lty=1, alpha=1){
if(is.null(alpha)){ alpha<- 1/(length(x)/500)
if(alpha>1) alpha<- 1
if(alpha<0.25) alpha<- 0.25 }

colb<- a.coladd.ade(col, 10)
cold<- a.coladd.ade(col, -50)
rgbc1<-col2rgb(colb, alpha = FALSE)
rgbc2<-col2rgb(cold, alpha = FALSE)
points(x, y, type='p', pch=pch, cex=cex, col=rgb(rgbc1[1], rgbc1[2], rgbc1[3], alpha*255, maxColorValue=255))
if(!is.null(z)) points(x, y, type='p', pch=pch-15  , cex=cex, col=a.coladd.ade(rgb(rgbc1[1], rgbc1[2], rgbc1[3], alpha*255, maxColorValue=255), -50))

rgbc<-rgbc2
# Lines
if(fitline==1 | fitline=='lm' | fitline=='lin' | fitline=='mean')   abline(lm(y~x), col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lwd=lwd, lty=lty)
if(fitline==2 | fitline=='lowess' | fitline=='lowes' | fitline=='low'){
loew<-loess(y[!is.na(x) & !is.na(y)]~x[!is.na(x) & !is.na(y)], span = span)
lines(loew$x[order(loew$x)],predict(loew)[order(loew$x)], col = rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lwd=lwd, lty=lty)
}
if(fitline==3 | fitline=='poly' | fitline=='polynom'){
myx<-seq(range(x, na.rm=TRUE)[1], range(x, na.rm=TRUE)[2], length.out=100)
newdat<-as.data.frame(cbind(myx, myx^2, myx^3, myx^4, myx^5, myx^6) )
names(newdat)<-c('x', 'x2', 'x3', 'x4', 'x5', 'x6')
pred<-predict(lm(y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6), na.action=na.exclude),newdat,interval="confidence")
points(newdat$x, pred[,1] , lwd=lwd, lty=lty , type='l', col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255))
}
}

#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=3){
plot(0, 0, type='s', xlab='', ylab='', xlim=xlim, ylim=ylim ,  axes=FALSE)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks), labels=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks))
if(!is.null(yticks) & length(yticks)==1) a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks), labels=pretty(seq(ylim[1], ylim[2], length.out =100), n = yticks))
if(!is.null(xticks) & length(xticks)>1)  a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=xticks, labels=xticks)
if(!is.null(yticks) & length(yticks)>1)  a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1, at=yticks, labels=yticks)
if(is.null(xticks)) a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
if(is.null(yticks)) a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
#abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
abline(v=v, h=h, col=lcol, lty=lty, lwd=lwd)
if(diag) abline(a=0, b=1, col=tcol, lty=lty, lwd=lwd)
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

if(!is.null(group)) legend(legendon, legend=g,  pt.cex=2, title=glab, pch=c(rep(pch, n)) ,  col = a.alpha.ade(col, 1), lty = c(rep(0, n)), box.col=tcol , box.lwd=1, text.col=tcol, bg=rgb(1,1,1),  text.width=max(strwidth(c(g, glab),font = 2)))
if(!is.null(group)) legend(legendon, legend=g,  pt.cex=2, title=glab,pch=c(rep(pch-15, n)) ,col = tcol, lty = c(rep(0, n)), box.col=tcol , box.lwd=1, text.col=tcol, bg=rgb(1,1,1,0),  text.width=max(strwidth(c(g, glab),font = 2)))
if(!is.null(z))  legend(legendon2, y.intersp=1.75, x.intersp=1.25, legend=vnames2,  pt.cex=z.ptcex, pch=c(rep(pch, nv2))  ,lwd=c(rep(0, nv2))  ,col = tcol, box.lwd=1, lty = c(rep(0, nv2)),  merge = TRUE, box.col=tcol, text.col=tcol, bg=rgb(1,1,1), title=zlab, text.width=max(strwidth(c(vnames2, zlab),font = 2)))
box(lwd=1, col=tcol)
}

}
################################################################################



#####  Style 6 #################################################################
if(wall==6){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)

# Points
plot.ade<-function(x,y,xlim, ylim, cex=1, pch=16, col=1 , lwd=2, lty=1, alpha=1){

colb<- a.coladd.ade(col,  10)
cold<- a.coladd.ade(col, -50)
rgbc1<-col2rgb(colb, alpha = FALSE)
rgbc2<-col2rgb(cold, alpha = FALSE)
points(x, y, type='p', pch=pch, cex=cex, col=rgb(rgbc1[1], rgbc1[2], rgbc1[3], alpha*255, maxColorValue=255))
if(!is.null(z)) points(x, y, type='p', pch=pch-15  , cex=cex, col=a.coladd.ade(rgb(rgbc1[1], rgbc1[2], rgbc1[3], alpha*255, maxColorValue=255), -50))

rgbc<-rgbc2

# Lines
if(fitline==1 | fitline=='lm' | fitline=='lin' | fitline=='mean')   abline(lm(y~x), col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lty=lty, lwd=lwd)
if(fitline==2 | fitline=='lowess' | fitline=='lowes' | fitline=='low'){
loew<-loess(y[!is.na(x) & !is.na(y)]~x[!is.na(x) & !is.na(y)], span = span)
lines(loew$x[order(loew$x)],predict(loew)[order(loew$x)], col = rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255), lty=lty, lwd=lwd)
}
if(fitline==3 | fitline=='poly' | fitline=='polynom'){
myx<-seq(range(x, na.rm=TRUE)[1], range(x, na.rm=TRUE)[2], length.out=100)
newdat<-as.data.frame(cbind(myx, myx^2, myx^3, myx^4, myx^5, myx^6) )
names(newdat)<-c('x', 'x2', 'x3', 'x4', 'x5', 'x6')
pred<-predict(lm(y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6), na.action=na.exclude),newdat,interval="confidence")
points(newdat$x, pred[,1] , lwd=lwd, type='l', lty=lty, col=rgb(rgbc[1], rgbc[2], rgbc[3], 255, maxColorValue=255))
}
}

#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd=3){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim ,  axes=FALSE)
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3, at=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks), labels=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks))
if(!is.null(xticks) & length(xticks)==1) a1<-axis(1, col=rgb(1,1,1), col.ticks=rgb(1,1,1),               lwd.ticks=1, at=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks), labels=pretty(seq(xlim[1], xlim[2], length.out =100), n = xticks))
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
if(diag) abline(a=0, b=1, col=tcol, lty=lty, lwd=lwd)
title(main)
box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=a.coladd.ade(bgcol, -35))
}

#  Legend  #
legens.ade<-function(ylims, g, xlim, ylim, lwd=3, pch, col, legendon){
n<- length(g)
if(!is.null(group)) legend(legendon, legend=g,  pt.cex=2,   title=glab,pch=c(rep(pch, n)) , col = a.alpha.ade(col, 1), lty = c(rep(0, n)), box.col=rgb(1,1,1) , box.lwd=3, text.col=tcol, bg=bgcol,  text.width=max(strwidth(c(g, glab),font = 2)))
if(!is.null(group)) legend(legendon, legend=g,  pt.cex=2,   title=glab,pch=c(rep(pch, n)) , col = a.alpha.ade(col, 1), lty = c(rep(0, n)), box.col=a.coladd.ade(bgcol, -35) ,               box.lwd=1, text.col=tcol, bg=rgb(1,1,1,0),  text.width=max(strwidth(c(g, glab),font = 2)))
if(!is.null(z))     legend(legendon2, y.intersp=1.75, x.intersp=1.25, legend=vnames2,  pt.cex=z.ptcex, pch=c(rep(pch, nv2))  ,lwd=c(rep(0, nv2))  ,col = tcol, box.lwd=3, lty = c(rep(0, nv2)),  merge = TRUE, box.col=rgb(1,1,1),              text.col=tcol, bg=bgcol,title=zlab, text.width=max(strwidth(c(vnames2, zlab),font = 2)))
if(!is.null(z))     legend(legendon2, y.intersp=1.75, x.intersp=1.25, legend=vnames2,  pt.cex=z.ptcex, pch=c(rep(pch, nv2))  ,lwd=c(rep(0, nv2))  ,col = tcol, box.lwd=1, lty = c(rep(0, nv2)),  merge = TRUE, box.col=a.coladd.ade(bgcol, -35), text.col=tcol, bg=rgb(1,1,1,0),title=zlab, text.width=max(strwidth(c(vnames2, zlab),font = 2)))
}
}
################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


################################################################################
# Plotaufruf
a.x<-as.numeric(eval(parse(text=paste("data$",x))))
a.y<-as.numeric(eval(parse(text=paste("data$",y))))

if(is.null(group) &  is.null(z)){
if(is.null(xlim)) xlim<- range(a.x, na.rm=TRUE)
if(is.null(ylim)) ylim<- range(a.y, na.rm=TRUE)
}

if(!is.null(group) &  is.null(z)){
if(is.null(xlim)) xlim<- range(a.x[!is.na(g)], na.rm=TRUE)
if(is.null(ylim)) ylim<- range(a.y[!is.na(g)], na.rm=TRUE)
}

if(!is.null(group) &  !is.null(z)){
if(is.null(xlim)) xlim<- range(a.x[!is.na(g) & !is.na(cex)], na.rm=TRUE)
if(is.null(ylim)) ylim<- range(a.y[!is.na(g) & !is.na(cex)], na.rm=TRUE)
}

if(is.null(group) &  !is.null(z)){
if(is.null(xlim)) xlim<- range(a.x[!is.na(cex)], na.rm=TRUE)
if(is.null(ylim)) ylim<- range(a.y[!is.na(cex)], na.rm=TRUE)
}



plot.box.ade( xlab=xlab, ylab=ylab, main=main, xlim=xlim, ylim=ylim, lwd=lwd)

if(is.null(group)){
plot.ade(a.x, a.y, xlim=xlim, ylim=ylim, cex=cex, pch=pch, lwd=lwd, col=col, lty=lty, alpha=alpha)
}


if(!is.null(group) | !is.null(z)){
if(legendon[1]!='none') legens.ade(ylim[1], vnames, xlim=xlim, ylim=ylim, lwd=lwd, pch=pch, col=col, legendon=legendon)
}

if(!is.null(group)){
for(i in 1:nlevels(g)){
subdat<-  subset(data, subset=(eval(parse(text=group))==levels(g)[i]))
if(length(cex)==length(g)){
if(length(pch)==nlevels(g)) plot.ade(eval(parse(text=paste("subdat$",x))), eval(parse(text=paste("subdat$",y))), xlim=xlim, ylim=ylim, lwd=lwd, cex=cex[which(g==levels(g)[i])], pch=pch[i], col=col[i], lty=lty, alpha )
if(length(pch)!=nlevels(g)) plot.ade(eval(parse(text=paste("subdat$",x))), eval(parse(text=paste("subdat$",y))), xlim=xlim, ylim=ylim, lwd=lwd, cex=cex[which(g==levels(g)[i])], pch=pch,    col=col[i], lty=lty, alpha )
}
if(length(cex)<length(g)){
if(length(pch)==nlevels(g)) plot.ade(eval(parse(text=paste("subdat$",x))), eval(parse(text=paste("subdat$",y))), xlim=xlim, ylim=ylim, lwd=lwd, cex=cex, pch=pch[i], col=col[i], lty=lty, alpha )
if(length(pch)!=nlevels(g)) plot.ade(eval(parse(text=paste("subdat$",x))), eval(parse(text=paste("subdat$",y))), xlim=xlim, ylim=ylim, lwd=lwd, cex=cex, pch=pch,    col=col[i], lty=lty, alpha )
}

}
}

################################################################################         



}
