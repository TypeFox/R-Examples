qq.ade <-
function(x, data=NULL, main='Q-Q Plot', xlab='Theoretical Quantiles', ylab='Sample Quantiles', xlim=NULL, ylim=NULL, lwd=1, cex=1,  pch=16, lty=1, xticks=NULL, yticks=NULL,  col=NULL, tcol=NULL,  bgcol=NULL, lcol=NULL, alpha=NULL, fitline=0, qline=TRUE, wall=0, v=NULL, h=NULL, diag=FALSE, band=FALSE, span=0.75){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr', 'plt',  'pin', 'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))


#####################################
# without Data.Frame or with
ismitdata<-FALSE
if(is.numeric(x)){
data<-NULL
data<-as.data.frame(x)
data$xmy<-x
if(is.null(xlab))  xlab<-gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(x)))
if(is.null(ylab))  ylab<-gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(y)))

x<-'xmy'
y<-'ymy'
ismitdata=TRUE
}

if(!ismitdata){
if(!is.null(data)){ if(!is.data.frame(data))  stop("(data) must be a data.frame!") }
}
if(is.null(xlab))  xlab<-x
#####################################

#####################################
# Errors

if(is.null(eval(parse(text=paste("data$",x)))))  stop('Variable x not found')
if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(lcol) & wall!=0)  lcol<- tcol
if(is.null(lcol) & wall==0)  lcol<- 1


a.y<-as.numeric(eval(parse(text=paste("data$",x))))

a.y <- sort(a.y)
a.x <- qnorm(ppoints(length(a.y)))[order(order(a.y))]


scatter.ade(a.x, a.y, group=NULL, z=NULL, data=NULL, vnames=NULL, main=main, xlab=xlab, ylab=ylab, glab=NULL, zlab=NULL, xlim=xlim, ylim=ylim, zlim=NULL, lwd=lwd, cex=cex,  pch=pch,  lty=lty, xticks=xticks, yticks=yticks, zticks=NULL, col=col, tcol=tcol,  bgcol=bgcol, lcol=lcol, alpha=alpha, fitline=fitline, wall=wall, v=v, h=h, diag=diag, span=span)


Nrun<-1000
if(is.numeric(band)){
Nrun<-band
band<-TRUE
}
if(is.logical(band)){
if(band){
a.m  <-mean(a.y, na.rm=TRUE)
a.sd <-sd(a.y, na.rm=TRUE)
a.n <-length(a.y)
ymax<-NULL
ymin<-NULL
for (i in 1:Nrun) {
y<-rnorm(a.n, a.m, a.sd)
y <- sort(y)
x<-qnorm(ppoints(length(y)))[order(order(y))]
if(is.null(ymax)) ymax<-y
if(is.null(ymin)) ymin<-y
ymax[y>ymax]<-y[y>ymax]
ymin[y<ymin]<-y[y<ymin]
}
points(x, ymax, col=lcol,  type='l', lty=2, lwd=2)
points(x, ymin, col=lcol,  type='l', lty=2, lwd=2)
}
}
y <- quantile(a.y[!is.na(a.y)], c(0.25, 0.75))
x <- qnorm(c(0.25, 0.75))
slope <- diff(y)/diff(x)
int <- y[1L] - slope * x[1L]
if(fitline==0 & qline) abline(int, slope, lwd=lwd, lty=lty, col=lcol)


}
