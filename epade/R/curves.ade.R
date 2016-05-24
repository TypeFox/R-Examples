curves.ade <-
function( x, y=NULL, group=NULL, data=NULL, vnames=NULL, main=NULL, xlab=NULL, ylab=NULL, legendon='topright', xlim=NULL, ylim=NULL, lwd=1, lwd2=1, cex=1,  pch=16,  lty=1, lty2=2, col=NULL, xticks=NULL, yticks=NULL, tcol=NULL,  bgcol=NULL, alpha=NULL, fitline=0, wall=0, v=NULL, h=NULL, diag=F, points=T){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr', 'plt',   'pin',  'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))

##############################
if(!is.character(x)){
xt<-deparse(substitute(x))
if(regexpr('~', xt)>=0){
y<-gsub('[~].*$', '', xt)
xpart<-gsub('^.*[~]', '', xt)
x<-gsub('[+].*$', '', xpart)
if(nchar(gsub('[^+]', '', xpart))==1) group<-gsub('^.*[+]', '', xpart)
}}
##############################


library(plotrix)
g<-group
if(is.data.frame(data)  &  is.character(group) &   !is.null(group) )   g<-eval(parse(text=paste("data$",group, sep='')))
if(!is.data.frame(data) & !is.character(group))    g<-group
if(is.data.frame(data) & is.character(x)){
x<-eval(parse(text=paste("data$",x, sep='')))
y<-eval(parse(text=paste("data$",y, sep='')))
}
if(!is.null(g)) g<-as.factor(g)

if(is.null(col) &  is.null(g)) col <- rgb(0.05,0.05,0.175)
if(!is.null(g)){
if(is.null(col)) col <- a.getcol.ade(nlevels(g))
}
lty<-rep(lty, 16)
if(is.null(alpha)){
alpha<-0.5
if(wall==0) alpha<-1
}

xlim2<- xlim

o<-order(x)
x<-x[o]
y<-y[o]
if(!is.null(g)) g<-g[o]

if(points) alpha2<-alpha
if(!points) alpha2<-0

scatter.ade(x=x, y=y, group=g, data=data, vnames=vnames, main=main, xlab=xlab, ylab=ylab, legendon=legendon, xlim=xlim2, ylim=ylim, lwd=lwd2, cex=cex,  pch=pch,  lty=lty2, col=col, xticks=xticks, yticks=yticks, tcol=tcol,  bgcol=bgcol, alpha=alpha2, fitline=fitline, wall=wall, v=v, h=h,  diag=diag)
if(is.null(g)) lines(x, y , type = "l", col=col, lty=lty, lwd=lwd)
for(k in 1:nlevels(g)){
lines(x[g==levels(g)[k]], y[g==levels(g)[k]] , type = "l", col=col[k], lty=lty[k], lwd=lwd)
}

################################################################################         
################################################################################
}
