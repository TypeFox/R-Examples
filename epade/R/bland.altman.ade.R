bland.altman.ade <-
function(x, y, data=NULL, ltext=TRUE, main='Bland-Altman Plot', xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, lwd=2, cex=1,  pch=16, lty=c(1,2,2), xticks=NULL, yticks=NULL,  col=NULL, tcol=NULL,  bgcol=NULL, lcol=c(4,2,2), alpha=NULL, fitline=1, wall=0, v=NULL, h=NULL, span=0.75){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr', 'plt', 'pin', 'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))
diag=FALSE

#####################################
# without Data.Frame or with
ismitdata<-FALSE
if(is.numeric(x)){
data<-NULL
data<-as.data.frame(x)
data$xmy<-x
data$ymy<-y
if(is.null(xlab))  xlab<-paste('Mean (', gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(x))), ', ', gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(y))),')', sep='')
if(is.null(ylab))  ylab<-paste(gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(x))), ' - ' , gsub('[(]{0}[A-Za-z0-9]*[$]', '' ,deparse(substitute(y))),sep='')

x<-'xmy'
y<-'ymy'
ismitdata=TRUE
}

if(!ismitdata){
if(!is.null(data)){ if(!is.data.frame(data))  stop("(data) must be a data.frame!") }
}
if(is.null(xlab))  xlab<-paste('Mean (', x, ', ', y,')', sep='')
if(is.null(ylab))  ylab<-paste(x, ' - ' , y,sep='')
#####################################

#####################################
# Errors

if(is.null(eval(parse(text=paste("data$",x)))))  stop('Variable x not found')
if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(lcol) & wall!=0)  lcol<- tcol
if(is.null(lcol) & wall==0)  lcol<- 1

################################################################################

a.x<-as.numeric(eval(parse(text=paste("data$",x))))
a.y<-as.numeric(eval(parse(text=paste("data$",y))))

a.diff<- (a.x-a.y)
a.mean<- (a.x+a.y)/2
dmean <- mean(a.diff, na.rm=T)
dsd <- sd(a.diff, na.rm=T)
h<-c(dmean, (dmean+1.96*dsd), (dmean-1.96*dsd), h)

if(is.null(ylim)){
ylim<-range(a.diff, na.rm=T)
if(ylim[1]>= dmean-1.96*dsd)  ylim[1]<- dmean-2*dsd
if(ylim[2]<= dmean+1.96*dsd)  ylim[2]<- dmean+2*dsd
}


scatter.ade(x=a.mean, y=a.diff, group=NULL, z=NULL, data=NULL, vnames=NULL, main=main, xlab=xlab, ylab=ylab, glab=NULL, zlab=NULL, xlim=xlim, ylim=ylim, zlim=NULL, lwd=lwd, cex=cex,  pch=pch,  lty=lty, xticks=xticks, yticks=yticks, zticks=NULL, col=col, tcol=tcol,  bgcol=bgcol, lcol=lcol, alpha=alpha, fitline=fitline, wall=wall, v=v, h=h, diag=diag, span=span)
if(!is.null(ltext)){
if(is.logical(ltext)){
if(ltext){
ltext<- c('Mean', '+1.96 SD', '-1.96 SD')
text(par('usr')[c(2,2,2)],  h, labels=ltext, adj=c(1.2,-0.5), col=tcol, font=2)
}
}
if(is.character(ltext)){
text(par('usr')[c(2,2,2)],  h, labels=ltext, adj=c(1.2,-0.5), col=tcol, font=2)
}


}
}
