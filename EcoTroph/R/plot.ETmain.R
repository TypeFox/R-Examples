plot.ETmain <-
function(x,scale1=NULL,scale2=NULL,scale3=NULL,legend.cex=NULL,ask=interactive(),...){
m<-x
if (is.null(legend.cex)) {legend.cex=0.8}
par(ask=ask)
options(warn=-1)
#the scale parameter can be log or not for the y axis
#library(RColorBrewer)

plot(m$biomass,title="Biomasses",scale=scale1,legend.cex=legend.cex,...)
plot(m$biomass_acc,title="Accessible Biomasses",legend.cex=legend.cex,scale=scale2,...)

for (pecheries in names(m$Y))
{
plot(m$Y[[paste(pecheries)]],title=paste(pecheries),scale=scale3,legend.cex=legend.cex,...)
}

##Plot "Total catch"

opar=par(no.readonly=TRUE)
par(mfrow=c(2,1))
par(mar=c(5,3,1,15))

m1<-m$ET_Main[,c("B","B_acc","Y_tot","P","Kin")]
y_max <- max(m1[,"Y_tot"])
y_min <- 0
plot(rownames(m1),m1[,"Y_tot"],col=1,axes=F, type='l',lwd=2,ylab=NA,xlab="TL",main="Fishing fleets cumulated catches ",ylim=c(y_min,y_max),xaxs = "i",xlim=c(2,5),...)

p <- 1
m_precedente <- rep(y_min,length(rownames(m1)))
m_cum <- m_precedente
for (pecheries in names(m$Y))
{
p <- p+1
toto <- m$Y[[paste(pecheries)]]
m_cum <- m_cum+ apply(toto,1,sum)
polygon(c(as.numeric(names(m_cum)),rev(as.numeric(names(m_cum)))),c(m_cum,rev(m_precedente)),col=p,xlim=c(2,5),border='NA')
m_precedente<-m_cum
}
axis(4,...)
legend(5.7,y_max,legend = names(m$Y), bg = 'gray90',col=c(2:p),pch=1,xpd=NA,cex=legend.cex)

plot(rownames(m1),m1[,"Y_tot"],col=2,axes=F, type='l',lwd=2,ylab=NA,xlab="TL",main="Total Catch Trophic Spectra",ylim=c(y_min,y_max),xaxs = "i",xlim=c(2,5),...)
axis(1,...)
axis(4,...)
par(opar)

##Summary plots
opar=par(no.readonly=TRUE)
par(mfrow=c(2,1))
par(mar=c(5,5,3,10))

plot(rownames(m1),m1[,"B"],log='y',col=2,type='l',ylim=c(0.0000004,max(m1)),xlim=c(2,5),xlab="TL",ylab=NA,lwd=2,main="Trophic spectra summary plots",...)
lines(rownames(m1),m1[,"B_acc"],type='l',col=3,lwd=2)
lines(rownames(m1),m1[,"Y_tot"],type='l',col=1,lwd=2)
legend(5.5,15.5,legend = colnames(m1)[-(4:5)], bg = 'gray90',col=c(2,3,1),pch=1,xpd=NA,cex=legend.cex)

m2<-m$ET_Main[,c("F_loss","Fish_mort")]
plot(rownames(m2),m2[,"F_loss"],col=4,pch=2,type='l',ylim=c(0.0000004,max(m2)),xlim=c(2,6),xlab="TL",ylab=NA,lwd=2,...)
lines(rownames(m2),m2[,"Fish_mort"],type='l',col=5,pch=3,lwd=2)
legend(6.7,0.4,legend = colnames(m2), bg = 'gray90',col=c(4,5),pch=1,xpd=NA,cex=legend.cex)
par(opar)

# Kinetic
#m=create.ETmain(ecopath_guinee)
opar=par(no.readonly=TRUE)
par(mfrow=c(1,1),mar=c(5, 4, 4, 2)+.1)
plot(row.names(m$ET_Main),m$ET_Main[,'Kin'],type='l',xlim=c(2,5),xlab="TL",ylab=NA,lwd=2,main="Kinetic",col='red',...)
par(opar)

# Selectivity
opar=par(no.readonly=TRUE)
plot(row.names(m$ET_Main),m$ET_Main[,'Selec'],type='l',xlim=c(2,5),xlab="TL",ylab=NA,lwd=2,main="Selectivity",col='red',...)
options(warn=-1)
par(opar)
}

