plot.Transpose <-
function(x,title=NULL,scale=NULL,legend.cex=NULL,...)    #the scale parameter can be log or not for the y axis
{                                                               
tab_Trans<-x
if(missing(title))
cat("title is missing\n")
if (is.null(legend.cex)) {legend.cex=0.8}
opar=par(no.readonly=TRUE)
par(mfrow=c(2,1))
m <- (tab_Trans[-1,])
m <- m[-c(32:51),]
tmp_vecteur=m[,1]

if(is.null(scale))
{
scale <- ''

y_max<-max(apply(m,1,sum))
y_min <- 0

par(mar=c(5,3,1,15))

plot(rownames(m),apply(m,1,sum),log=scale,col=2,type='l',axes=F,xlab="TL",main=paste(title,": by group"),ylab=NA,ylim=c(y_min,y_max),xaxs = "i",xlim=c(2,5),...)

m_cum <- apply(m,1,cumsum)

for(i in length(rownames(m_cum)):1) 
{
polygon(c(2,as.numeric(colnames(m_cum)),5),c(y_min,m_cum[i,],y_min),col=i,xlim=c(2,5),border='NA')
polygon(c(2,2,5,5),c(min(m_cum[i,]),y_min,y_min,min(m_cum[i,])),col='white',border='NA')
}
axis(4,...)
legend(5.7,y_max,legend = colnames(m), bg = 'gray90',col=seq(1:ncol(m)),pch=1,xpd=NA,cex=legend.cex)
   
plot(rownames(m),apply(m,1,sum),log=scale,col=2,type='l',bg='gray',axes=F,main=paste(title,": trophic spectra"),ylab=NA,xlab="TL",ylim=c(y_min,y_max),xlim=c(2,5),...)
axis(1,...)
axis(4,...)
par(opar)

}

else
{
scale <- 'y'   


y_max<-max(apply(m,1,sum))
y_min <- sum(apply(m, 1, sum))*0.0001

par(mar=c(5,3,1,15))

#plot(rownames(m),apply(m,1,sum),log=scale,col=2,type='l',axes=F,ylab=title,xlab="TL",main=title,ylim=c(y_min,y_max),xaxs = "i",xlim=c(2,5))
plot(rownames(m),apply(m,1,sum),log=scale,col=2,type='l',axes=F,xlab="TL",main=paste(title,": by group"),ylab=NA,ylim=c(y_min,y_max),xaxs = "i",xlim=c(2,5),...)

m_cum <- apply(m,1,cumsum)

for(i in length(rownames(m_cum)):1) 
{
polygon(c(2,as.numeric(colnames(m_cum)),5),c(y_min,m_cum[i,],y_min),col=i,xlim=c(2,5),border='NA')
polygon(c(2,2,5,5),c(min(m_cum[i,]),y_min,y_min,min(m_cum[i,])),col='white',border='NA')
}
axis(4,...)
legend(5.7,y_max,legend = colnames(m), bg = 'gray90',col=seq(1:ncol(m)),pch=1,xpd=NA,cex=legend.cex)
 
plot(rownames(m),apply(m,1,sum),log=scale,col=2,type='l',bg='gray',axes=F,main=paste(title,": trophic spectra"),ylab=NA,xlab="TL",ylim=c(y_min,y_max),xlim=c(2,5),...)
axis(1,...)
axis(4,...)
par(opar)  
}
}
