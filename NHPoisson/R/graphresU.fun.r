graphresU.fun <-
function(unires, posE, Xvariables=NULL, namXv=NULL, flow=0.5, tit='',
addlow=TRUE,histWgraph=TRUE, plotDisp=c(2,2),indgraph=FALSE)
{
tt<-c(1:length(unires))
if (is.null(Xvariables)==TRUE) nXv<-0
else nXv<-dim(Xvariables)[2]
Xvariablest<-Xvariables[posE,]
n<-length(unires)

if (histWgraph==TRUE)	dev.new(record=TRUE)


par(mfrow=plotDisp)
nplot<-plotDisp[1]*plotDisp[2]

plot(tt, unires, cex = 1, xlab = "index", 
ylab = "uniform residuals", type = "n")
points(tt, unires, cex = 0.3,pch=16)
if (addlow==TRUE)
{aux<-lowess(tt,unires,f=flow)
lines(aux$x,aux$y)
}

iXv<-1
igraf<-1

while (iXv<=nXv)
{
if ((igraf-nplot*floor(igraf/nplot))==0) #it check if the number of performed plots is multiple of nplot
{
mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)
}
plot(Xvariablest[,iXv], unires, cex = 1, xlab = namXv[iXv], 
ylab = "uniform residuals", type = "n")
points(Xvariablest[,iXv], unires, cex = 0.3,pch=16)
if (addlow==TRUE)
{aux<-lowess(Xvariablest[,iXv],unires,f=flow)
lines(aux$x,aux$y)
}
igraf<-igraf+1
iXv<-iXv+1
}

mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)


if (indgraph==TRUE) par(mfrow=c(1,1)) else par(mfrow=c(2,2))

res<-unires
lagres<-c(NA, res[1:(length(res)-1)])
plot(lagres,res,ylab='residuals (t)', xlab='residuals(t-1)', cex=0.5)
obreg<-lm(res~lagres, na.action=na.exclude)
lines(lagres,predict(obreg))
corr<-cor.test(res,lagres,na.action='remove')
title(sub=paste(' Pearson coef.: ',round(corr$estimate,2),'. P-value: ',round(corr$p.value,3),sep=''), cex=0.7)
mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)

auxacf<-acf(res, plot=T, xlab='lag',main='')
if (indgraph==TRUE) mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)


ret<-10*log10(n) 
a<-NULL
indret<-c(1:ret)
for (i in indret)
{
a[i]<-Box.test(res,type='Ljung-Box',lag=i)$p.value
}
plot(indret, a, pch=16, xlab='lag', ylab='p-values', cex=0.5, ylim=c(0,1))
abline(h=0.05,col='red')
title(sub=paste('Ljung-Box p-values ',sep=''), cex=0.7)
if (indgraph==TRUE) mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)


qqPlot(unires,dist='unif')
ksres<-ks.test(unires, y = "punif", min = 0, max = 1)
title(sub=paste('KS (uniform) p-value: ',round(ksres$p.value,3),sep=''),cex=0.7)
if (indgraph==TRUE) mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)


return(NULL)
}
