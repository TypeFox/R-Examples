graphResCov.fun <-
function(Xvar, nint,mlePP, h=NULL, typeRes='Pearson',
namX=NULL, histWgraph=TRUE, plotDisp=c(2,2), tit='')
{
Xvar<-as.matrix(Xvar)
n<-dim(Xvar)[1]
if (is.null(tit)) tit<-mlePP@tit


mXres<-NULL
mXm<-NULL
mXpc<-NULL

nXv<-dim(Xvar)[2]

if (histWgraph==TRUE)	dev.new(record=TRUE)

par(mfrow=plotDisp)
nplot<-plotDisp[1]*plotDisp[2]

auxX<-graphResX.fun(X=Xvar[,1], nint=nint, mlePP=mlePP,
typeRes=typeRes, namX=namX[1])
mXres<-cbind(mXres, auxX$Xres)
mXm<-cbind(mXm, auxX$Xm)
mXpc<-cbind(mXpc, auxX$pc)

iXv<-2
igraf<-1
while (iXv<=nXv)
{
if ((igraf-nplot*floor(igraf/nplot))==0) #it check if the number of performed plots is multiple of 4
{
	mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=0.8)
	mtext(paste(typeRes, " residuals  ", sep=' '), outer = TRUE, line = -3,cex=0.7)
}

auxX<-graphResX.fun(X=Xvar[,iXv], nint=nint,mlePP=mlePP,
typeRes=typeRes, namX=namX[iXv])
mXres<-cbind(mXres, auxX$Xres)
mXm<-cbind(mXm, auxX$Xm)
mXpc<-cbind(mXpc, auxX$pc)

iXv<-iXv+1
igraf<-igraf+1
}

mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)
mtext(paste(typeRes, " residuals. Number of intervals:  ",nint, sep=' '), outer = TRUE, line = -3,cex=0.7)

return(list(mXres=mXres,mXm=mXm,mXpc=mXpc,nint=nint, mlePP=mlePP))

}
