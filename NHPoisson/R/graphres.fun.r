graphres.fun <-
function(objres=NULL, typeRes='Raw', t=NULL, res=NULL, 
lint=NULL, posE=NULL, fittedlambda=NULL, typeI='Disjoint',
Xvariables=NULL, namXv=NULL,  histWgraph=TRUE, plotDisp=c(2,2),
addlow=FALSE,lwd=2,tit=NULL, flow=0.5, xlegend='topleft',legcex=0.5)
{

if (is.null(objres)&(is.null(res)|is.null(t)))
stop ('Argument objres or pair of arguments (t,res) must be specified')

if (is.null(objres)=='FALSE')
{
if (is.null(typeRes)) typeRes<-'Pearson'
if (typeRes=='Raw') {res<-objres$RawRes}
else {res<-objres$ScaRes$ScaRes
typeRes<-objres$ScaRes$typeRes}
typeI<-objres$typeI
            t<-objres$mlePP@t
if (typeI=='Disjoint') 
{
t<-t[objres$pm]
Xvariables<-Xvariables[objres$pm,]
}
lint<-objres$lintV
fittedlambda<-objres$fittedlambda
posE<-objres$mlePP@posE

}
if (is.null(Xvariables)) nXv<-0
else
{
Xvariables<-as.matrix(Xvariables)
            nXv<-dim(Xvariables)[2]
}

if (length(lint)==1) lint<-rep(lint,length(t))

if (histWgraph==TRUE)	dev.new(record=TRUE)

par(mfrow=plotDisp)
nplot<-plotDisp[1]*plotDisp[2]


ic<-NULL
if (typeRes=='Pearson')
{
ic<-2/lint**0.5
}
if ((typeRes=='Raw')&(is.null(fittedlambda)==FALSE))
{
ic<-2*(fittedlambda/lint)**0.5
}

if (is.null(ic))
{
limysup<-max(res, na.rm=TRUE)
limyinf<-min(res, na.rm=TRUE)
plot(t, res, cex=1, xlab = "time", ylim=c(limyinf, limysup),
ylab = paste (typeRes,"residuals", sep=' '), type='n')
}
else
{
limysup<-max(res, ic, na.rm=TRUE)
limyinf<-min(res, -ic, na.rm=TRUE)
plot(t, res, cex=1, xlab = "time", ylim=c(limyinf, limysup),
ylab = paste (typeRes,"residuals", sep=' '), type='n')
lines(t,ic, col='red')
lines(t,-ic, col='red')
}

if (typeI=='Overlapping') 
{
points(t, res, cex = 0.3,pch=16,col='grey')
points(t[posE],res[posE] , cex = 0.3,pch=16)
legend(x=xlegend, legend=c('Residuals in the occurrence times'), col=c('black'), pch=16, cex=legcex,bty='n')
}
else
points(t, res, cex = 0.3,pch=16)

if (addlow==TRUE)
{
indna<-(is.na(res)==FALSE)
aux<-lowess(t[indna],res[indna],f=flow)
lines(aux$x,aux$y, lwd=lwd)
}

iXv<-1
igraf<-1
while (iXv<=nXv)
{
if ((igraf-nplot*floor(igraf/nplot))==0) #it check if the number of performed plots is multiple of 4
{
mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=0.8)
mtext(paste(typeI, typeRes, "residuals  ", sep=' '), outer = TRUE, line = -3,cex=0.7)
}
plot(Xvariables[,iXv], res, cex = 1, xlab = namXv[iXv], 
ylab = paste (typeRes,"residuals", sep=' '), type = "n")
if (typeI=='Overlapping') 
{
points(Xvariables[,iXv], res, cex = 0.3,pch=16,col='grey')
points(Xvariables[posE,iXv],res[posE] , cex = 0.3,pch=16)
legend(x=xlegend, legend=c('Residuals in the occurrence times'), col=c('black'), pch=16, cex=legcex,bty='n')
}
else  points(Xvariables[,iXv], res, cex = 0.3,pch=16)

if (addlow==TRUE)
{
indna<-((is.na(res)==FALSE)&(is.na(Xvariables[,iXv])==FALSE))
aux<-lowess(Xvariables[indna,iXv],res[indna],f=flow)
lines(aux$x,aux$y, lwd=lwd)
}
igraf<-igraf+1
iXv<-iXv+1
}
mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)
mtext(paste(typeI, typeRes, "residuals  ", sep=' '), outer = TRUE, line = -3,cex=0.7)

return(NULL)
}
