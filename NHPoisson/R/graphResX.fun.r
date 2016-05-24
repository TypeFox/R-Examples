graphResX.fun <-
function(X,nint,  mlePP, typeRes='Pearson',h=NULL,namX=NULL)
{

lambdafit<-mlePP@lambdafit
      n<-length(lambdafit)
posE<-mlePP@posE
inddat<-mlePP@inddat
      if (is.null(inddat)) inddat<-rep(1,n)

X[inddat==0]<-NA
lambdafit[inddat==0]<-NA

pc<-quantile(X,cumsum(rep(1/nint, (nint-1))),na.rm=TRUE)
int<-rep(0,length(X))
men.fun<-function(i,X)
{
int<<-int+(X>=i)
}
aux<-apply(as.matrix(pc),MARGIN=1,FUN=men.fun,X)
Xm<-tapply(X, INDEX=int, FUN=mean, na.rm = TRUE)
lintV<-tapply(inddat, INDEX=int, FUN=sum, na.rm=TRUE)
indice<-rep(0,n)
indice[inddat==0]<-NA

if (typeRes=='Raw')
{
indice[posE]<-1
Xsumalfit<-tapply(lambdafit, INDEX=int, FUN=mean, na.rm = TRUE)
ic<-2*(Xsumalfit/lintV)**0.5
}
else
{
if (is.null(h)) 
{
h<-1/lambdafit**0.5
typeRes<-'Pearson'
}
if (is.null(typeRes)) stop('Please indicate argument type 
of residuals')

indice[posE]<-1*h[posE]
lambdafit<-lambdafit*h
Xsumalfit<-tapply(lambdafit, INDEX=int, FUN=mean, na.rm = TRUE)
ic<-2/lintV**0.5
}

Xemplambda<-tapply(indice, INDEX=int, FUN=mean, na.rm = TRUE)

Xres<-Xemplambda-Xsumalfit

limysup<-max(Xres, ic, na.rm=TRUE)
limyinf<-min(Xres, -ic, na.rm=TRUE)

plot(Xm, Xres, xlab = namX,
ylab = paste (typeRes,"residuals", sep=' '), pch=16, cex = 0.3,
ylim=c(limyinf, limysup))
lines(Xm,ic, col='red')
lines(Xm,-ic, col='red')

if (length(Xres)<nint)
{
Xres<-c(Xres, rep(NA, nint-length(Xres)))
Xm<-c(Xm, rep(NA, nint-length(Xm)))

}
return(list(Xres=Xres,Xm=Xm,pc=pc, typeRes=typeRes, namX=namX, lambdafit=lambdafit, posE=posE))

}
