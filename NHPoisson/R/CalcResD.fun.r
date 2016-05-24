CalcResD.fun <-
function(mlePP, h=NULL, nint=NULL,lint=NULL,  
typeRes=NULL,modSim='FALSE')
{
n<-length(mlePP@lambdafit)
inddat<-mlePP@inddat
inddat[inddat==0]<-NA
posE<-mlePP@posE
if (is.null(h)) 
{
h<-1/mlePP@lambdafit**0.5
typeRes<-'Pearson'
}
if (is.null(typeRes)) stop('Please indicate argument typeRes')

if (!is.null(nint)& !is.null(lint))
stop('Error: only one of nint and lint must be supplied')
if (is.null(nint)& is.null(lint))
stop('Error: one of nint and lint must be supplied')
if (is.null(lint)) lint<-ceiling(n/nint)
if (is.null(nint)) nint<-ceiling(n/(lint))

indice<-rep(0,n)
indice[posE]<-1*h[posE]
indice<-indice*inddat
indiceR<-rep(0,n)
indiceR[posE]<-1
indiceR<-indiceR*inddat

lambdafit<-mlePP@lambdafit*h*inddat
lambdafitR<-mlePP@lambdafit*inddat

int<-floor(c(0:(n-1))/lint)
emplambda<-tapply(indice, INDEX=int, FUN=mean, na.rm = TRUE)
emplambdaR<-tapply(indiceR, INDEX=int, FUN=mean, na.rm = TRUE)
sumalfit<-tapply(lambdafit, INDEX=int, FUN=mean, na.rm = TRUE)
sumalfitR<-tapply(lambdafitR, INDEX=int, FUN=mean, na.rm = TRUE)
lintV<-tapply(inddat, INDEX=int, FUN=sum, na.rm=TRUE)

ultlint<-n-(nint-1)*lint
pm1<-floor(lint/2)
pm<-pm1+c(0,cumsum(rep(lint,(nint-2))) )
pm<-c(pm, pm[length(pm)]+ceiling(ultlint/2))

if (modSim==FALSE)
{
cat(fill=T)
cat('Number of intervals to calculate the disjoint residuals: ', nint, ' of length: ',lint, fill=T)
if(lint*nint!=n) cat(' except the las one of length ',ultlint, fill=T)
cat(fill=T)
}
ScaRes<-emplambda-sumalfit
RawRes<-emplambdaR-sumalfitR
return(list(RawRes=RawRes,ScaRes=list(ScaRes=ScaRes,typeRes=typeRes),
emplambda=emplambdaR,fittedlambda=sumalfitR, 
lintV=lintV,nint=nint, lint=lint,pm=pm,typeI='Disjoint',h=h,mlePP=mlePP))

}
