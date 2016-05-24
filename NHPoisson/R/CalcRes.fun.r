CalcRes.fun <-
function(mlePP, lint, h=NULL, typeRes=NULL)
{

n<-length(mlePP@lambdafit)
t<-mlePP@t
tit<-mlePP@tit
if (is.null(h)) 
{
h<-1/mlePP@lambdafit**0.5
typeRes<-'Pearson'
}
if (is.null(typeRes)) stop('Please indicate argument typeRes')

inddat<-mlePP@inddat
inddat[inddat==0]<-NA
posE<-mlePP@posE
lambdafit<-mlePP@lambdafit*h*inddat

lambdafitR<-mlePP@lambdafit*inddat
indice<-rep(0,n)
indice[posE]<-1*h[posE]
indice<-indice*inddat
indiceR<-rep(0,n)
indiceR[posE]<-1
indiceR<-indiceR*inddat

iini<-1
ifin<-lint
posmed<-floor(lint/2)+1

emplambda<-NULL
emplambda[1:(posmed-1)]<-NA
emplambda[posmed]<-mean(indice[iini:ifin], na.rm=T)
sumalfit <- NULL
sumalfit[1:(posmed - 1)] <- NA
sumalfit[posmed] <- mean(lambdafit[iini:ifin], na.rm=T)
emplambdaR<-NULL
emplambdaR[1:(posmed-1)]<-NA
emplambdaR[posmed]<-mean(indiceR[iini:ifin], na.rm=T)
sumalfitR <- NULL
sumalfitR[1:(posmed - 1)] <- NA
sumalfitR[posmed] <- mean(lambdafitR[iini:ifin], na.rm=T)
lintV<- NULL
lintV[1:(posmed - 1)] <- NA
lintV[posmed] <- sum(inddat[iini:ifin],na.rm=T)

j <- posmed+1
while((ifin < n))
{
iini<-iini+1
ifin<-ifin+1
emplambda[j]<-mean(indice[iini:ifin], na.rm=T)
sumalfit[j] <- mean(lambdafit[iini:ifin], na.rm=T)
emplambdaR[j]<-mean(indiceR[iini:ifin], na.rm=T)
sumalfitR[j] <- mean(lambdafitR[iini:ifin], na.rm=T)
lintV[j]<-sum(inddat[iini:ifin], na.rm=TRUE)
j<-j+1
}

sumalfit[j:n] <- NA
emplambda[j:n]<-NA
sumalfitR[j:n] <- NA
emplambdaR[j:n]<-NA
lintV[j:n]<-NA

ScaRes<-emplambda-sumalfit
RawRes<-emplambdaR-sumalfitR
return(list(RawRes=RawRes,ScaRes=list(ScaRes=ScaRes,typeRes=typeRes),
emplambda=emplambdaR,fittedlambda=sumalfitR,lintV=lintV,lint=lint,
typeI='Overlapping',h=h,mlePP=mlePP))

}
