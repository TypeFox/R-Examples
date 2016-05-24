`VarianceRacfARp` <-
function(phi,lags, MaxLag, n){
psi<-ARToMA(phi, MaxLag-1)
PMAX<-max(lags)
X<-matrix(rep(c(psi,0),PMAX)[1:(PMAX*MaxLag)],ncol=PMAX)
X<-X*outer((1:MaxLag),(1:PMAX),">=")
X<-X[,lags]
id<-matrix(rep(c(1,rep(0,MaxLag)),MaxLag)[1:(MaxLag^2)],nrow=MaxLag,ncol=MaxLag)
(id-X%*%solve(InformationMatrixARp(phi, lags))%*%t(X))/n
}

