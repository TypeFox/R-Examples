`VarianceRacfARz` <-
function(zeta, lags, MaxLag, n){
PMAX<-max(lags)
p<-max(lags)
zta<-numeric(p)
zta[lags]<-zeta
phi<-PacfToAR(zta)
psi<-ARToMA(phi, MaxLag-1)
PMAX<-length(phi)
X<-matrix(rep(c(-psi,0),PMAX)[1:(PMAX*MaxLag)],ncol=PMAX)
X<-X*outer((1:MaxLag),(1:PMAX),">=")
J<-Jacobian(zta)[,lags]
X<-X%*%J
id<-matrix(rep(c(1,rep(0,MaxLag)),MaxLag)[1:(MaxLag^2)],nrow=MaxLag,ncol=MaxLag)
(id-X%*%solve(InformationMatrixARz(zeta,lags))%*%t(X))/n
}

