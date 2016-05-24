`InformationMatrixARz` <-
function(zeta,lags){
    pmax<-max(lags)
    z<-numeric(pmax)
    z[lags]<-zeta
    Iar<-InformationMatrixAR(PacfToAR(z))
    if (pmax == 1)
        return(Iar)
    J<-Jacobian(z)
    Iz<-t(J)%*%Iar%*%J
    q<-length(lags)
    Izeta <-matrix(rep(0,q*q),nrow=q, ncol=q)
    for (i in 1:q)
        for (j in 1:q)
            Izeta[i,j]=Iz[lags[i],lags[j]]
    Izeta
}

