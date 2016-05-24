`GetB` <-
function(phi){
 p<-length(phi)
 if (p == 1)
    a<-phi^2
 else {
    a<-p*as.vector(acf(phi,lag.max=p,type="covariance",demean=FALSE,plot=FALSE)$acf)
    for (i in 1:(p-1))
        a<-c(a,p*as.vector(acf(c(phi[-(1:i)],rep(0,i)),lag.max=p,type="covariance",demean=FALSE,plot=FALSE)$acf)[1:(p-i)])
    }
 FromSymmetricStorageUpper(a)
}

