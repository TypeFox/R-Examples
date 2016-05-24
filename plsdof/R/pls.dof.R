

pls.dof=function(pls.object,n,y,K,m,DoF.max){
    TT<-pls.object$TT
    Yhat<-pls.object$Yhat[,2:(m+1)]
    TK=matrix(,m,m)
    KY<-krylov(K,K%*%y,m)
    #KK<-diag(n)
    lambda<-eigen(K)$values
    tr.K<-vector(length=m)
    for (i in 1:m){
        #KK<-K%*%KK
        tr.K[i]<-sum(lambda^i)
        #tr.K[i]<-tr(KK)
    }
    BB=t(TT)%*%KY
    BB[row(BB)>col(BB)]=0
    b<-t(TT)%*%y
    DoF=vector(length=m)
    Binv<-backsolve(BB,diag(m))
    tkt<-rep(0,m)
    ykv<-rep(0,m)
    KjT<-array(dim=c(m,n,m))
    dummy<-TT
    for (i in 1:m){
        dummy<-K%*%dummy
        KjT[i,,]<-dummy
    }
    trace.term=rep(0,m)
    for (i in 1:m){
            Binvi<-Binv[1:i,1:i,drop=FALSE]
            ci<-Binvi%*%b[1:i]
            Vi<-TT[,1:i,drop=FALSE]%*%t(Binvi)
            trace.term[i]<-sum(ci*tr.K[1:i])
            ri<-y-Yhat[,i]
            for (j in 1:i){
                KjTj=KjT[j,,]
                tkt[i]<-tkt[i]+ci[j]*tr(t(TT[,1:i,drop=FALSE])%*%KjTj[,1:i,drop=FALSE])
                ri<-K%*%ri
                ykv[i]<-ykv[i]+ sum(ri*Vi[,j])
            }
            }
    DoF<-trace.term + 1:m - tkt + ykv
    DoF[DoF>DoF.max]=DoF.max
    sigmahat=sqrt(pls.object$RSS[-1]/(n-DoF))
        return(list(DoF=DoF,sigmahat=sigmahat))





}