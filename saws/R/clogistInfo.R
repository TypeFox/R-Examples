`clogistInfo` <-
function(n,m,x,beta,h){
    ## estimate score vector and information
    ## matrix for a single set
    nbeta<-length(beta)
    ## LOGLIKP[i] gives loglik 
    ##      with beta[i]<-beta[i]+h
    ## LOGLIKM[i] gives loglik
    ##      with beta[i]<-beta[i]-h
    x<-as.matrix(x)
    LOGLIKP<-rep(0,nbeta)
    LOGLIKM<-rep(0,nbeta)
    LOGLIK<-clogistLoglike(n,m,x,beta)
    all.not.equal<-function(x) var(x)>0
    ## if all the jth column of covariates within a set 
    ## are equal (i.e. are fixed) then score and info 
    ## equal zero  for the
    ## corresponding row and/or column, so no need to 
    ## calculate them
    notfixed<-(1:nbeta)[apply(x,2,all.not.equal)]
    rinfo<-matrix(0,nbeta,nbeta)
    for (i in notfixed){
        betatemp<-beta
        betatemp[i]<-beta[i]+h
        LOGLIKP[i]<-clogistLoglike(n,m,x,betatemp)
        betatemp[i]<-beta[i]-h
        LOGLIKM[i]<-clogistLoglike(n,m,x,betatemp)
        rinfo[i,i]<- LOGLIKP[i] + LOGLIKM[i]-
            2*LOGLIK
    }
    if (length(n[notfixed])==1)
    return(list(score=(LOGLIKP-LOGLIKM)/(2*h),
             info=rinfo/(h**2),
             loglik=LOGLIK ) )
    nfc<- 0
    for (i in notfixed[-1]){
        nfc<-nfc+1
        for (j in notfixed[1:nfc]){
            betatemp<-beta
            betatemp[i]<-beta[i]+h
            betatemp[j]<-beta[j]-h
            rinfo[i,j]<-(LOGLIKP[i] + LOGLIKM[j]-
                clogistLoglike(n,m,x,betatemp) - LOGLIK)
            rinfo[j,i]<-rinfo[i,j]
        }
    }
    list(score=(LOGLIKP-LOGLIKM)/(2*h),info=rinfo/(h**2),loglik=LOGLIK)
}

