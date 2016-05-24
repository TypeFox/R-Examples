`Qh.times.dhyper` <-
function(h,n1,n0,M,SM,T0,use.ranks=TRUE){
    K<-n1+n0-M
    ## for these sets of permutations, here are numbers of zeros and 
    ## ones in each group
    k1p<-h
    k0p<-K-h
    m1p<-n1-k1p
    m0p<-n0-k0p
    ## now do the chopping
    if (m0p/n0 >= m1p/n1){
        k0c<- 0
        k1c<- k1p-floor(n1*k0p/n0)
    }
    else {
        k0c<- k0p-floor(n0*k1p/n1)
        k1c<- 0
    }
    Kstar<-k0c+k1c
    m0c<- m0p
    m1star<- m1p
    hstar<- k1c
    m1star<-n1-h
    n1star<-m1star+hstar
    n0star<-m0c+k0c
    if (use.ranks) S0Kstar<- -(Kstar-1)/2
    else S0Kstar<-0
    VM<-(M-1)*var(SM)*var(c(rep(1,n1star),rep(0,n0star)))
    Nstar<-M+Kstar
    SbarNstar<- (Kstar*S0Kstar + sum(SM))/Nstar
    ## Note if ZNstar has n1star ones and n0star=Nstar-n1star zeros then 
    ## Var(ZNstar) = (n0star*n1star)/(Nstar*(Nstar-1))
    VNstar<- (Nstar-1)*var(c(rep(S0Kstar,Kstar),SM))* (Nstar-n1star)*n1star/(Nstar*(Nstar-1))
    Zstat<- (T0*sqrt(VNstar) - hstar*S0Kstar+n1star*SbarNstar - (n1star-hstar)*mean(SM) )/sqrt(VM)
    Qhhat<- c(pnorm( Zstat ))
    return(dhyper(h,K,M,n1)*Qhhat)
}

