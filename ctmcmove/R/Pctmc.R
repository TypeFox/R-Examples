Pctmc <- function(Q,t){
    ##
    ## Pctmc
    ## Author: Ephraim Hanks
    ## Last Update: 20160120
    ##
    ## Function to compute the probability transition matrix
    ##  of a CTMC for a given time lag
    ##
    ## Q = either the Rate matrix or the Infinitessimal Generator of the CTMC
    ## t = the time lag desired
    ##
    if(diag(Q)[1]==0){
        N=dim(Q)[1]
        R=Q
        rvals=apply(R,1,sum)
    }
    if(diag(Q)[1]>0){
        N=dim(Q)[1]
        rvals=diag(Q)
        R=-Q
        diag(R)=0
    }
    if(diag(Q)[1]<0){
        N=dim(Q)[1]
        rvals=-diag(Q)
        diag(R)=0
    }
    r=max(rvals)
    Phat=R/r
    diag(Phat) <- 1-rvals/r
    rt=r*t
    K=ceiling(max(c(rt+5*sqrt(rt),20)))
    P=diag(N)*dpois(0,rt)
    Pk=diag(N)
    for(k in 1:K){
        Pk=Pk%*%Phat
        P=P+dpois(k,rt)*Pk
    }
    P
}
