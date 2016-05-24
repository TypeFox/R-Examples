`getsurv` <-
function(times,icfit,nonUMLE.method="interpolation"){
   
    ntimes<-length(times)

    method<-match.arg(nonUMLE.method,c("interpolation","left","right"))

    ## function to get survival given L,R, and p
    ## does not use LRin attributes, so values exactly on intmap values may represent 
    ## the limit approaching the intmap value
    getsurvOneStratum<-function(L,R,p){
        Sout<-rep(NA,ntimes)
        mle<-rep(TRUE,ntimes)
        S<-c(1-cumsum(p))
        nonUMLE.func<-switch(method,
            interpolation=function(i,Time){
                if (R[i]==Inf) stop("cannot interpolate when R[i]=Inf")
                ## fix May 9, 2013: had S[i-1] did not work when i==1
                ## replace with 1
                if (i==1){
                    1 + ((Time - L[i])/(R[i]-L[i]))*(S[i]-1)
                } else {
                    S[i-1] + ((Time - L[i])/(R[i]-L[i]))*(S[i]-S[i-1])}},
            left=function(i,Time){ ifelse(i<=1,1,S[i-1]) },
            right=function(i,Time){ S[i]})

        k<-length(p)
        for (i in 1:ntimes){
            if (any(times[i]==R)){ Sout[i]<-S[times[i]==R]
            } else if (times[i]<=L[1]){ Sout[i]<-1
            } else if (times[i]>=R[k]){ Sout[i]<-0
            } else {
                if  (times[i]>L[k]){ 
                    Sout[i]<-nonUMLE.func(k,times[i])
                    mle[i]<-FALSE
                } else {
                    ## iLup is the index of the smallest L endpoint
                    ## larger than times[i]
                    iLup<-min((1:k)[L>=times[i]])
                    if (R[iLup-1]<=times[i] | L[iLup]==times[i]) Sout[i]<-S[iLup-1]
                    else {
                        Sout[i]<- nonUMLE.func(iLup-1,times[i])
                        mle[i]<-FALSE
                    }
                }
           }
        }
        out<-list(S=Sout,times=times,unique.mle=mle,nonUMLE.method=method)
        out
    }
    
    if (is.null(icfit$strata)){
        #stop("icfit should have strata element")
        # instead of giving an error assume only one strata
        icfit$strata<-c(length(icfit$pf))
    }
    nstrata<-length(icfit$strata)
    strata<-icfit$strata
    cnt<-1
    for (i in 1:nstrata){
        ## fix bug, earlier had: I<-cnt:strata[i]
        I<-cnt:(cnt+strata[i]-1)
        p<-icfit$pf[I]
        L<-icfit$intmap[1,I]
        R<-icfit$intmap[2,I]
        cnt<-cnt+strata[i]
        # if there is more than one strata, make a list
        if (i==1) OUT<-list(getsurvOneStratum(L,R,p))
        else OUT<-c(OUT,list(getsurvOneStratum(L,R,p)))
    }
    if (nstrata>1){
        strataNames<-names(strata)
        names(strataNames)<-1:nstrata
        OUT<-c(OUT,list(strataNames=strataNames))
    }
    OUT
}
