wspoissonTest<-function(x,w, nullValue=NULL, 
    alternative = c("two.sided", "less", "greater"), 
    conf.level = 0.95, 
    midp=FALSE, nmc=10^5,
    wmtype=c("max","mean","minmaxavg","tcz"),
    mult=1,useSeed=TRUE,seed=104011){
    if (midp & useSeed) set.seed(seed)
    alternative<-match.arg(alternative)
    # for one-sided intervals use either (1-conf.level) or conf.level quantile
    # for two-sided use (1-conf.level)/2 and   1- (1-conf.level)/2
    # Use alpha and 1-alpha for both types, so use  
    alpha<-1-conf.level
    if (alternative=="two.sided") alpha<-alpha/2 

    # set conf limits to extremes
    lower<-0
    upper<-Inf

    y<-sum(w*x)
    v<-sum(x*w^2)
    wmtype<-match.arg(wmtype)
    if (wmtype=="max"){
        ## for Gamma method of Fay and Feuer 
        wm<-max(w)
        ystar<- y+wm
        vstar<- v+wm^2
    } else if (wmtype=="mean"){
        ## 
        wm<-mean(w)
        ystar<-y+wm
        vstar<-v+wm^2
    } else if (wmtype=="minmaxavg"){
        wm<-mean(c(min(w),max(w)))
        ystar<-y+wm
        vstar<-v+wm^2
    } else if (wmtype=="tcz"){
        ## when midp=FALSE, this gives the Tiwari, Clegg, Zou, 2006 method
        ## denoted G4 in Ng, Filardo, and Zheng
        ystar<- y + mean(w)
        vstar<- v + mean(w^2)
    }

    if (midp){
        B<-rbinom(nmc,1,.5)
        if (y>0){  GL<-rgamma(nmc,y^2/v,scale=v/y)
        } else GL<-rep(0,nmc)
        #GU<-rgamma(nmc,(y+wm)^2/(v+wm^2),scale=(v+wm^2)/(y+wm))
        GU<-rgamma(nmc,(ystar)^2/(vstar),scale=(vstar)/(ystar))
        T<- B*GL+(1-B)*GU
        ## To do: see which type of quantile matches 
        ## p-value the best
        if (is.null(nullValue)){
            pAL<-pAG<-NA
        } else {
            pAL<- length(T[T<=nullValue])/nmc
            pAG<- length(T[T>=nullValue])/nmc
        }
        if (alternative=="two.sided"){
            ci<-quantile(T,probs=c(alpha,1-alpha))
        } else if (alternative=="less"){
            ci<-c(0,quantile(T,probs=1-alpha))
        } else if (alternative=="greater"){
            ci<-c(quantile(T,probs=alpha),Inf)
        }
    } else {
    # Gamma method of Fay and Feuer, 1997, Stat in Med, 791-801
        if (alternative=="two.sided" | alternative=="greater"){
            # remember alpha<-alpha/2 if alternative='two.sided'
            if (y>0) lower<- qgamma(alpha,y^2/v,scale=v/y)
        }
        if (alternative=="two.sided" | alternative=="less"){
            # remember alpha<-alpha/2 if alternative='two.sided'
            #upper<-qgamma(1-alpha,(y+wm)^2/(v+wm^2),scale=(v+wm^2)/(y+wm))
            upper<-qgamma(1-alpha,(ystar)^2/(vstar),scale=(vstar)/(ystar))

        }
        if (is.null(nullValue)){
            pAL<-pAG<-NA
        } else {
            # invert confidence intervals to get one-sided p-values
            pAG<- pgamma(nullValue,y^2/v,scale=v/y)
            #pAL<- 1-pgamma(nullValue,(y+wm)^2/(v+wm^2),scale=(v+wm^2)/(y+wm))
            pAL<- 1-pgamma(nullValue,(ystar)^2/(vstar),scale=(vstar)/(ystar))
        }
        ci<-c(lower,upper)
    } 
    attr(ci,"conf.level")<-conf.level
    if (alternative=="two.sided"){
        p.value<-min(1,2*pAG,2*pAL)    
    } else if (alternative=="less"){
        p.value<-pAL
    } else if (alternative=="greater"){
        p.value<-pAG
    }
    names(y) <- "Weighted Sum"
    if (mult!=1) names(y)<-paste0("Weighted Sum per ",mult)
    k<-length(x)
    names(k)<-"number of summands"
    w<-k*w/sum(w)

    parms<-c(var(w),mult)
    names(parms)<-c("standardized Variance of Weights","Multiplier")
    if (mult==1) parms<-parms[1]
    if (!is.null(nullValue)) names(nullValue) <- "Weighted Sum"
    if (midp){
        method<-"Confidence Distribution method for Weighted Sum of Poissons"
    } else {
        method<-"Gamma Method for Weighted Sum of Poissons"
    }
    if (wmtype=="max"){
        ## for Gamma method of Fay and Feuer 
        ## since original Gamma method no need to add descriptor
    } else if (wmtype=="mean"){
        method<-paste(method,"[with wm=mean(w) modification]") 
    } else if (wmtype=="minmaxavg"){
        method<-paste(method,"[with wm=mean(c(min(w),max(w))) modification]") 
    } else if (wmtype=="tcz"){
        method<-paste(method,"[with Tiwari-Clegg-Zou modification]") 
    }


    dname<- paste("x=",deparse(substitute(x)), "and w=", deparse(substitute(w)))

    attr(ci, "conf.level") <- conf.level
    rval <- list(statistic = k, parameter = parms, p.value = p.value, 
        conf.int = ci*mult, estimate = y*mult, null.value = nullValue, 
        alternative = alternative, method = method, data.name = dname)
    class(rval) <- "htest"
    return(rval)
}


#xfive<-c(0,8,63,112,262,295)
#nfive<-c(327,30666,123419,149919,104088,34392)
#ntotal<-c(319933,931318,786511,488235,237863,61313)

#set.seed(1031)
#x<-wspoissonTest(xfive,10^5*ntotal/(nfive*sum(ntotal)),nullValue=170.70168,midp=FALSE)

#wspoissonTest(xfive,10^5*ntotal/(nfive*sum(ntotal)),nullValue=170.70168,midp=FALSE)
#wspoissonTest(xfive,10^5*ntotal/(nfive*sum(ntotal)),nullValue=170.70168,midp=FALSE,wmtype="minmaxavg")
#wspoissonTest(xfive,10^5*ntotal/(nfive*sum(ntotal)),nullValue=170.70168,midp=FALSE,wmtype="tcz")
#wspoissonTest(xfive,10^5*ntotal/(nfive*sum(ntotal)),nullValue=170.70168,midp=FALSE,wmtype="mean")



