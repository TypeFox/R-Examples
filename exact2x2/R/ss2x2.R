

ss2x2 <-function(p0,p1,power=.80,n1.over.n0=1,sig.level=0.05,alternative=c("two.sided","one.sided"),
paired=FALSE,strict=FALSE,tsmethod=NULL,nullOddsRatio=1,errbound=10^-6,print.steps=FALSE, approx=FALSE){
    # use same calls as with power2x2
    # so we can use named calls, change variable names
    SIG.LEVEL<-sig.level
    ALT<-alternative
    PAIRED<-paired
    STRICT<-strict
    TSMETHOD<-tsmethod
    NULLOR<-nullOddsRatio
    EB<-errbound

    ## create root function for uniroot.integer


    rootfunc<-function(n0){
        n0<-ceiling(n0)
        n1<-ceiling(n0*n1.over.n0)
        calPOWER<-power2x2(p0,p1,n0,n1,sig.level=SIG.LEVEL,
            alternative=ALT,paired=PAIRED,strict=STRICT,
            tsmethod=TSMETHOD,nullOddsRatio=NULLOR,errbound=EB)$power
        if (print.steps) print( paste("n0=",n0," n1=",n1," power=",calPOWER,sep=""))
        power-calPOWER
    }
    alternative<-match.arg(alternative)
    tside <- switch(alternative, one.sided = 1, two.sided = 2)
    ## to get starting value
    ## calculate sample size based on normal approximation
    ## See Fleiss, 1981 p. 45 eqn 3.18
    Ca<- qnorm(sig.level/tside)
    Cb<- qnorm(power)
    r<-n1.over.n0
    Pbar<- (p0 + r*p1)/(r+1)
    Qbar<-1-Pbar
    mprime<- (Ca*sqrt((r+1)*Pbar*Qbar) - Cb*sqrt(r*p0*(1-p0)+p1*(1-p1)) )^2/(r*(p0-p1)^2 )
    m<- (mprime/4)*(1+sqrt(1 + (2*(r+1))/(mprime*r*abs(p0-p1)) ) )^2
    N0.approx<-m
    ## simple normal approximation, no continuity correction
    #Za<- qnorm(1-sig.level/tside)
    #Zb<- qnorm(power)
    #V0<-p0*(1-p0)
    #V1<-p1*(1-p1)
    #delta<-p0-p1
    #N0.approx<- ((V0 + V1/n1.over.n0)*(Za+Zb)^2) / delta^2
    # start low, need to underestimate
    N0.start<-ceiling(0.9*N0.approx)-1

    if (approx){
        ## create pretty output using power.htest class
        if (is.null(TSMETHOD)) TSMETHOD<-"NULL (see exact2x2 help)"
        if (paired==FALSE){
            METHOD<-"Approximate Power for Fisher's Exact Test"
        } else if (paired==TRUE){
            METHOD<-"Approximate Power for McNemar's Exact Test"
        }        
        out<-list(power=power,p0=p0,p1=p1,
              n0=N0.approx,n1=n1.over.n0*N0.approx,
              sig.level=sig.level,alternative=alternative,strict=strict,tsmethod=tsmethod,
              nullOddsRatio=nullOddsRatio, note=paste("errbound=",errbound), 
              method = METHOD)
        class(out)<-"power.htest"
    } else {
        ## 2015/02/26 fix problem if starting value too large
        ## check starting value to make sure it is less 
        ## than the target power
        power.approx<-power2x2(p0,p1,N0.start,
            ceiling(N0.start*n1.over.n0),sig.level=SIG.LEVEL,
            alternative=ALT,paired=PAIRED,strict=STRICT,
            tsmethod=TSMETHOD,nullOddsRatio=NULLOR,
            errbound=EB)$power
        if (power.approx>power){
            while (power.approx>power){
                N0.start<-ceiling(N0.start/2)
                power.approx<-power2x2(p0,p1,N0.start,
                    ceiling(N0.start*n1.over.n0),
                    sig.level=SIG.LEVEL,
                    alternative=ALT,paired=PAIRED,
                    strict=STRICT,
                    tsmethod=TSMETHOD,nullOddsRatio=NULLOR,
                    errbound=EB)$power
            }
        }

        if (print.steps) print(paste("starting calculation at n0=",N0.start, 
               " n1=",ceiling(N0.start*n1.over.n0)))
        if (strict & alternative=="two.sided") warning("Exact sample size based on monotonicity of power function. This is not guaranteed when strict=TRUE and alternative='two.sided' ")

        N0<-uniroot.integer(rootfunc,lower=N0.start,upper=Inf, step.power=3,print.steps=FALSE)$root
        N1<-ceiling(n1.over.n0*N0)
        out<-power2x2(p0,p1,N0,N1,sig.level=SIG.LEVEL,
                alternative=ALT,paired=PAIRED,strict=STRICT,
                tsmethod=TSMETHOD,nullOddsRatio=NULLOR,errbound=EB)
    }
    out    
}

#ss2x2(.5,.99,power=.8,print.steps=TRUE,approx=FALSE)