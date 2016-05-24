meldMC<-function(T1,T2, nullparm=NULL, parmtype=c("difference","oddsratio","ratio"),
    conf.level=0.95, 
    alternative=c("two.sided","less","greater"),
    dname="",estimate1=NA, estimate2=NA){

    ptype<-match.arg(parmtype)
    if (ptype=="difference"){
        g<-function(T1,T2){ T2-T1 }
        if (is.null(nullparm)) nullparm<-0
    } else if (ptype=="ratio"){
        g<-function(T1,T2){ T2/T1  }
        if (is.null(nullparm)) nullparm<-1
    } else if (ptype=="oddsratio"){
        g<-function(T1,T2){ T2*(1-T1)/( T1*(1-T2) ) }
        if (is.null(nullparm)) nullparm<-1
    }

    lowerLimit<- g(1,0)
    upperLimit<-g(0,1)

    alt<-match.arg(alternative)

    lower<-upper<-NA

    if (alt=="two.sided"){
        dolo<-dohi<-TRUE
        alpha<-(1-conf.level)/2  
    } else if (alt=="less"){
        ## alt=less so lower interval is lowest possible, 
        ## do not calculate
        dolo<-FALSE
        lower<- lowerLimit
        dohi<-TRUE
        alpha<-1-conf.level
    } else if (alt=="greater"){
        # alt=greater so upper interval is highest possible, do not calculate
        dolo<-TRUE
        dohi<-FALSE
        upper<- upperLimit
        alpha<-1-conf.level
    } else stop("alternative must be 'two.sided', 'less', or 'greater' ")

    Beta<-g(T1,T2)
    oBeta<-Beta[order(Beta)]
    NMC<-length(Beta)

    if (dolo){
        # see Efron and Tibshirani (1993) Introduction to the Bootstrap,  p. 160, bottom
        k<-floor((NMC+1)*alpha)
        if (k==0) warning("increase nmc, confidence limits may be anti-conservative")
        lower<-oBeta[k]
        pg<-length(Beta[Beta<=nullparm])/NMC
    }
    if (dohi){
        # see Efron and Tibshirani (1993) Introduction to the Bootstrap,  p. 160, bottom
        k<-floor((NMC+1)*alpha)
        if (k==0) warning("increase nmc, confidence limits may be anti-conservative")
        upper<-oBeta[NMC+1-k]
        pl<-length(Beta[Beta>=nullparm])/NMC

    }



    # we do not get estimates from betaParms
    estimate<-g(estimate1,estimate2)
    names(estimate)<-ptype

    if (alt=="two.sided"){
        p.value<- min(1,2*pl,2*pg)
    } else if (alt=="less"){
        p.value<- pl
    } else if (alt=="greater"){
        p.value<- pg
    }
    ci<-c(lower,upper)
    attr(ci,"conf.level")<-conf.level
    #dname<-paste("sample 1:(",x1,"/",n1,"), sample 2:(",x2,"/",n2,")",sep="")
    #dname<-dname
    #method<-paste("exact melded test for two binomials")
    method<-"melded test"
    stat<-estimate1
    parm<-estimate2
    names(stat) <- "estimate 1"
    names(parm) <- "estimate 2"
    names(nullparm)<-paste("Null",ptype)

    structure(list(statistic = stat, parameter = parm, 
        p.value = p.value, 
        conf.int = ci, estimate = estimate, null.value = nullparm, 
        alternative = alt, method = method, 
        data.name = dname), class = "htest")

}








betaMeldTest <-
function(betaParms1,
    betaParms2,nullparm=NULL, 
    parmtype=c("difference","oddsratio","ratio"),
    conf.level=0.95, conf.int=TRUE,
    alternative=c("two.sided","less","greater"),eps=10^-8,
    dname="",
    estimate1=NA, estimate2=NA){

    # check betaParms1 and betaParms2 to make sure they
    # are in the correct format
    checkBetaParms<-function(betaParm){
        betaParmNames<-sort(c("alower",
            "blower","aupper","bupper"))
        if (is.list(betaParm)){
            if (!all(sort(names(betaParm))==betaParmNames)){
                stop("list must have named elements, 
                   'alower', 'blower', 'aupper', and 'bupper' ")
            } else if (any(betaParm$alower<0) |
                       any(betaParm$aupper<0) |
                       any(betaParm$blower<0) |
                       any(betaParm$bupper<0) ) {
                stop("betaParms cannot have elements less than 0")
            }
        } else stop("betaParm should be a list")
    }
    checkBetaParms(betaParms1)
    checkBetaParms(betaParms2)


    aL1<-betaParms1$alower
    aU1<-betaParms1$aupper
    bL1<-betaParms1$blower
    bU1<-betaParms1$bupper


    aL2<-betaParms2$alower
    aU2<-betaParms2$aupper
    bL2<-betaParms2$blower
    bU2<-betaParms2$bupper



    ptype<-match.arg(parmtype)
    if (ptype=="difference"){
        g<-function(T1,T2){ T2-T1 }
        if (is.null(nullparm)) nullparm<-0
    } else if (ptype=="ratio"){
        g<-function(T1,T2){ T2/T1  }
        if (is.null(nullparm)) nullparm<-1
    } else if (ptype=="oddsratio"){
        g<-function(T1,T2){ T2*(1-T1)/( T1*(1-T2) ) }
        if (is.null(nullparm)) nullparm<-1
    }

    lowerLimit<- g(1,0)
    upperLimit<-g(0,1)



    # Assume 
    # X1 ~ binom(n1,p1)
    # X2 ~ binom(n2,p2)
    #
    # with g(p1,p2)= p2-p1  (difference)
    #   or g(p1,p2)= p2/p1  (ratio)
    #   or g(p1,p2)= p2(1-p1)/(p1(1-p2))  (oddsratio)
    #
    #  Want to test with D=nullparm
    #    greater:   H0: g(p1,p2) <= D
    #               H1: g(p1,p2) > D             
    # or less:      H0: g(p1,p2) >= D
    #               H1: g(p1,p2) < D   
    #
    #
    # or two.sided, which has pvalue= min(1, 2*pg, 2*pl)
    #      where pg is p-value associated with "greater" alt Hyp
    #            pl is p-value associated with "less"  alt Hyp
    #   
    #    for greater (calculate lower CL) we use   
    # T1 ~ Beta(aU1,bU1)
    # T2 ~ Beta(aL2,bL2)
    #    and p-value is calculated under null: Pr[ g(T1,T2) <= D ]
    #
    #    for less (calculate upper CL) we use 
    # T1 ~ Beta(aL1,bL1)
    # T2 ~ Beta(aU2,bU2)
    #    and p-value is calculated under null: Pr[ g(T1,T2) >= D ]

    if (ptype=="difference"){ 
        # Pr[ g(T1,T2) >=D] = Pr[ T2-T1 >= D] = Pr[ T1 <= T2 - D]
        # = \int F1(t2 - D) f2(t2)  dt2
        funcLess<-function(t2,D){
            pbeta(t2-D,aL1,bL1)*dbeta(t2,aU2,bU2)
        }
        # Pr[ g(T1,T2) <=D] = Pr[ T2-T1 <= D] = Pr[ T2 <= T1 + D]
        # = \int F2(t1 + D) f1(t1)  dt1
        funcGreater<-function(t1,D){
            pbeta(t1+D,aL2,bL2)*dbeta(t1,aU1,bU1)
        }
    } else if (ptype=="ratio"){
        # Pr[ g(T1,T2) >=D] = Pr[ T2/T1 >= D] = Pr[ T1 <= T2/D]  
        # = \int F1(t2/D) f2(t2)  dt2
        funcLess<-function(t2,D){
            pbeta(t2/D,aL1,bL1)*dbeta(t2,aU2,bU2)
        }
        # Pr[ g(T1,T2) <=D] = Pr[ T2/T1 <= D] = Pr[ T2 <= T1*D]
        # = \int F2(t1 * D) f1(t1)  dt1
        funcGreater<-function(t1,D){
            pbeta(t1*D,aL2,bL2)*dbeta(t1,aU1,bU1)
        }
    } else if (ptype=="oddsratio"){
        # Pr[ g(T1,T2) >=D] = Pr[ T2(1-T1)/T1(1-T2) >= D] 
        # = Pr[ T2(1-T1) >= T1(1-T2)D] = Pr[ T2 >= T1*(T2 + (1-T2)D) ]
        # = Pr[ T1 <= T2/{ T2 + (1-T2)D } ]  
        # = \int F1(t2/(t2+(1-t2)D) f2(t2)  dt2
        funcLess<-function(t2,D){
            pbeta(t2/(t2+(1-t2)*D),aL1,bL1)*dbeta(t2,aU2,bU2)
        }
        # Pr[ g(T1,T2) <=D] = Pr[ T2(1-T1)/T1(1-T2) <= D] 
        # = Pr[ T2(1-T1) <= T1(1-T2)D] = Pr[ T2(1-T1) + T1*T2*D <= T1*D  ]
        # = Pr[ T2( (1-T1) + T1*D) <= T1*D ]
        # = Pr[ T2 <= T1*D/{(1-T1) + T1*D} ]
        # = \int F2(t1*D/(1-t1+t1*D)) f1(t1)  dt1
        funcGreater<-function(t1,D){
            pbeta(t1*D/(1-t1+t1*D),aL2,bL2)*dbeta(t1,aU1,bU1)
        }
    }


    # p-value functions 
    pGreater<-function(delta){
        ## for the integrate function to work well, 
        ## pick values that 
        ## make sense  in funcGreater
        ## recall funcGreater is 
        ##   pbeta2( W[t] )* dbeta1(t) 
        ##       where pbeta2(.)=pbeta(.,x2,n2-x2+1)
        ##             dbeta1(.)=dbeta(.,x1+1,n1-x1)
        ## and W[t] is different depending on the parmtype
        ##    difference: W[t] = t + D 
        ##       ratio:   W[t] = t*D
        ##   odds ratio:  W[t] = t*D/(1-t+t*D)
        ##
        ## First, choose LowerInt=a and UpperInt=b so 
        ## that \int_a^b  dbeta1(t) dt = 1- eps/2 
        LowerInt<-qbeta(eps/4,aU1,bU1)
        UpperInt<-qbeta(1-eps/4,aU1,bU1)
        ##  Second, choose a2 so that pbeta2( W[a2] )= eps/2
        ##       or   qbeta2( eps/2) = W[a2] 
        q<- qbeta(eps/2, aL2,bL2)
        if (ptype=="difference"){
            a2<- q-delta
        } else if (ptype=="ratio"){
            a2<-q/delta
        } else if (ptype=="oddsratio"){
           # solve q = t*D/(1-t+t*D)   for t
           a2<- q/(delta+q-delta*q)
        }
        LowerInt<-max(a2,LowerInt)

        pout<-rep(0,length(delta))
        for (i in 1:length(delta)){
            if (LowerInt<UpperInt){
                pout[i]<-integrate(funcGreater,LowerInt,
                   UpperInt,D=delta[i])$value
            }
        }
        ## since we underestimate the integral 
        ## (assuming perfect integration) by at most eps, 
        ## add back eps to get conservative p-value
        pout<-pout+eps
        pout
    }
    pLess<-function(delta){
        ## for the integrate function to work well, 
        ## pick values that 
        ## make sense  in funcLess
        ## recall funcLess is 
        ##   pbeta1( W[t] )* dbeta2(t) 
        ##       where pbeta1(.)=pbeta(.,x1,n1-x1+1)
        ##             dbeta2(.)=dbeta(.,x2+1,n2-x2)
        ## and W[t] is different depending on the parmtype
        ##    difference: W[t] = t - D 
        ##       ratio:   W[t] = t/D
        ##   odds ratio:  W[t] = t/(t+(1-t)*D)
        ##
        ## First, choose LowerInt=a and UpperInt=b so 
        ## that \int_a^b  dbeta2(t) dt = 1- eps/2 
        LowerInt<-qbeta(eps/4,aU2,bU2)
        UpperInt<-qbeta(1-eps/4,aU2,bU2)
         
        ##  Second, choose a1 so that pbeta1( W[a1] )= eps/2
        ##       or   qbeta1( eps/2) = W[a1] 
        q<- qbeta(eps/2, aL1,bL1)
        if (ptype=="difference"){
            a2<- q+delta
        } else if (ptype=="ratio"){
            a2<-q*delta
        } else if (ptype=="oddsratio"){
           # solve q = t/(t+(1-t)*D)   for t
           a2<- q*delta/(1-q+delta*q)
        }
        LowerInt<-max(a2,LowerInt)
        pout<-rep(0,length(delta))
        for (i in 1:length(delta)){
            if (LowerInt<UpperInt){
                pout[i]<-integrate(funcLess,LowerInt,
                       UpperInt,D=delta[i])$value
            }
        }
        ## since we underestimate the integral 
        ## (assuming perfect integration) by at most eps, 
        ## add back eps to get conservative p-value
        pout<-pout+eps
        pout
    }
    alt<-match.arg(alternative)

    lower<-upper<-NA

    if (alt=="two.sided"){
        dolo<-dohi<-TRUE
        alpha<-(1-conf.level)/2  
    } else if (alt=="less"){
        ## alt=less so lower interval is lowest possible, 
        ## do not calculate
        dolo<-FALSE
        lower<- lowerLimit
        dohi<-TRUE
        alpha<-1-conf.level
    } else if (alt=="greater"){
        # alt=greater so upper interval is highest possible, 
        # do not calculate
        dolo<-TRUE
        dohi<-FALSE
        upper<- upperLimit
        alpha<-1-conf.level
    } else stop("alternative must be 'two.sided', 'less', or 'greater' ")

    


    # we do not get estimates from betaParms
    estimate<-g(estimate1,estimate2)
    names(estimate)<-ptype

    if (dolo){
        ## Take care of special cases, 
        ## when x2=0 T2 is point mass at 0
        ## when x1=n1 T1 is a point mass at 1
        if (aL2==0 & bU1>0){
            if (conf.int) lower<- g( qbeta(1-alpha,aU1,bU1), 0 )
            if (ptype=="difference"){ 
                pg<- 1- pbeta(-nullparm,aU1,bU1)
            } else pg<-1
        } else if (bU1==0 & aL2>0){
            if (conf.int) lower<- g( 1, qbeta(alpha,aL2,bL2) )
            if (ptype=="difference"){
                # Jul 15 2014: Fix error
                #pg<- 1-pbeta(1+nullparm,aL2,bL2)
                pg<- pbeta(1+nullparm,aL2,bL2)
            } else if (ptype=="ratio"){
                pg<-pbeta(nullparm,aL2,bL2)
            } else if (ptype=="oddsratio"){
                pg<-pbeta(nullparm/(1-nullparm),aL2,bL2)
            }
        } else if (aL2==0 & bU1==0){
            if (conf.int) lower<- g( 1, 0 )
            pg<-1
        } else {
            if (conf.int){ 
                rootfunc<-function(delta){
                    pGreater(delta)-alpha
                }
                if (upperLimit==Inf){
                    # uniroot cannot take Inf as an upper limit
                    # find T such that rootfunc(T) has 
                    # opposite sign as
                    # rootfunc(lowerLimit)
                    for (i in 1:50){
                        upperLimit<-2^i
                        if (sign(rootfunc(lowerLimit))!=
                            sign(rootfunc(upperLimit))){
                            break()
                        }
                    }
                }
                if (upperLimit==2^50){
                    warning("lower conf limit appears to 
                              be larger than 
                              2^50=approx=10^16, set to 2^50")
                    lower<-2^50
                } else {
                    lower<-uniroot(rootfunc,
                         c(lowerLimit,upperLimit))$root
                }
            }
            pg<-pGreater(nullparm)
            #reset upperLimit
            upperLimit<-g(0,1)
        }
    }
    if (dohi){
        ## Take care of special cases, 
        ## when bU2=0 T2 is point mass at 1
        ## when aL1=0 T1 is a point mass at 0
        if (aL1==0 & bU2==0){
            if (conf.int) upper<-g(0,1)
            pl<-1
        } else if (aL1==0){
            if (conf.int) upper<- g(0,qbeta(1-alpha,aU2,bU2))
            if (ptype=="difference"){
               pl<-1-pbeta(nullparm, aU2,bU2)
            } else if (ptype=="ratio"){
               pl<- 1
            } else if (ptype=="oddsratio"){
               pl<-1
            }
        } else if (bU2==0){
            if (conf.int) upper<- g(qbeta(alpha,aL1,bL1), 1)
            if (ptype=="difference"){
               pl<-pbeta(1+nullparm, aL1, bL1)
            } else if (ptype=="ratio"){
               pl<-pbeta(1/nullparm,aL1,bL1)
            } else if (ptype=="oddsratio"){
               pl<- 1
            }
        } else {
            if (conf.int){
                rootfunc<-function(delta){
                    pLess(delta)-alpha
                }
                if (upperLimit==Inf){
                    # uniroot cannot take Inf as an upper limit
                    # find T such that rootfunc(T) has 
                    # opposite sign as
                    # rootfunc(lowerLimit)
                    for (i in 1:50){
                        upperLimit<-2^i
                        if (sign(rootfunc(lowerLimit))!=
                            sign(rootfunc(upperLimit))){
                            break()
                        }
                    }
                }
                if (upperLimit==2^50){
                    warning("upper conf limit appears to be
                              larger than 2^50=approx=10^15, 
                              set to Inf")
                    upper<-Inf
                } else {
                    upper<-uniroot(rootfunc,
                        c(lowerLimit,upperLimit))$root
                }
            }
            pl<-pLess(nullparm)
        }
    }
    if (alt=="two.sided"){
        p.value<- min(1,2*pl,2*pg)
    } else if (alt=="less"){
        p.value<- pl
    } else if (alt=="greater"){
        p.value<- pg
    }
    ci<-c(lower,upper)
    attr(ci,"conf.level")<-conf.level
    #dname<-paste("sample 1:(",x1,"/",n1,"), 
    # sample 2:(",x2,"/",n2,")",sep="")
    #dname<-dname
    #method<-paste("exact melded test for two binomials")
    method<-"melded test using two beta CDs"
    stat<-estimate1
    parm<-estimate2
    names(stat) <- "estimate 1"
    names(parm) <- "estimate 2"
    names(nullparm)<-paste("Null",ptype)

    structure(list(statistic = stat, parameter = parm, 
        p.value = p.value, 
        conf.int = ci, estimate = estimate, 
        null.value = nullparm, 
        alternative = alt, method = method, 
        data.name = dname), class = "htest")

}



betaMeldTestMidp.mc <-
function(betaParms1,
    betaParms2,nullparm=NULL, 
    parmtype=c("difference","oddsratio","ratio"),
    conf.level=0.95, conf.int=TRUE,
    alternative=c("two.sided","less","greater"),
    dname="",
    estimate1=NA, estimate2=NA, nmc=10^6){

    # check betaParms1 and betaParms2 to make sure they
    # are in the correct format
    checkBetaParms<-function(betaParm){
        betaParmNames<-sort(c("alower",
            "blower","aupper","bupper"))
        if (is.list(betaParm)){
            if (!all(sort(names(betaParm))==betaParmNames)){
                stop("list must have named elements, 
                   'alower', 'blower', 'aupper', and 'bupper' ")
            } else if (betaParm$alower<0 |
                       betaParm$aupper<0 |
                       betaParm$blower<0 |
                       betaParm$bupper<0 ) {
                stop("betaParms cannot have elements less than 0")
            }
        } else stop("betaParm should be a list")
    }
    checkBetaParms(betaParms1)
    checkBetaParms(betaParms2)


    aL1<-betaParms1$alower
    aU1<-betaParms1$aupper
    bL1<-betaParms1$blower
    bU1<-betaParms1$bupper


    aL2<-betaParms2$alower
    aU2<-betaParms2$aupper
    bL2<-betaParms2$blower
    bU2<-betaParms2$bupper

    # sample in a balanced way, so that we always get nmc from aL1, bL1 and nmc from aU1, bU1.
    # similarly for sample 2

    # rrbeta(n,a,b) is rbeta(n,a,b) but allow limits at 0 for parameters, so that 
    #      a=0, b>0 gives rep(0,n)  
    #      a>0, b=0 gives rep(1,n)
    rrbeta<-function(n,a,b){
        if (length(a)>1 | length(b)>1) stop("rewrite rrbeta for vector parameters")
        if (a==0 & b==0){ 
            out<-rep(NA,n)
        } else if (a==0){ 
            out<-rep(0,n)
        } else if (b==0){ 
            out<-rep(1,n)
        } else {
            out<-rbeta(n,a,b)
        }
        out
    }


    T1<-sample(c(rrbeta(ceiling(nmc/2),aL1,bL1),
                 rrbeta(ceiling(nmc/2),aU1,bU1)), replace=FALSE) 
    T2<-sample(c(rrbeta(ceiling(nmc/2),aL2,bL2),
                 rrbeta(ceiling(nmc/2),aU2,bU2)), replace=FALSE) 

    mcout<-meldMC(T1,T2, 
        nullparm=nullparm,
        parmtype=parmtype,
        conf.level=conf.level,
        alternative=alternative,
        dname=dname, estimate1=estimate1, estimate2=estimate2)   
    mcout
}


# Test function
#x1<-9
#n1<-15
#x2<-13
#n2<-13
#betaMeldTestMidp.mc(betaParms1=list(alower=x1,blower=n1-x1+1,aupper=x1+1,bupper=n1-x1),
#    betaParms2=list(alower=x2,blower=n2-x2+1,aupper=x2+1,bupper=n2-x2),nullparm=NULL, 
#    parmtype=c("difference","oddsratio","ratio"),
#    conf.level=0.95, conf.int=TRUE,
#    alternative=c("two.sided","less","greater"),
#    dname="",
#    estimate1=NA, estimate2=NA, nmc=10^6)








bpcp2sampControl<-function(Delta=0, stype="km", eps=10^-8, nmc=10^6, method="mm.mc", seed=391291){
    if (!(stype=="km" | stype=="mue")) stop("stype should be 'km' or 'mue', see bpcp help")
    if (eps<0 | eps>.1) stop("eps not in reasonable range, see help for betaMeldTest")
    if (Delta<0) stop("Delta is width of grouped confidence intervals, must be greater than 0")
    is.wholenumber <-  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if (!is.wholenumber(nmc) | nmc<0) stop("nmc should be a non-negative integer (but can be numeric)")
    if (!(method=="mm.mc" | method=="mc.mc")) stop("method should be 'mm.mc' (method of moments for 1 sample, Monte Carlo for melding) or 'mc.mc' (MC for both) ") 
    list(Delta=Delta,stype=stype, eps=eps, nmc=nmc, method=method, seed=seed)
}



bpcp2samp<-function(time,status,group, testtime, 
    parmtype=c("difference","oddsratio","ratio"),
    nullparm=NULL,
    alternative=c("two.sided","less","greater"), conf.level=0.95, 
    midp=FALSE, 
    control=bpcp2sampControl()){

    # so we get the same answer on the same data set
    # for simulations use control=bpcp2sampControl(seed=NULL)
    if (!is.null(control$seed)) set.seed(control$seed)

    if (class(group)=="factor"){
        # use order of levels, but only keep the levels with data 
        ug<- levels(group)
        ug<- ug[ug %in% group]
    } else ug<-sort(unique(group))

    ## argument checking
    if (length(ug)!=2) stop("group does not have 2 levels") 
    if ((length(group[group==ug[1]])<1) | (length(group[group==ug[2]])<1)){ 
        stop("should have at least one observation in each group") }
    if (length(time)!=length(status) | length(time)!=length(group) ) stop("length of time, status and group should be the same")
    if (length(testtime)>1 | !is.numeric(testtime[1]) | testtime[1]<0) stop("testtime must be a vector of length 1")
    ## end of argument checking

    ## to avoid problems with same name of argument and 
    ## value, rename some objects
    Nullparm<-nullparm
    CL<-conf.level
    alt<-match.arg(alternative)
    
    ptype<-match.arg(parmtype)
    if (ptype=="difference"){
      Dname<-paste("S(",testtime,";group=",ug[2],
                 ")-S(",testtime,";group=",ug[1],")",sep="")
    } else if (ptype=="ratio"){
      Dname<-paste("S(",testtime,";group=",ug[2],")/
                    S(",testtime,";group=",ug[1],")",sep="")
    } else if (ptype=="oddsratio"){
      Dname<-paste("odds[S(",testtime,";group=",ug[2],
                ")]/odds[S(",testtime,";group=",ug[1],
                 ")]",sep="")
    }
    
    # regardless of method, need estimate, so run 
    # fastest (nmc=0) for single sample
    I<-group==ug[1]
    fit1<-bpcp(time[I],status[I],
              Delta=control$Delta,
              stype=control$stype, midp=midp)
    I<-group==ug[2]
    fit2<-bpcp(time[I],status[I],
              Delta=control$Delta, 
              stype=control$stype, midp=midp)
            
    getList<-function(fit,testtime){
      # get list of results at testtime
      i<- (fit$L<testtime & testtime<fit$R) |
        (fit$L<=testtime & fit$Lin & testtime<fit$R) | 
        (fit$L<testtime  & testtime<=fit$R & fit$Rin) |
        (fit$L<=testtime & fit$Lin & testtime<=fit$R & fit$Rin)
      if (length(i[i])!=1) stop("at testtime not picking one interval from fit")
      out<-list(estimate=fit$surv[i],
                betaParms=list(
                  alower=fit$betaParms$alower[i],
                  blower=fit$betaParms$blower[i],
                  aupper=fit$betaParms$aupper[i],
                  bupper=fit$betaParms$bupper[i]))
      out
    }
    list1<-getList(fit1,testtime)
    list2<-getList(fit2,testtime)

    
    if (control$method=="mm.mc"){  
      if (midp){
        testout<-betaMeldTestMidp.mc(betaParms1=list1$betaParms,
                              betaParms2=list2$betaParms,
                              nullparm=Nullparm,
                              parmtype=ptype,
                              conf.level=CL,
                              alternative=alt,
                              dname=Dname,
                              estimate1=list1$estimate,
                              estimate2=list2$estimate,
                              nmc=control$nmc)  
      } else {
        testout<-betaMeldTest(betaParms1=list1$betaParms,
                              betaParms2=list2$betaParms,
                              nullparm=Nullparm,
                              parmtype=ptype,
                              conf.level=CL,
                              alternative=alt,
                              eps=control$eps,
                              dname=Dname,
                              estimate1=list1$estimate,
                              estimate2=list2$estimate)  
      }
    } else if (control$method=="mc.mc"){
        if (midp){
            I<-group==ug[1] 
            x1<-kmgw.calc(time[I],status[I],keepCens=TRUE)
            mc1<-bpcp.mc(x1,nmc=control$nmc,
                    testtime=testtime,DELTA=control$Delta,
                    midp=TRUE)
            I<-group==ug[2]
            x2<-kmgw.calc(time[I],status[I],keepCens=TRUE)
            mc2<-bpcp.mc(x2,nmc=control$nmc,
                 testtime=testtime,DELTA=control$Delta,
                 midp=TRUE)
            testout<-meldMC(c(mc1$Smc$Slo,mc1$Smc$Shi), 
                            c(mc2$Smc$Slo,mc2$Smc$Shi),
                              nullparm=Nullparm,
                              parmtype=ptype,
                              conf.level=CL,
                              alternative=alt,
                              dname=Dname,
                              estimate1=list1$estimate,
                              estimate2=list2$estimate) 

        } else {
            # midp=FALSE
            I<-group==ug[1] 
            x1<-kmgw.calc(time[I],status[I],keepCens=TRUE)
            mc1<-bpcp.mc(x1,nmc=control$nmc,
                 testtime=testtime,DELTA=control$Delta, 
                 midp=FALSE)
            I<-group==ug[2]
            x2<-kmgw.calc(time[I],status[I],keepCens=TRUE)
            mc2<-bpcp.mc(x2,nmc=control$nmc,
                     testtime=testtime,DELTA=control$Delta,
                     midp=FALSE)
            if (alt=="two.sided"){
                dolo<-dohi<-TRUE
                alpha<-(1-conf.level)/2  
            } else if (alt=="less"){
                ## alt=less so lower interval is lowest 
                ## possible, do not calculate
                dolo<-FALSE
                dohi<-TRUE
                alpha<-1-conf.level
            } else if (alt=="greater"){
                # alt=greater so upper interval is 
                # highest possible, do not calculate
                dolo<-TRUE
                dohi<-FALSE
                alpha<-1-conf.level
            } else stop("alternative must be 'two.sided', 'less', or 'greater' ")

            if (dolo){
                # take Shi from group 1, and Slo from group 2
                testout.lo<-meldMC(c(mc1$Smc$Shi), 
                              c(mc2$Smc$Slo),
                              nullparm=Nullparm,
                              parmtype=ptype,
                              conf.level=1-alpha,
                              alternative="greater",
                              dname=Dname,
                              estimate1=list1$estimate,
                              estimate2=list2$estimate) 
            } 
            if (dohi){
                # take Slo from group 1, and Shi from group 2
                testout.hi<-meldMC(c(mc1$Smc$Slo), 
                              c(mc2$Smc$Shi),
                              nullparm=Nullparm,
                              parmtype=ptype,
                              conf.level=1-alpha,
                              alternative="less",
                              dname=Dname,
                              estimate1=list1$estimate,
                              estimate2=list2$estimate) 
            } 
            if (alt=="two.sided"){
                testout<-testout.lo
                # two.sided p-value is min(1,2*plower,2*pupper)
                testout$p.value<-min(1,2*testout.lo$p.value,
                                       2*testout.hi$p.value)
                testout$conf.int<-c(testout.lo$conf.int[1],
                                    testout.hi$conf.int[2])
                testout$alternative<-alt
            } else if (alt=="less"){
                testout<-testout.hi
            } else if (alt=="greater"){
                testout<-testout.lo
            }
        }
    } else stop("control()$method should equal 'mm.mc' or 'mc.mc' ")

    testout$method<-"Two-Sample Melded BPCP Test"
    if (midp) testout$method<- "Two-Sample Melded BPCP Test (mid-p version)"
    names(testout$statistic)<-paste("S(",testtime,";group=",ug[1],")",sep="")
    names(testout$parameter)<-paste("S(",testtime,";group=",ug[2],")",sep="")
    testout

}

#kmgw.calc(1:10,c(1,1,0,1,0,1,1,0,1,1),keepCens=TRUE)

#bpcp2samp(c(1,3,5,4,9,2,7.1,9.5,8,10),rep(1,10),c(rep(1,6),rep(2,4)),5,control=bpcp2sampControl(method="mc.mc"),midp=FALSE)
#bpcp2samp(c(1,3,5,4,9,2,7.1,9.5,8,10),rep(1,10),c(rep(1,6),rep(2,4)),5,control=bpcp2sampControl(method="mm.mc"),midp=FALSE)
#bpcp2samp(c(1,3,5,4,9,2,7.1,9.5,8,10),rep(1,10),c(rep(1,6),rep(2,4)),5,control=bpcp2sampControl(method="mc.mc"),parmtype="oddsratio",midp=FALSE)
#bpcp2samp(c(1,3,5,4,9,2,7.1,9.5,8,10),rep(1,10),c(rep(1,6),rep(2,4)),5,control=bpcp2sampControl(method="mm.mc"),parmtype="oddsratio",midp=FALSE)
#bpcp2samp(c(1,3,5,4,9,2,7.1,9.5,8,10),rep(1,10),c(rep(1,6),rep(2,4)),5,control=bpcp2sampControl(method="mc.mc"),parmtype="ratio",midp=FALSE)
#bpcp2samp(c(1,3,5,4,9,2,7.1,9.5,8,10),rep(1,10),c(rep(1,6),rep(2,4)),5,control=bpcp2sampControl(method="mm.mc"),parmtype="ratio",midp=FALSE)





fixtdiff<-function(time,status,group, testtime, 
    trans=c("identity","cloglog","log"),
    varpooled=TRUE, correct=FALSE, doall=FALSE){

 
    if (class(group)=="factor"){
        # use order of levels, but only keep the levels with data 
        ug<- levels(group)
        ug<- ug[ug %in% group]
    } else ug<-sort(unique(group))

    ## argument checking
    if (length(ug)!=2) stop("group does not have 2 levels") 
    if ((length(group[group==ug[1]])<1) | (length(group[group==ug[2]])<1)){ 
        stop("should have at least one observations in each group") }
    if (length(time)!=length(status) | length(time)!=length(group) ) stop("length of time, status and group should be the same")
    if (length(testtime)>1 | !is.numeric(testtime[1]) | testtime[1]<=0) stop("testtime must be positive and a vector of length 1")
    ## end of argument checking

    xp<-kmgw.calc(time,status)
    I<-group==ug[1]
    n1<-length(group[I])
    x1<-kmgw.calc(time[I],status[I])
    I<-group==ug[2]
    n2<-length(group[I])
    x2<-kmgw.calc(time[I],status[I])

    getSVar<-function(fit,testtime){
        # get list of results at testtime
        time<-c(0,fit$time)
        i<-max((0:length(fit$time))[time<=testtime])
        if (i==0){
             out<-list(surv=1, var=0)
        } else {
             out<-list(surv=fit$KM[i],var=fit$gw[i])
        }
        out
    }
    svp<-getSVar(xp,testtime)
    sv1<-getSVar(x1,testtime)
    sv2<-getSVar(x2,testtime)

    trans<-match.arg(trans)

    if ((trans!="identity" | !varpooled) & correct) stop("correction only defined for trans='identity' with varpooled ")

    ## The order of the following if statements matter. 
    ##   if (doall) then
    ##        -The first if statement creates X2all
    ##        -each if statment following creates X2all<-c(X2all,X2)
    ##        -last if statment does X2<-c(X2all,X2)
    ##    so add methods to the middle
    ##    
    ## In Klein et al 2007, Stat in Med, 4505-4519
    ## Write Greenwood variance as Vhat(Shat(t)) = Shat(t)^2 sigma^2(t)
    ##  so sigma^2(t)  = Vhat/Shat^2

    sigma2.1<- sv1$var/sv1$surv^2
    sigma2.2<- sv2$var/sv2$surv^2
    sigma2<- svp$var/svp$surv^2


    if (doall | (trans=="cloglog" & varpooled)){
        if (sv1$surv==1 | sv1$surv==0 | sv2$surv==1 | sv2$surv==0 | svp$surv==1 | svp$surv==0 | svp$var==0){ 
            X2<-NA
        } else {
            ## Klein et al, eq 8
            X2<- (n1*n2/(n1+n2)^2)*(log(-log(sv2$surv)) - log(-log(sv1$surv)))^2/
                (sigma2/ (log(svp$surv)^2 ) )
        } 
        if (doall) X2all<-X2
    } 
    if (doall | (trans=="cloglog" & !varpooled)){
        if (sv1$surv==1 | sv1$surv==0 | sv2$surv==1 | sv2$surv==0 | svp$surv==1 | svp$surv==0 | svp$var==0){ 
            X2<-NA
        } else {
            ## Klein et al, eq 3
            X2<- (log(-log(sv2$surv)) - log(-log(sv1$surv)))^2/
                (sigma2.1/ (log(sv1$surv)^2 )  + sigma2.2/ (log(sv2$surv)^2 ))
        } 
        if (doall) X2all<-c(X2all,X2)
    } 

    if (doall | (trans=="identity" & varpooled)){
        if (correct){
            # see Fleiss et al 3rd edition, pp. 50-55, rough correction reduces to 
            # usual correction when no censoring, not sure if it makes sense with censoring
            Zpos<- (abs(sv2$surv-sv1$surv)-.5*(1/n1+1/n2))/sqrt(svp$var*(n1+n2)*(1/n1+1/n2))
            X2<-Zpos^2
        } else {
            X2<- ((n1*n2)/(n1+n2)^2)*(
                   (sv2$surv - sv1$surv)^2/(svp$var ) )
        }
        if (doall) X2all<-c(X2all,X2)
    }
    if (doall | (trans=="identity" & !varpooled)){
        if (correct){
          stop("no continuity correction with varpooled=FALSE")
        } else {
            # Klein et al, eq 1
            # note sv$var = Vhat =surv^2 * sigma2
            X2<-   (sv2$surv - sv1$surv)^2/(sv1$var + sv2$var ) 
        }
        if (doall) X2all<-c(X2all,X2)
    }

    if (doall | (trans=="log" & varpooled)){
        if (sv1$surv==0 | sv2$surv==0  | svp$var==0 | svp$surv==1 | svp$surv==0 ){ 
            X2<-NA
        } else {
            ## Klein et al, eq 7
            X2<- (n1*n2/((n1+n2)^2))*(log(sv2$surv) - log(sv1$surv))^2/( sigma2) 
        } 
        if (doall) X2all<-c(X2all,X2)
    } 
    if (doall | (trans=="log" & !varpooled)){
        if (sv1$surv==0 | sv2$surv==0  | svp$var==0 | svp$surv==1 | svp$surv==0 ){ 
            X2<-NA
        } else {
            ## Klein et al, eq 2
            X2<-   (log(sv2$surv) - log(sv1$surv))^2/
                       ( sigma2.1 + sigma2.2) 
        } 
        if (doall) X2<-c(X2all,X2)
    } 


    Z<- sign(sv2$surv - sv1$surv)*sqrt(X2)            
    plo<- pnorm(Z)
    phi<- 1- pnorm(Z)
    p2<- pmin(1,2*plo,2*phi)
    if (doall){ pnames<-paste(rep(c("cloglog","identity","log"),each=2),
                        rep(c("pooled","unpooled"),3),sep=":")
    } else pnames<-ifelse(varpooled,
                       paste(trans,"pooled",sep=":"),
                       paste(trans,"unpooled",sep=":"))
    names(plo)<-names(phi)<-names(p2)<-pnames 
    list(plo=plo,phi=phi,p2=p2)
}




