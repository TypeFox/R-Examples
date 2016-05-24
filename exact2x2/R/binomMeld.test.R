binomMeld.test <-
function(x1,n1,x2,n2,nullparm=NULL, 
    parmtype=c("difference","oddsratio","ratio"),
    conf.level=0.95, conf.int=TRUE,
    alternative=c("two.sided","less","greater"),eps=10^-8){

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
    # T1 ~ Beta(x1+1,n1-x1)
    # T2 ~ Beta(x2,n2-x2+1)
    #    and p-value is calculated under null: Pr[ g(T1,T2) <= D ]
    #
    #    for less (calculate upper CL) we use 
    # T1 ~ Beta(x1,n1-x1+1)
    # T2 ~ Beta(x2+1,n2-x2)
    #    and p-value is calculated under null: Pr[ g(T1,T2) >= D ]

    if (ptype=="difference"){ 
        # Pr[ g(T1,T2) >=D] = Pr[ T2-T1 >= D] = Pr[ T1 <= T2 - D]
        # = \int F1(t2 - D) f2(t2)  dt2
        funcLess<-function(t2,D){
            pbeta(t2-D,x1,n1-x1+1)*dbeta(t2,x2+1,n2-x2)
        }
        # Pr[ g(T1,T2) <=D] = Pr[ T2-T1 <= D] = Pr[ T2 <= T1 + D]
        # = \int F2(t1 + D) f1(t1)  dt1
        funcGreater<-function(t1,D){
            pbeta(t1+D,x2,n2-x2+1)*dbeta(t1,x1+1,n1-x1)
        }
    } else if (ptype=="ratio"){
        # Pr[ g(T1,T2) >=D] = Pr[ T2/T1 >= D] = Pr[ T1 <= T2/D]  
        # = \int F1(t2/D) f2(t2)  dt2
        funcLess<-function(t2,D){
            pbeta(t2/D,x1,n1-x1+1)*dbeta(t2,x2+1,n2-x2)
        }
        # Pr[ g(T1,T2) <=D] = Pr[ T2/T1 <= D] = Pr[ T2 <= T1*D]
        # = \int F2(t1 * D) f1(t1)  dt1
        funcGreater<-function(t1,D){
            pbeta(t1*D,x2,n2-x2+1)*dbeta(t1,x1+1,n1-x1)
        }
    } else if (ptype=="oddsratio"){
        # Pr[ g(T1,T2) >=D] = Pr[ T2(1-T1)/T1(1-T2) >= D] 
        # = Pr[ T2(1-T1) >= T1(1-T2)D] = Pr[ T2 >= T1*(T2 + (1-T2)D) ]
        # = Pr[ T1 <= T2/{ T2 + (1-T2)D } ]  
        # = \int F1(t2/(t2+(1-t2)D) f2(t2)  dt2
        funcLess<-function(t2,D){
            pbeta(t2/(t2+(1-t2)*D),x1,n1-x1+1)*dbeta(t2,x2+1,n2-x2)
        }
        # Pr[ g(T1,T2) <=D] = Pr[ T2(1-T1)/T1(1-T2) <= D] 
        # = Pr[ T2(1-T1) <= T1(1-T2)D] = Pr[ T2(1-T1) + T1*T2*D <= T1*D  ]
        # = Pr[ T2( (1-T1) + T1*D) <= T1*D ]
        # = Pr[ T2 <= T1*D/{(1-T1) + T1*D} ]
        # = \int F2(t1*D/(1-t1+t1*D)) f1(t1)  dt1
        funcGreater<-function(t1,D){
            pbeta(t1*D/(1-t1+t1*D),x2,n2-x2+1)*dbeta(t1,x1+1,n1-x1)
        }
    }


    # p-value functions 
    pGreater<-function(delta){
        ## for the integrate function to work well, pick values that 
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
        ## First, choose LowerInt=a and UpperInt=b so that \int_a^b  dbeta1(t) dt = 1- eps/2 
        LowerInt<-qbeta(eps/4,x1+1,n1-x1)
        UpperInt<-qbeta(1-eps/4,x1+1,n1-x1)
        ##  Second, choose a2 so that pbeta2( W[a2] )= eps/2
        ##       or   qbeta2( eps/2) = W[a2] 
        q<- qbeta(eps/2, x2,n2-x2+1)
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
                pout[i]<-integrate(funcGreater,LowerInt,UpperInt,D=delta[i])$value
            }
        }
        ## since we underestimate the integral (assuming perfect integration) by at most eps, 
        ## add back eps to get conservative p-value
        pout<-pout+eps
        pout
    }
    pLess<-function(delta){
        ## for the integrate function to work well, pick values that 
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
        ## First, choose LowerInt=a and UpperInt=b so that \int_a^b  dbeta2(t) dt = 1- eps/2 
        LowerInt<-qbeta(eps/4,x2+1,n2-x2)
        UpperInt<-qbeta(1-eps/4,x2+1,n2-x2)
         
        ##  Second, choose a1 so that pbeta1( W[a1] )= eps/2
        ##       or   qbeta1( eps/2) = W[a1] 
        q<- qbeta(eps/2, x1,n1-x1+1)
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
                pout[i]<-integrate(funcLess,LowerInt,UpperInt,D=delta[i])$value
            }
        }
        ## since we underestimate the integral (assuming perfect integration) by at most eps, 
        ## add back eps to get conservative p-value
        pout<-pout+eps
        pout
    }
    lower<-upper<-NA
    alt<-match.arg(alternative)
    if (alt=="two.sided"){
        dolo<-dohi<-TRUE
        alpha<-(1-conf.level)/2  
    } else if (alt=="less"){
        ## alt=less so lower interval is lowest possible, do not calculate
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

    


    # use MLE, maybe other estimates may make sense but do not use them now
    estimate<-g(x1/n1,x2/n2)
    if (ptype=="difference"){
        names(estimate)<-"difference (p2-p1)"
    } else if (ptype=="ratio"){
        names(estimate)<-"ratio (p2/p1)"
    } else if (ptype=="oddsratio"){
        names(estimate)<-"odds ratio {p2(1-p1)}/{p1(1-p2)}"
    }

    if (dolo){
        ## Take care of special cases, when x2=0 T2 is point mass at 0
        ## when x1=n1 T1 is a point mass at 1
        if (x2==0 & x1<n1){
            if (conf.int) lower<- g( qbeta(1-alpha,x1+1,n1-x1), 0 )
            if (ptype=="difference"){ 
                pg<- 1- pbeta(-nullparm,x1+1,n1-x1)
            } else pg<-1
        } else if (x1==n1 & x2>0){
            if (conf.int) lower<- g( 1, qbeta(alpha,x2,n2-x2+1) )
            if (ptype=="difference"){
                ## Aug 22, 2014: Fixed following line
                ## WRONG: pg<- 1-pbeta(1+nullparm,x2,n2-x2+1)
                pg<- pbeta(1+nullparm,x2,n2-x2+1)
            } else if (ptype=="ratio"){
                pg<-pbeta(nullparm,x2,n2-x2+1)
            } else if (ptype=="oddsratio"){
                ## Aug 22, 2014: Fixed following line
                ## WRONG:pg<-pbeta(nullparm/(1-nullparm),x2,n2-x2+1)
                pg<-1
            }
        } else if (x2==0 & x1==n1){
            if (conf.int) lower<- g( 1, 0 )
            pg<-1
        } else {
            if (conf.int){ 
                rootfunc<-function(delta){
                    pGreater(delta)-alpha
                }
                if (upperLimit==Inf){
                    # uniroot cannot take Inf as an upper limit
                    # find T such that rootfunc(T) has opposite sign as
                    # rootfunc(lowerLimit)
                    for (i in 1:50){
                        upperLimit<-2^i
                        if (sign(rootfunc(lowerLimit))!=sign(rootfunc(upperLimit))){
                            break()
                        }
                    }
                }
                if (upperLimit==2^50){
                    warning("lower conf limit appears to be larger than 2^50=approx=10^16, set to 2^50")
                    lower<-2^50
                } else {
                    lower<-uniroot(rootfunc,c(lowerLimit,upperLimit))$root
                }
            }
            pg<-pGreater(nullparm)
            #reset upperLimit
            upperLimit<-g(0,1)
        }
    }
    if (dohi){
        ## Take care of special cases, when x2=n2 T2 is point mass at 1
        ## when x1=0 T1 is a point mass at 0
        if (x1==0 & x2==n2){
            if (conf.int) upper<-g(0,1)
            pl<-1
        } else if (x1==0){
            if (conf.int) upper<- g(0,qbeta(1-alpha,x2+1,n2-x2))
            if (ptype=="difference"){
               pl<-1-pbeta(nullparm, x2+1,n2-x2)
            } else if (ptype=="ratio"){
               pl<- 1
            } else if (ptype=="oddsratio"){
               pl<-1
            }
        } else if (x2==n2){
            if (conf.int) upper<- g(qbeta(alpha,x1,n1-x1+1), 1)
            if (ptype=="difference"){
                ## Aug 22, 2014: Fixed following line
                ## WRONG:pl<-pbeta(1+nullparm, x1, n1-x1+1)
                pl<-pbeta(1-nullparm, x1, n1-x1+1)
            } else if (ptype=="ratio"){
               pl<-pbeta(1/nullparm,x1,n1-x1+1)
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
                    # find T such that rootfunc(T) has opposite sign as
                    # rootfunc(lowerLimit)
                    for (i in 1:50){
                        upperLimit<-2^i
                        if (sign(rootfunc(lowerLimit))!=sign(rootfunc(upperLimit))){
                            break()
                        }
                    }
                }
                if (upperLimit==2^50){
                    warning("upper conf limit appears to be larger than 2^50=approx=10^15, set to Inf")
                    upper<-Inf
                } else {
                    upper<-uniroot(rootfunc,c(lowerLimit,upperLimit))$root
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
    dname<-paste("sample 1:(",x1,"/",n1,"), sample 2:(",x2,"/",n2,")",sep="")
    method<-paste("melded binomial test for",ptype)
    stat<-x1/n1
    parm<-x2/n2
    names(stat) <- "proportion 1"
    names(parm) <- "proportion 2"
    names(nullparm)<-ptype

    structure(list(statistic = stat, parameter = parm, 
        p.value = p.value, 
        conf.int = ci, estimate = estimate, null.value = nullparm, 
        alternative = alt, method = method, 
        data.name = dname), class = "htest")

}
