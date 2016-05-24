

nimDiffOR<-function(p,delta=0.10, q=0.20){
    or<-((q+delta)*(1-q))/(q*(1-q-delta))
    # this is like NiM3 in Rohmel and Kieser,  2013: 2335-2348
    # with q=.2, and delta=.1 and using difference from t>q 
    I<-p>q
    if (delta>=0){
        out<-pmin(1,p+delta,or*p/(1+(or-1)*p))
        out[I]<-pmin(1,p[I]+delta)
    } else {
        out<-pmax(0,p+delta,or*p/(1+(or-1)*p))
        out[I]<-pmax(0,p[I]+delta)
    }
    ## comuputer rounding may give values outside possible range 
    out[out<0]<-0
    out[out>1]<-1
    out
}


nimOR<-function(p,delta=.1,q=.2){
    or<-((q+delta)*(1-q))/(q*(1-q-delta))
    out<-or*p/(1+(or-1)*p)
    ## comuputer rounding may give values outside possible range 
    out[out<0]<-0
    out[out>1]<-1
    out
}

nimDiff<-function(p,delta=.1, q=NULL){
    # include q even though do not need it to avoid errors in program
    out<-pmin(1,p+delta)
    out[out<0]<-0
    out
}



nicqControl<-function(rdig=5,slowint=FALSE,mint=100,interr=10^-3,
    epsilon=10^(-4),alpha=0.025,tau.conf.level=0.95){
    # add some checks later
    if (tau.conf.level>1 | tau.conf.level<0) stop("tau.conf.level should be between 0 and 1")
    list(rdig=rdig,slowint=slowint,mint=mint,interr=interr,
        epsilon=epsilon,alpha=alpha,
        tau.conf.level=tau.conf.level)
}

getfx2<-function(i,n1,n2,g,slowint=FALSE,mint=100,interr=10^-3){
    # the boundary class of distributions has
    # F2(t)  = g(F1(t))
    #    letting u=F1(t), we get 
    # F2(F1inv(u)) = g(u)
    intfunc<-function(x,x2=0){
        dbinom(x2,n2,g(x))* dbeta(x,i,n1-i+1)
    }
    fx2<-rep(NA,n2+1)
    if (slowint){
     for (j in 0:n2){
        # since the intfunc may be near zero over 
        # most of the range of 
        # integration, first calculate intfunc 
        # over many places to 
        # find the maximum, then integrate from 
        # 0 to max, and from max to 1
        ui<-1:mint/(mint+1)
        intvalues<-intfunc(ui,x2=j)
        maxu<- ui[intvalues==max(intvalues)][1]
        if (is.na(maxu)) stop("maxu is missing, search 
                getfx2 for maxu to debug program")

        int1<-integrate(intfunc,lower=0,upper=maxu,x2=j)$value
        # if you get errors in integrate, use following 2 
        # lines for debugging
        #int2<-try(integrate(intfunc,lower=maxu,
        #     upper=1,x2=j)$value)
        #if (class(int2)=="try-error") browser()
        int2<-integrate(intfunc,lower=maxu,upper=1,x2=j)$value
        fx2[j+1]<-int1+int2
      }
    } else {
     for (j in 0:n2){
        #fx2[j+1]<-integrate(intfunc,lower=0,
        #       upper=1,x2=j)$value
        m<-mint
        ## the following is the same as taking the average of 
        ##     sum( (1/m)*intfunc(1:m/m,x2=j) )      and 
        ##     sum( (1/m)*intfunc((0:(m-1))/m,x2=j) ) 
        fx2[j+1]<- sum( (1/m)*intfunc((1:(m-1))/m,x2=j) )+0.5*sum( (1/m)*intfunc(c(0,m)/m,x2=j) )
      }
    }
    # to get rid of small integration errors, standardize fx2 to sum to one
    if (abs(sum(fx2)-1)>interr){
      stop("integration error greater than interr, increase interr or use slowint=TRUE")
    }
    fx2<-fx2/sum(fx2)
    fx2
}

x1<-getfx2(i=7,n1=20,n2=30,g=nimDiff,slowint=FALSE,mint=100,interr=10^-3)
x2<-getfx2(i=7,n1=20,n2=30,g=nimDiffOR,slowint=FALSE,mint=100,interr=10^-3)

nicqGetpvalue<-function(X2T1i,n1,n2,g,i0,alternative="less",control=nicqControl()){
    fx2<-getfx2(i0,n1,n2,g,
        slowint=control$slowint,
        mint=control$mint,
        interr=control$interr)
    X2<-0:n2
    if (alternative=="less"){
        pvals<- cumsum(fx2)
    } else {
        pvals<- rev(cumsum(rev(fx2)))
    }    
    p.value<-pvals[X2==X2T1i]
    p.value
}

#nicqGetpvalue(X2T1i=4,n1=10,n2=20,g=function(x){ nimDiffOR(x,q=.2,delta=.1)},i0=2,alternative="less",control=nicqControl())




nicqCalc<-function(X2T1i,g,i0,n1,n2,delta0=.10,q0=.20,
    conf.int=TRUE,
    conf.level=0.95,
    conf.sided=c("two.sided","one.sided"), 
    alternative="less",
    control=nicqControl()){
    ## to keep from recursive calls, rename variables before calling nicq
    q<-q0
    conf.sided<-match.arg(conf.sided)
    ALT<-alternative
    I0<-i0
    p.value<-nicqGetpvalue(X2T1i,n1,n2,g=function(p){
                  g(p,q=q0,delta=delta0)},
                  i0=I0,
                  alternative=ALT,
                  control=control)

    if (conf.int){
        alpha<-1-conf.level
        if (conf.sided=="two.sided") alpha<-alpha/2
        epsilon<-control$epsilon
        func<-function(Delta,alt="less"){
            alpha - nicqGetpvalue(X2T1i,n1,n2,g=function(p){
                  g(p,q=q0,delta=Delta)},
                  i0=I0,
                  alternative=alt,
                  control=control)
        }
        # the range for delta is: 
        #          -q <= delta <= 1-q
        ci<-c(-q,1-q)
        if (conf.sided=="two.sided" | ALT=="less"){
            if (func(1-q-epsilon,alt="less")<0){
                ci[2]<-1-q
            } else if (func(-q+epsilon,alt="less")>0){
                # all values of delta0 give p-values less than alpha
                # so upper limit is lowest possible
                ci[2]<- -q+epsilon
            } else { 
                ci[2]<-uniroot(func,c(-q+epsilon,1-q-epsilon),alt="less")$root
            }
        }
        if (conf.sided=="two.sided" | ALT=="greater"){
            if (func(-q+epsilon,alt="greater")<0){
                ci[1]<--q
            } else if (func(1-q-epsilon,alt="greater")>0) {
               # all values of delta0 give p-values less than alpha
               # so lower limit is highest possible
                ci[1]<- 1-q-epsilon
            } else { 
                ci[1]<-uniroot(func,c(-q+epsilon,1-q-epsilon),alt="greater")$root
            }
        }
        attr(ci,"conf.level")<-conf.level
        attr(ci,"sided")<-conf.sided
    } else {
        ci<-NULL
    }

    names(X2T1i)<- "x2=x_2(t_1(i))"
    parm<-c(q,i0,n1,n2)
    names(parm)<-c("q","i","n1","n2")
    est<- c(X2T1i/n2,i0/n1,X2T1i/n2-i0/n1)

    names(est)<-c("x2/n2","i/n1","x2/n2-i/n1")

    names(delta0)<-paste0("F2(tau)-F1(tau) \n [where F1(tau)=",q0,"]") 
    #names(delta0)<-"F2(tau)-F1(tau)"  
    attr(delta0,"quantile")<-q0
 
    #names(delta0)<-"difference in proportions (at qth control quantile)"    

    METHOD<-"Non-inferiority test for variable margin using control quantiles"
    METHOD<-paste0(METHOD,"\n with F1(tau)=",q)
    out<-list(statistic=X2T1i,
        parameter=parm,
        p.value=p.value,
        conf.int=ci,
        estimate=est,
        null.value=delta0,
        alternative=alternative,
        method=METHOD)
    class(out)<-"htest"
    out
}

#nicqCalc(X2T1i=100,g=nimDiffOR,i0=20,n1=100,n2=300,delta0=.10,q0=.20,
#    conf.sided=c("two.sided","one.sided"),conf.level=0.95,alternative="less",
#    control=nicqControl())




getimaxpower<-function(n1,n2,g,alpha=0.05,rdig=5,maxprop=1,alt="less",...){
    #  find the i that gives the maximum power when F1=F2
    a<-function(u){u}
    ii<-1:round(n1*maxprop)
    power<-rep(NA,length(ii))
    for (i in 1:length(ii)){
        fx2null<-getfx2(ii[i],n1,n2,g,...)
        if (alt=="less"){
            pvalues<-cumsum(fx2null)
        } else if (alt=="greater"){
            # do cumsum from top down, so rev, then rev to get back
            pvalues<-rev(cumsum(rev(fx2null)))
        } else stop("alt should be 'less' or 'greater' ")
        fx2alt<-getfx2(ii[i],n1,n2,a,...)
        reject<-rep(FALSE,n2+1)
        reject[pvalues<=alpha]<-TRUE
        power[i]<- sum(fx2alt[reject])
    }
    names(power)<-ii
    # in case there are very small differences in power (perhaps due to 
    # computer rounding, eliminate them by rounding to the nearest rdig digits
    rpow<-round(power,rdig)
    imaxpow<-min(ii[rpow==max(rpow)])
    #list(power=power,imaxpow=imaxpow)
    imaxpow
}





nicqTest<-function(x,delta0,q,g=nimDiffOR,yc=NULL,nc=NULL,nt=NULL,ic="prop",
    z=NULL,status=NULL,ties=c("cons","approx"),
    alternative=c("less","greater"),
    conf.level=0.95,
    conf.int=TRUE,
    conf.sided=c("two.sided","one.sided"), 
    gname=NULL,
    control=nicqControl()){
    alternative <- match.arg(alternative)
    ties <- match.arg(ties)
    # first figure out the input format
    # we can enter the data in 4 ways:
    #
    #   input type 1)
    #      x= vector of failures in both groups
    #      z= vector of group membership, 1=control, 2=test
    #      status= 1 is observed, 0 is censored, if
    #          NULL set to all 1s
    #
    #   input type 2)
    #      x= vector of failures on test group
    #      yc = vector of failures in control group
    #
    #   input type 3)
    #      x= number of failures in test group that have 
    #          occured by the 
    #         ic^th failure in control group
    #      ic= rank of the failure of control group upon 
    #           which to base test 
    #          (ignores issues of ties). 
    #         (can also be "prop" or "maxpower", then it 
    #         is calculated)
    #      nc=number in control group
    #      nt=number in test group
    #
    input.type<-0
    if (!is.null(z)){
        #################
        # input type 1:
        #################
        input.type<-1
        # make sure z is the correct format
        if (length(z)!=length(x)) stop("z must be the same length as x")
        if (!all(z==1 | z==2)) stop("elements of z must be 1 (for control) or 2 (test)")
        # create yt and yc
        yt<- x[z==2]
        if (!is.null(yc)) stop("yc input must be NULL if z is not NULL")
        yc<-x[z==1]
    } else if (!is.null(yc)){
        #################
        # input type 2:
        #################
        input.type<-2
        yt<-x
    }
    tauhat<-NA
    if (input.type==1 | input.type==2){
        # calc nc and nt
        nc<-length(yc)
        nt<-length(yt)
        if (ic=="prop"){ 
            ic<- ceiling(nc*q) 
        } else if (ic=="maxpower"){
            warning("when using ic='maxpower' the confidence intervals on the difference are suspect")
            g0<-function(p){ g(p,delta=delta0,q=q) } 
            ic<-getimaxpower(nc,nt,g0,alpha=control$alpha,rdig=control$rdig,alt=alternative)
        }
        ycic<- sort(yc)[ic]
        if (floor(nc*q)<=ic & ic<=ceiling(nc*q)){ tauhat<-ycic }
        # check for ties at ycic
        mc<-length(yc[yc==ycic])
        mt<-length(yt[yt==ycic])
        if (mc>1 | mt>0){
            if (ties=="cons"){
                if (alternative=="less"){
                    x2t1i<-length(yt[yt<=ycic])
                } else {
                    x2t1i<-length(yt[yt<ycic])
                }
            } else if (ties=="approx"){ 
               x1alm1<- length(yc[yc<ycic])
               x2t1i.temp<- length(yt[yt<ycic])+(ic-x1alm1)*mt/(mc+1)
               special.rounding<-function(x,alt=alternative){
                   # since x is an iteger+a ratio of integers
                   # there should not be any machine level errors
                   if (x-floor(x)==.5){
                       if (alt=="less"){
                           out<-ceiling(x)
                       } else {
                           out<-floor(x)
                       }
                   } else {
                       out<-round(x)
                   }   
                   out
               }
               x2t1i<-special.rounding(x2t1i.temp)
            }
        }  else {
            x2t1i<-length(yt[yt<=ycic])
        }     
        # check for censoring problems
        # i.e., censoring before ycic
        if (!is.null(status)){
            if (is.null(z)) stop("if status is not NULL then z must be not NULL")
            if (any(x<=ycic & status==0)) stop("test undefined with censored values at or before the ic^th failure in the control group")
        }
        # estimate the qth quantile of the control group
        if (is.null(status)){
            statusc<-rep(1,nc)
        } else {
            statusc<-status[z==1]
        }
        # get estimate of tau and confidence interval, 
        # use bpcp R package
        # since bpcp cannot use negative values, 
        # if there are any then 
        # trick bpcp by adding a constant and then 
        # subtracting it at the end
        Add<-0
        if (min(yc)<=0){
            Add<- ceiling(abs(min(yc)))+1   
        }
        ycAdd<-yc+Add
        bout<-bpcp(ycAdd,statusc, 
             alpha=1-control$tau.conf.level)
        tauVector<-c(quantile(bout,1-q)[2:4]-Add,
                control$tau.conf.level)
    }  else if (is.null(yc) & is.null(z)){
        #################
        # input type 3:
        #################
        # tau cannot be estimated if do not input yc
        tauVector<-rep(NA,4)
        #warning("Input style assumes no ties in failures")
        if (length(x)!=1) stop("if yc=NULL and z=NULL then 
                    x should be a one-dimensional vector")
        # copied from help in is.integer
        is.wholenumber <-function(x, 
           tol = .Machine$double.eps^0.5){
            abs(x - round(x)) < tol
        }

        integer.check<-function(x){
              if (is.null(x) | !is.numeric(x)){
                    stop(paste0("if yc=NULL and z=NULL then ",
                      deparse(substitute(x)),
                      " should be numeric"))
              }
              if (!is.wholenumber(x) | x<1){
                    stop(paste0("if yc=NULL and z=NULL then ",
                       deparse(substitute(x)),
                      " should be a positive integer"))
              }
              round(x)
        }
        nc<-integer.check(nc)
        nt<-integer.check(nt)

        ## get numeric ic if it is "prop" or "maxpower"
        ## otherwise it should be an integer
        if (ic=="prop"){ 
            ic<- ceiling(nc*q) 
        } else if (ic=="maxpower"){
            warning("when using ic='maxpower' the 
                 confidence intervals on the 
                 difference are suspect")
            g0<-function(p){ g(p,delta=delta0,q=q) } 
            ic<-getimaxpower(nc,nt,g0,
                    alpha=control$alpha,
                    rdig=control$rdig,
                    alt=alternative)
        } else {
            ic<-integer.check(ic)
            if (ic!=ceiling(nc*q)) warning("ic not equal 
                    to ceiling(nc*q), may give difference in
                    proportions outside of 
                    confidence interval.")
        }
        x<-integer.check(x)
        x2t1i<-x

    }  else {
        stop("arguments not defined correctly, 
              check z,yc,ic,nc,nt")
    }
    
    names(tauVector)<-c("tau","lower CL",
        "upper CL","conf.level")

    #######################################
    #   Main Calculation Function
    #######################################
    CSIDE<-match.arg(conf.sided)
    CINT<- conf.int
    CLEVEL<- conf.level
    ALT<-alternative
    out<-nicqCalc(x2t1i,g,ic,nc,nt,delta0,q,conf.int=CINT,
          conf.level=CLEVEL,conf.sided=CSIDE,alternative=ALT)

    out$estimate<-c(out$estimate,tauVector)
    # add details to out$method
    if (is.null(gname)){
        gname<-deparse(substitute(g))
    }
    if (input.type==1 | input.type==2){
        if (ties=="cons"){
            tieStatement<-" and a conservative adjustment for ties"
        } else if (ties=="approx"){
            tieStatement<-" and an approximate adjustment for ties"
        }
        if (length(unique(c(yt,yc)))==nc+nt){
            ## no ties
            tieStatement<-"  (there are no ties in responses)"
        }

    } else { tieStatement<-" and assuming no tied responses" }

    out$method<-paste0(out$method," with variable margin function=",gname,
         tieStatement)

    class(out)<-c("nicq","htest")
    out
}

print.nicq<-function(x,...){
    printHtest<-getS3method("print","htest")
    est<-x$estimate
    x$estimate<-NULL
    printHtest(x,...)
    call<-match.call(expand.dots=TRUE)
    if (is.null(call$digits)){
        digs<-getOption("digits")
    } else {
        digs<-call$digits
    }
    cat("proportions and difference:\n")
    print(est[1:3],digits=digs)
    cat("\n")
    if (!is.na(est[4])){
        cat(paste0("qth quantile of control with ",
           100*est["conf.level"],"% confidence interval:\n"))
        print(est[4:6],digits=digs)
    }
    invisible(x)
}



powerNicqTest<-function(n1=NULL,n2=NULL,power=NULL,sig.level=0.025,
     n2.over.n1=1, q=.20,delta0=.10,
     alternative=c("less","greater"),gnull=nimDiffOR,galt=function(x){x},
     minn=5,maxn=10^5,...){

    alternative<-match.arg(alternative)

    if (is.null(power)){
         if (is.null(n2) & !is.null(n1)){ 
             n2<- ceiling(n1*n2.over.n1)
         } else if (is.null(n1)){
            stop("if power=NULL then n1 must be supplied")
         }
    } else if (!is.null(power) & !is.null(n1)){
        stop("cannot supply both power and n1")
    }

    g0<-function(p){ gnull(p,delta=delta0,q=q) } 
    #####################################
    # Calculate power given n1 and n2
    #####################################
    calcpower<-function(n1,n2){
        i<- ceiling(q*n1)
        fx2null<-getfx2(i,n1,n2,g0)
        if (alternative=="less"){
            pvals<- cumsum(fx2null)
        } else {
           pvals<- rev(cumsum(rev(fx2null)))
        } 
        reject<- pvals<=sig.level
        fx2alt<-getfx2(i,n1,n2,galt)
        power<- sum(fx2alt[reject])
        power
    }


    if (is.null(power)){
        n1<-ceiling(n1)
        n2<-ceiling(n2)
        power<- calcpower(n1,n2)
    } else if (is.null(n1)){
        rootfunc<-function(N1){
            power - calcpower(N1,ceiling(n2.over.n1*N1))
        }
        urout<-uniroot.integer(rootfunc,c(minn,maxn),...)
        n1<-urout$root
        n2<-ceiling(n2.over.n1*n1)
        # f.root=nominal.power - calcpower
        # so calcpower= nominal.power-f.root
        power<- power - urout$f.root
    }

    #NOTE <- paste0("q=",q)
    #NOTE<-""
    METHOD <- c("Non-inferiority control quantile test power")
    structure(list(n1 = n1, n2=n2, delta0 = delta0, q=q, sig.level = sig.level, 
        power = power, alternative = alternative, 
        method = METHOD), class = "power.htest")

}
#powerNicqTest(n1=345)
#powerNicqTest(n1=346)
## find minimum n1 that have power greater than 0.80
## takes 13 iterations to find n1=346
## so do not run it here
## powerNicqTest(power=.80,print.steps=TRUE)
