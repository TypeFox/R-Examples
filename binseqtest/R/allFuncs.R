

#########################################################
# general functions needed 
#########################################################


unirootDiscrete<-function (f, interval, lower = min(interval), upper = max(interval), 
    tol=10^-5,pos.side = FALSE, print.steps = FALSE,
    maxiter = 1000, ...) 
{
    iter <- 0
    if (!is.numeric(lower) || !is.numeric(upper) || lower >= 
        upper) 
        stop("lower < upper  is not fulfilled")
    step.power<-ceiling(log2(((upper-lower)/tol)+1))
    N<-2^step.power

    #getx<-function(i,n=N,Lower=lower,Upper=upper){ Lower+(Upper-Lower)*(i)/(n)}
    xlow <- lower
    f.low <- f(lower, ...)
    iter <- iter + 1
    if (print.steps) {
        print(paste("x=", xlow, " f(x)=", f.low))
    }

    xup<-upper
    f.up<-f(xup,...)
    iter<-iter+1
    if (print.steps) {
        print(paste("x=", xup, " f(x)=", f.up))
    }

    if (sign(f.low*f.up)!=-1) stop("f(lower) must be opposite sign as f(upper)")

    step.power<-step.power-1
    i<-2^step.power
    xmid<-lower + (i/N)*(upper-lower)
    #xmid<-getx(i)
    f.mid<-f(xmid,...)
    if (print.steps) {
        print(paste("x=", xmid, " f(x)=", f.mid))
    }

    while (step.power>-1){ t
        step.power<-step.power-1   
        if (sign(f.low*f.mid)==1){
            f.low<-f.mid
            xlow<-xmid
            i<-i+2^step.power
        }
        else{
            f.up<-f.mid
            xup<-xmid
            i<-i-2^step.power
        }
        xmid<-lower + (i/N)*(upper-lower)
        #xmid<-getx(i)
        f.mid<-f(xmid,...)
        iter<-iter+1
        if (print.steps) {
            print(paste("x=", xmid, " f(x)=", f.mid))
        } 
    }

    f.pos<-max(f.low,f.up)
    xpos<-ifelse(f.pos==f.low,xlow,xup) 
    f.neg<-min(f.low,f.up)
    xneg<-ifelse(f.neg==f.low,xlow,xup)

    if (pos.side) {
        root <- ifelse(f.mid > 0, xmid, xpos)
        f.root <- ifelse(f.mid > 0, f.mid, f.pos)
    }
    else {
        root <- ifelse(f.mid < 0, xmid, xneg)
        f.root <- ifelse(f.mid < 0, f.mid, f.neg)
    }
    list(iter = iter, f.root = f.root, root = root)
}




#########################################################
#
#   abparms class
#
#########################################################

# validity check for abparms objects
validAbparms<-function(object){
    out<-TRUE
    Nk<-object@Nk
    a<-object@a
    b<-object@b
    # only give one error message at a time, start with simplest

    

    if (length(Nk)!=length(a) & length(Nk)!=length(b)){
        out<-"length of Nk, a, and b should be the same"
    } else if (any(is.na(Nk))){
        out<-"All values of Nk must be non-missing"
    } else if (any(Nk<1) | ((length(Nk)>1) && any(diff(Nk)<=0))){
        out<- "Nk values should be greater than 0 and increasing"
    } else if (any(!is.na(a))){
        # check 'a' parameters (note a==NA means do not stop there)
        acheck<-a[!is.na(a)]
        Ncheck<-Nk[!is.na(a)]
        if (any(acheck<0) | any(acheck>=Ncheck)) out<-" 'a' values should be greater than or equal to 0 and less than Nk, use NA for not stopping"
        # all 'a' should be increasing
        if (length(acheck)>1 && any(diff(acheck)<=0)){
            out<-" 'a' values should be increasing (NA allowed)"
        }    
    } else if (any(!is.na(b))){
        # check 'b' parameters (note b==NA means do not stop there)
        bcheck<-b[!is.na(b)]
        Ncheck<-Nk[!is.na(b)]
        if (any(bcheck>Ncheck) | any(bcheck<=0) ) out<-" 'b' values should be less than or equal Nk and greater than 0, use NA for not stopping"
        # all 'b' can be increasing or decreasing, but cannot increase faster than Nk
        if (length(bcheck)>1 && any(diff(bcheck)>diff(Ncheck))){
            out<-" 'b' values cannot increase faster than Nk (NA allowed)"
        }    
    } else if (any(!is.na(a) & !is.na(b))) {
        nomiss<- !is.na(a) & !is.na(b)
        if (any(b[nomiss]-a[nomiss]<0)){
            out<-" 'b' values should be greater than 'a' values"
        }
    } else {
        # no errors    
    }
    out
}

setClass("abparms",
    representation(Nk="numeric",a="ANY",b="ANY",binding="character"),
    prototype(Nk=50,a=as.numeric(NA),b=as.numeric(NA),binding="both"),
    validity=validAbparms)

#########################################################
#
#   bound class
#
#########################################################

validBound<-function(object){
    out<-TRUE
    ## check that probability sum to 1, should work with any theta
    check<-function(object,theta=.41231){
        checksum<-sum( object@K* theta^object@S * (1-theta)^(object@N-object@S) )
        checksum
    }
    #if (length(object@K)==0 | length(object@S)==0 | length(object@N)==0 | length(object@UL)==0
    if (abs(1-check(object))> 10^-8) out<-"sum of probabilities not within 10^-8 of 1"
    out
}
setClass("bound",
    representation(S="numeric",N="numeric",K="numeric",order="numeric",UL="character"),
    contains="abparms",
    validity=validBound)


#########################################################
#
#   boundNBF class
#    (boundary with non-binding futility stopping bound)
#
#########################################################
validBoundNBF<-function(object){
    out<-TRUE
    ## check that probability sum to 1, should work with any theta
    check<-function(object,theta=.41231){
        checksum<-sum( object@Ks* theta^object@Ss * (1-theta)^(object@Ns-object@Ss) )
        checksum
    }
    if (abs(1-check(object))> 10^-8) out<-"sum of probabilities not within 10^-8 of 1"
    out
}
setClass("boundNBF",
    representation(Ss="numeric",Ns="numeric",Ks="numeric",orders="numeric",ULs="character"),
    validity=validBoundNBF, contains="bound")


##############################################################
#
#   boundEst class
#     (boundary with Statistics, i.e., with CI, p-values, MUE)
##############################################################
validBoundEst<-function(object){
    out<-TRUE
    out
}
setClass("boundEst",
    representation(estimate="numeric",lower="numeric",upper="numeric",
            tsalpha="numeric",theta0="numeric",
            plower="numeric",pupper="numeric",pval="numeric"),
    contains="bound",
    validity=validBoundEst)


##############################################################
#
#   boundNBFEst class
#     (boundary NBF with Statistics, i.e., with CI, p-values, MUE)
#
#  Stopped using roxygen2 because it does not really support S4 methods at this time
##############################################################

validBoundNBFEst<-function(object){
    out<-TRUE
    out
}
setClass("boundNBFEst",
    representation(estimate="numeric",lower="numeric",upper="numeric",
            tsalpha="numeric",theta0="numeric",
            pvals="matrix",pval="numeric"),
    contains="boundNBF",
    validity=validBoundNBFEst)


#########################################################
#
#   some functions
#
#########################################################




abBindBothCalcK<-function(object){
    ## take abparms object and convert to bound object
    # NOTE: treats object@binding as "both" even if it is "lower" or "upper"
    # Since N is slot in bound object, use Nk for abparms object
    #    but really it is just N
    Nk<-object@Nk
    a<-object@a
    b<-object@b
    Nmax<-max(Nk)
    # same as choose(n,0:n), but slightly faster, I think
    binomCoef<-function(n){
        exp(lchoose(n,0:n))
    }

    binomCoefVector<-function(start,forwardN){
        # start = vector of counts start out with
        # forwardN= number of samples until next stop, i.e., Nk[i]-Nk[i+1]
        startN<-length(start)
        cnt<-rep(0,startN-1+(forwardN+1))
        j<-1
        bc<-binomCoef(forwardN)
        for (i in 1:startN){
            cnt[i:(i+forwardN)]<-cnt[i:(i+forwardN)]+start[i]*bc
        }
        cnt
    }
    ## test it, start from 1, go forward 4
    ## then do not stop, go forward 3 more
    ## should equal going forward 7
    #binomCoefVector(binomCoef(4),3)
    #binomCoef(7)

    # function for calculating which S values stop and which continue 
    # at each Nk[i]
    # upperLower="upper" if upper and "lower" if lower boundary
    Svalues<-function(a,b,Smin,Smax){
        if (is.na(a) & is.na(b)){ 
            Scontinue<-Smin:Smax
            Sstop<-NA
            upperLower<-NA
        } else if (is.na(a)){
            Scontinue<- Smin:(b-1)
            Sstop<-b:Smax
            upperLower<-rep("upper",length(Sstop)) 
       } else if (is.na(b)){
            Scontinue<- (a+1):Smax
            Sstop<- Smin:a 
            upperLower<-rep("lower",length(Sstop))
       } else {
            Scontinue<- (a+1):(b-1)
            Sstop<-c(Smin:a,b:Smax)
            upperLower<-c(rep("lower",length(Smin:a)),rep("upper",length(b:Smax)))
       }
       list(stop=Sstop,continue=Scontinue,upperLower=upperLower)
    }

    # initialize parameters that change each iteraction
    cnt<-1
    Kcontinue<-1
    Scontinue<-0
    
    Nsteps<-diff(c(0,Nk))

    upperLowerTemp<-Ktemp<-Stemp<-Ntemp<-rep(NA,sum(Nk)+Nmax)

    for (i in 1:length(Nsteps)){
        # for each iteration, take Nsteps[i] more samples
        # min and max number of possible successes at the end is Smin and Smax
        # Scontinue is vector of success values at start of iteration
        Smin<-min(Scontinue)
        Smax<-max(Scontinue)+ Nsteps[i]
        # for last step, stop for everything
        if (i==length(Nsteps)){ s<-list(stop=Smin:Smax, continue=NA, upperLower=rep("end",length(Smin:Smax)))
        }else s<-Svalues(a[i],b[i],Smin,Smax)
        # Send is possible S values at end of iteration
        Send<-Smin:Smax
        bcv<- binomCoefVector(Kcontinue,Nsteps[i])
        if (any(!is.na(s$stop))){
            newcnt<-cnt+length(s$stop)
            ii<- cnt:(newcnt-1)
            Stemp[ii]<- s$stop
            upperLowerTemp[ii]<-s$upperLower
            Ntemp[ii]<- rep(Nk[i],length(s$stop))
            istop<- Send %in% s$stop
            if (any(istop)) Ktemp[ii]<- bcv[istop]
            Scontinue<- s$continue
            if (any(!istop)) Kcontinue<- bcv[!istop]
            cnt<-newcnt
        } else {
            Scontinue<- s$continue
            Kcontinue<- bcv
        }
    }
    S<-Stemp[!is.na(Stemp)]
    N<-Ntemp[!is.na(Stemp)]
    K<-Ktemp[!is.na(Stemp)]
    UL<-upperLowerTemp[!is.na(Stemp)]
    # do stagewise ordering:
    #
    #   for lower:  order by N (lowest to highest) then by (S/N) (lowest to highest)
    #        equivalently order by (-Nmax+N) + (.5)*(S/N) (lowest to highest)
    #   for upper:  order by N (highest to lowest), then by (S/N) (lowest to highest)
    #        equivalently order by -(-Nmax+N) + (.5)*(S/N)  (lowest ot highest)
    #   for end: order by S/N (lowest to highest)
    #
    order<-rep(NA,length(S))
    osign<- rep(NA,length(S))
    osign[UL=="lower"]<-1
    osign[UL=="upper"]<- -1
    osign[UL=="end"]<- 0
    orderingStat<-         osign*(-Nmax +N) + (.5)*(S/N)
    # in case object is not of class abparms, just take that part of it 
    if (class(object)!="abparms"){
        object<-new("abparms",Nk=object@Nk,a=object@a,b=object@b,binding=object@binding)
    }
    B<-new("bound",object,S=S,N=N,K=K,order=rank(orderingStat),UL=UL)
    #B<-list(abparms=object,S=S,N=N,K=K,order=rank(orderingStat),UL=UL)
    B
}


abtoBound<-function(from){
    object<-from
    if (!(class(object) %in% c("abparms","bound","boundEst","boundNBF","boundNBFEst"))) stop("need to input object that contains abparms object")
    ## take abparms object and convert to bound object
    binding<-object@binding
    if (binding=="both"){
        B<-abBindBothCalcK(object)
    } else if (binding=="upper"){
        b2sided<-abBindBothCalcK(object)
        # binding=upper means that lower is not binding
        # set a to NA
        object@a<- as.numeric(rep(NA,length(object@a)))
        bu<-abBindBothCalcK(object)
        # boundNBF objects have extra bound for s=superiority
        # where you do not stop for futility, elements are
        # Ss, Ns, orders, ULs
        B<-new("boundNBF",b2sided,Ss=bu@S,Ns=bu@N,Ks=bu@K,orders=bu@order,ULs=bu@UL)
    } else if (binding=="lower"){
        b2sided<-abBindBothCalcK(object)
        # binding=lower means that upper is not binding
        # set b to NA
        object@b<- as.numeric(rep(NA,length(object@b)))
        bl<-abBindBothCalcK(object)
        # boundNBF objects have extra bound for s=superiority
        # where you do not stop for futility, elements are
        # Ss, Ns, orders, ULs
        B<-new("boundNBF",b2sided,Ss=bl@S,Ns=bl@N,Ks=bl@K,orders=bl@order,ULs=bl@UL)
    }
    B
}



pCalc<-function(S,N,K,order,theta0=.5,alternative="two.sided",ponly=FALSE){
    # calculate the density at each boundary point, use exp(log()+log()) to avoid overflow
    probs<- exp( log(K) + S*log(theta0) + (N-S)*log(1-theta0) )

    # to get p-values, you need cumsums using the order
    # so need to put the probs in order calculate the cumsums,
    #  then put back to the original order
    pvalue<-function(probs, order){
        n<-length(order)
        # o lists the index of the smallest, next smallest, etc of `order'
        o<- order(order)
        pcumsum.ordered<- cumsum(probs[o])
        # i is original order, sorted the same way as probs was sorted
        i<- (1:n)[o]
        ## get back to original order
        pcumsum<- pcumsum.ordered[order(i)]
        # check: 
        #data.frame(probs,order,o,probs[o],pcumsum.ordered,pcumsum)
        pcumsum
    }  
    # calculate p-values
    if (ponly){
        #    only calculate the type of p-value you need
        if (alternative=="less"){ 
            pval<- pvalue(probs, order)
        } else if (alternative=="greater"){
            pval<- pvalue(probs, -order)
        } else if (alternative=="two.sided"){
            plower<- pvalue(probs, order)
            pupper<- pvalue(probs, -order)
            # two-sided p-value 
            pval<- pmin(1, 2*pmin(plower, pupper) )
        }
        out<-pval
    } else {
    #   calculate both type of p-values
        plower<- pvalue(probs, order)
        pupper<- pvalue(probs, -order)
        if (alternative=="less"){ 
            pval<- plower
        } else if (alternative=="greater"){
            pval<- pupper
        } else if (alternative=="two.sided"){
            pval<- pmin(1, 2*pmin(plower, pupper) )
        }
        out<-list(plower=plower,pupper=pupper,pval=pval)
    }
    out
}



ciCalc<-function(S,N,K,order,type="upper",alpha=0.025){
      if (type=="two.sided") stop("function only calculates one side at a time")
      # reverse order for lower intervals, otherwise calculations are the same
      if (type=="upper"){ 
          # o lists the index of the smallest, next smallest, etc of `order'
          o<- order(order)
      } else {
          # reverse order for lower Conf limits
          o<- order(-order)
      }
      n<-length(order)
      So<-S[o]
      No<-N[o]
      Ko<-K[o]
      
      # Confidence interval is the root of the following function
      rootfunc<- function(theta, alpha, ii){
          # ii = 1:i
          sum( Ko[ii]*theta^So[ii]*(1-theta)^(No[ii]-So[ii]) ) - alpha 
      }

      ci<-rep(NA,n)
      for (i in 1:(n-1)){
          ci[i]<-uniroot(rootfunc, interval = c(0,1), alpha = alpha,ii=1:i)$root
      }
      # last interval is defined as 1 for upper and 0 for lower 
      if (type=="upper"){
          ci[n]<-1
      } else {
          ci[n]<-0
      }
      # use iorig to get back to original order 
      iorig<- order((1:n)[o])
      CI<- ci[iorig]
      # check: 
      #d<-data.frame(order=order,CI=CI)
      #d[order(d$order),]
      CI
}

#ciCalc(Bs@S,Bs@N,Bs@K,Bs@order,type="lower",alpha=.025)



getTSalpha<-function(tsalpha=NULL,alternative=NULL,conf.level=NULL){
    if (is.null(tsalpha) & is.null(alternative) & is.null(conf.level)){
        stop("Must supply either tsalpha or (alternative and conf.level)")
    } else if (is.null(tsalpha) & !is.null(alternative) & !is.null(conf.level)){
        if (!is.null(alternative) && (alternative!="greater" & alternative!="less" & alternative!="two.sided")){
            stop("alternative should be 'greater', 'less', or 'two.sided' ") 
        }
        if (!is.null(conf.level) && (conf.level>=1 | conf.level<.5)) stop("conf.level should be .5<= conf.level <1") 
        if (alternative=="two.sided"){ 
            onesidedAlpha<- (1-conf.level)/2
            TSalpha<-rep(onesidedAlpha,2)
        } else if (alternative=="less"){
            # alternative is less
            # so you want the CI to be one-sided on the upper side in order to 
            # show significantly less than something
            # that means error on the lower side is 0
            TSalpha<-c(0,1-conf.level)
        } else if (alternative=="greater"){
            # see reasoning for 'less' except opposite
            TSalpha<-c(1-conf.level,0)
        }
        names(TSalpha)<-c("alphaLower","alphaUpper")
        TSalpha
    } else if (!is.null(tsalpha)){
        # tsalpha overrides alternative and conf.level
        TSalpha<-tsalpha
        if (length(TSalpha)!=2) stop("tsalpha should be of length 2")
        if (sum(TSalpha)>1 | sum(TSalpha)<=0 | any(TSalpha<0) | any(TSalpha>.5)){
            stop("tsalpha should represent the errors allowed on either side. We require
                  each element, a,  should have 0<=a<0.5, and the sum should be less than 1")
        }
        if (is.null(names(TSalpha))){ names(TSalpha)<-c("alphaLower","alphaUpper")
        } else if (   !identical(names(TSalpha),c("alphaLower","alphaUpper"))   ){
            stop("names for tsalpha should be either missing or 'alphaLower' and 'alphaUpper' in that order")
        }
    } else {
        stop("either supply tsalpha or (conf.level and alternative)")
    }
    TSalpha
}

#check getTSalpha function
#getTSalpha(NULL,NULL,NULL)
#getTSalpha(c(.025,.05))
#getTSalpha(c(alphaUpper=.025,alphaLower=.05))
#getTSalpha(c(Lower=.025,Upper=.05))
#getTSalpha(alternative="less",conf.level=.95)
#getTSalpha(alternative="greater",conf.level=.95)
#getTSalpha(alternative="two.sided",conf.level=.95)

getAlternative<-function(tsalpha){
    x<-tsalpha
    if (is.null(names(x))){ names(x)<-c("alphaLower","alphaUpper")
    } else if (   !identical(names(x),c("alphaLower","alphaUpper"))   ){
       stop("names for tsalpha should be either missing or 'alphaLower' and 'alphaUpper' in that order")
    }
    if (sum(x)==0 | sum(x)>1){
        stop("invalid TSalpha")
    } else if (x["alphaLower"]==0){
        out<-"less"
    } else if (x["alphaUpper"]==0){
        out<-"greater"
    } else {
        out<-"two.sided"
    } 
    out
}
#getAlternative(c(.025,.05))
#getAlternative(c(.025,0))
#getAlternative(c(0,.05))



analyzeBound<-function(object,theta0=.5,stats="all",alternative="two.sided",conf.level=0.95,tsalpha=NULL,...){
    B<-object
    TSalpha<-getTSalpha(tsalpha,alternative,conf.level)
    alternative<-getAlternative(TSalpha)

    if (stats=="all"){
        pval<-pCalc(B@S,B@N,B@K,B@order,theta0=theta0,alternative=alternative)
        # median unbiased estimator is average of lower and upper CL at .5 level
        Estimate<- .5*ciCalc(B@S,B@N,B@K,B@order,type="upper",alpha=0.5)+
               .5*ciCalc(B@S,B@N,B@K,B@order,type="lower",alpha=0.5)    
        if (TSalpha["alphaLower"]>0){
            lower<-ciCalc(B@S,B@N,B@K,B@order,type="lower",alpha=TSalpha["alphaLower"])
        } else {
            lower<-rep(0,length(B@N))
        }
        if (TSalpha["alphaUpper"]>0){
            upper<-ciCalc(B@S,B@N,B@K,B@order,type="upper",alpha=TSalpha["alphaUpper"])
        } else {
            upper<-rep(1,length(B@N))
        }
        out<-new("boundEst",B,estimate=Estimate,lower=lower,upper=upper,
               tsalpha=TSalpha,theta0=theta0,
               plower=pval$plower,pupper=pval$pupper,pval=pval$pval)
    } else if (stats=="pval"){
        out<-pCalc(B@S,B@N,B@K,B@order,theta0=theta0,alternative=alternative,ponly=TRUE)
    }
    out
}

analyzeBoundNBF<-function(object,theta0=.5,stats="all",alternative="two.sided",
        conf.level=0.95,tsalpha=NULL,cipMatch=TRUE,...){
    # for nonbinding futility bounds, first do the usual analysis for a bound
    # then for the binding bounds and end bounds 
    # replace the p-values and the confidence limit on the opposite side
    # (e.g., for binding on the upper, replace the lower confidence limit) 

    if (object@binding=="both") stop("use analyzeBound when binding='both'")

    
    TSalpha<-getTSalpha(tsalpha,alternative,conf.level)
    alternative<-getAlternative(TSalpha)


    out<-analyzeBound(object, theta0, stats, tsalpha=TSalpha)
    b<-length(out@N)
    # create replacement indicator, rpi, for values to replace
    if (object@binding=="lower"){
        rpi<- out@UL!="upper"       
    } else if (object@binding=="upper"){
        rpi<- out@UL!="lower"       
    } 
    # part of the boundary that needs to be replaced
    Nrp<-object@N[rpi]
    Srp<-object@S[rpi]

    # create a new bound with the superiority boundaries from the NBF
    Bbind<-new("bound",
            Nk=object@Nk,
            a=object@a,
            b=object@b,
            binding=object@binding,
            N=object@Ns,
            S=object@Ss,
            K=object@Ks,
            order=object@orders,
            UL=object@ULs)
    # now do analysis on superiority boundaries
    tempTSalpha<-TSalpha
    if (object@binding=="upper"){
        tempTSalpha["alphaUpper"]<-0
        if (TSalpha["alphaLower"]<=0) stop("binding='upper' but tsalpha[1]=0")
    } else if (object@binding=="lower"){
        tempTSalpha["alphaLower"]<-0
        if (TSalpha["alphaUpper"]<=0) stop("binding='lower' but tsalpha[2]=0")
   }
    outBbind<-analyzeBound(Bbind, theta0, stats, tsalpha=tempTSalpha)
    
    # get unique id for each boundary point
    # basically use N.S for each point, but need to find the 
    # factor of 10 to divide S by before adding onto N
    ndigits<- ceiling(log10( max(c(object@Ss,object@S)) ))
    Sfactor<- 10^ndigits
    # to make sure there is no fuzz problems (i.e., so we can use == or %in% later), 
    # round to nearest ndigits
    NSnew<- round(object@Ns + object@Ss/Sfactor,ndigits)
    NS<- round(Nrp + Srp/Sfactor,ndigits)
    RPI<- NSnew %in% NS
    # check that the N and S match up when you pick out values from each boundary
    if (!all( out@N[rpi]==outBbind@N[RPI] &    out@S[rpi]==outBbind@S[RPI])) stop("matching for suppority part not correct, try creating boundary using built in functions, 
             or maybe need to rewrite the function analyzeBoundNBF, email package maintainer")

    # only change plower for binding='lower' and pupper for binding='upper'
    if (object@binding=="lower"){
        # binding=lower means stop at lower bound for sure
        # stopping at lower bound means plower<alpha
        # first find max of superiority or end using futility 
        pmaxSE<-max(out@plower[rpi])
        out@plower[rpi]<-outBbind@plower[RPI]
        # now find max of superiority or end not using futility 
        pmaxSEnbf<-max(out@plower[rpi])
        # for the p-values already on the futility boundary
        # adjust them so that they still make sense with the ordering
        # i.e., the p-value is larger on the futility boundary than
        # at other places
        plowerF<-out@plower[!rpi]    
        out@plower[!rpi]<- pmaxSEnbf +  (plowerF-pmaxSE)*( (1-pmaxSEnbf)/(1-pmaxSE) ) 
        # change upper Conf Limit at upper and end of bounary 
        # so the p-values will match the CI at those parts of the boundary
        if (cipMatch) out@upper[rpi]<-outBbind@upper[RPI]
    } else  if (object@binding=="upper"){
        # first find max of superiority or end using futility 
        pmaxSE<-max(out@pupper[rpi])
        out@pupper[rpi]<-outBbind@pupper[RPI]
        # now find max of superiority or end not using futility 
        pmaxSEnbf<-max(out@pupper[rpi])
        # for the p-values already on the futility boundary
        # adjust them so that they still make sense with the ordering
        # i.e., the p-value is larger on the futility boundary than
        # at other places
        pupperF<-out@pupper[!rpi]    
        out@pupper[!rpi]<- pmaxSEnbf +  (pupperF-pmaxSE)*( (1-pmaxSEnbf)/(1-pmaxSE) ) 
        # change upper Conf Limit at upper and end of boundary
        # so the p-values will match the CI at those parts of the boundary
        if (cipMatch) out@lower[rpi]<-outBbind@lower[RPI]
    } else stop("binding should be 'lower' or 'upper' ")
    if (alternative=="two.sided"){
        out@pval<- pmin(2*out@plower,2*out@pupper,1)
    } else if (alternative=="less"){
        out@pval<-out@plower
    } else if (alternative=="greater"){
        out@pval<-out@pupper
    } 
    if (stats=="pval"){
        # only output vector of pvalues not whole object 
        out<-out@pval
    }
    out
}

setGeneric("analyze",analyzeBound)

setMethod("analyze",
    signature(object = "bound"),
    analyzeBound)
setMethod("analyze",
    signature(object = "boundNBF"),
    analyzeBoundNBF)

setMethod("analyze",
    signature(object = "abparms"),
    function (object, theta0 = 0.5, stats = "all", alternative = "two.sided", 
        conf.level = 0.95, tsalpha = NULL,cipMatch=TRUE,...) 
    {
        B<-abtoBound(object)
        analyze(B, theta0=theta0, stats=stats, 
             alternative=alternative, conf.level=conf.level,
             tsalpha=tsalpha, cipMatch=cipMatch,...)

    }
)





designOBF<-function(Nmax, theta0=.5, k=Inf, tsalpha=NULL, alternative="two.sided", 
    conf.level=0.95,binding="both"){
    if (k<Nmax){
        # find k stopping times, always round up so last is always Nmax
        N<-ceiling( Nmax*(1:k)/k )
    } else {
        N<- 1:Nmax
    }  
    TSalpha<-getTSalpha(tsalpha,alternative,conf.level)
    alternative<- getAlternative(TSalpha)
    # first check that Nmax is large enough to reject without any early looks
    pupper<-dbinom(Nmax,Nmax,theta0)
    plower<-dbinom(0,Nmax,theta0)
    noUpperBound<- pupper>=TSalpha["alphaLower"]
    noLowerBound<- plower>=TSalpha["alphaUpper"]
    if (noUpperBound & noLowerBound) stop("Nmax too small to reject even without any early looks")
    # for two-sided boundaries, it may be that even with N=Nmax you will never 
    # stop early for one side
    #############################################################################
    # For alternative="greater" (alphaUpper=0) or "less" (alphaLower=0) 
    # then just get one-sided OF boundary
    #
    # for alternative="two.sided" we need to get two-sided boundary
    #   - we use the one-sided boundaries (elements of TSalpha) as the starting values
    #   - since we use the one-sided boundaries as starting values, we 
    #     know that the upper error (amount of alpha spent on the upper)
    #     when using the lower boundary must be less than without using 
    #     the lower boundary, so the one-sided boundary is the closest in 
    #     boundary (smallest Cupper possible for that one-sided alpha). 
    #     Similarly the one-sided lower boundary is the largest Clower possible.
    #  
    ############################################################################
    # to get the tolerance needed for Cupper and Clower, 
    # find out the maximum difference between values
    # and multiply by .9
    # round so that computer calculation differences (not real differences)
    # do not show up 
    ddd<-c((1:Nmax*theta0) - floor( 1:Nmax*theta0 ),
           ceiling( 1:Nmax*theta0 )- (1:Nmax*theta0)) 
    diffs<- (sort(unique( round(ddd,9))))
    if (any(diffs>0)){
        TOL<-max(10^-8,.9*min(abs(diffs[diffs>0])))
    } else {
        # if all diffs are 0 or 1, then use TOL=.9 to avoid computer ties
        TOL<- .9
    } 
    # can be problems with machine level error and then taking ceiling 
    # so that ceiling(12) becomes ceiling(12.0000000000001)=13
    # so need to round first
    rdigit<- ceiling(-log10(.Machine$double.eps^.5))
    # if alternative="two.sided" and noUpperBound, there is no need to find the upper 
    # bound because we know there will not be one
    if (alternative=="greater" | (alternative=="two.sided" & !noUpperBound)){
    # range for Cupper: 0 < Cupper < Nmax*(1-theta0)
         # getb:get boundary for a specific cutoff value, on the upper side
         getb<-function(Cupper){
             b<-ceiling(round(N*theta0 + Cupper,rdigit))
             ### first remove  b>N
             b[b>N]<-NA
             ### now remove values that cannot get to, i.e., 
             ### make boundary minimal
             ### need to find values where N[j]-b[j] = N[j+1]-b[j+1] 
             ### and set b[j+1] to NA
             diffb<-c(NA,diff(N-b))
             b[!is.na(diffb) & diffb==0]<-NA
             # no need to keep values with NA for both a and b
             keep<-  !is.na(b)
             a<-rep(NA,length(b[keep]))
             class(a)<-"numeric"
             ab<-new("abparms",Nk=N[keep],a=a,b=b[keep],binding=binding)
             ab
          }
         rootfuncCu<-function(Cu){
             ab<-getb(Cu )
             B<-abtoBound(ab)
             pval<-pCalc(B@S,B@N,B@K,B@order,theta0=theta0,alternative="greater",ponly=TRUE)
             # we want the largest S at Nmax, say (Smax,Nmax) to reject, 
             # but we do not want (Smax-1,Nmax) to reject
             # so use unirootDiscrete so it just barely rejects at (Smax,Nmax)
             nmax<-max(B@N)
             smax<-max(B@S[B@N==nmax])
             pval[B@S==smax & B@N==nmax] - TSalpha["alphaLower"]
           }
          # we want the rootfuncCu to give negatives values, i.e., pval[LastReject]<alpha
          # so use pos.side=FALSE
          #     to debug: print.steps=TRUE
          CupperOneSided<- unirootDiscrete(rootfuncCu,
                lower=0,
                upper=Nmax*(1-theta0),
                tol=TOL,
                pos.side=FALSE,print.steps=FALSE)$root
          if (alternative=="greater" | noLowerBound){
              # one more time, get the final boundary 
              ab<-getb(CupperOneSided)
              B<-abtoBound(ab)
          }
    } 
    if (alternative=="less" | (alternative=="two.sided" & !noLowerBound)){
     # range for Clower: -Nmax*theta0 < Clower < 0
        geta<-function(Clower){
            a<-floor(round(N*theta0 + Clower,rdigit))
            ### first remove a<0 
            a[a<0]<-NA
            ### now remove values that cannot get to, i.e., make boundary minimal
            # need to find values where a[j]=a[j+1] and set a[j+1] to NA
            diffa<-c(NA,diff(a))
            a[!is.na(diffa) & diffa==0]<-NA
            # no need to keep values with NA for both a and b
            keep<-  !is.na(a) 
            b<-rep(NA,length(a[keep]))
            class(b)<-"numeric"
            ab<-new("abparms",Nk=N[keep],a=a[keep],b=b,binding=binding)
            ab
         }
         rootfuncCl<-function(Cl){
             ab<-geta(Cl)
             B<-abtoBound(ab)
             pval<-pCalc(B@S,B@N,B@K,B@order,theta0=theta0,alternative="less",ponly=TRUE)
             # we want the smallest S at Nmax, say (Smin,Nmax) to reject, 
             # but we do not want (Smin+1,Nmax) to reject
             # so use unirootDiscrete so it just barely rejects at (Smin,Nmax)
             nmax<-max(B@N)
             smin<-min(B@S[B@N==nmax])
             pval[B@S==smin & B@N==nmax] - TSalpha["alphaUpper"]
           }
          # we want the rootfuncCl to give negatives values, i.e., pval[LastReject]<alpha
          # so use pos.side=FALSE
          #     to debug: print.steps=TRUE
          ClowerOneSided<- unirootDiscrete(rootfuncCl,
                lower=-Nmax*(theta0),
                upper=0,
                tol=TOL,
                pos.side=FALSE,print.steps=FALSE)$root
        if (alternative=="less" | noUpperBound){
            # one more time, get the final boundary 
            ab<-geta(ClowerOneSided)
            B<-abtoBound(ab)
        }
    } 
    if (alternative=="two.sided" & (!noUpperBound & !noLowerBound)){
        ##############################################################
        # For alternative=two.side, we already found CupperOneSided and 
        #   ClowerOneSided   that are the extremes for each side.
        #   Let the two.sided Cupper be CupperTwoSided and similarly 
        #   ClowerTwoSided
        #     We know that
        #          CupperOneSided <= CupperTwoSided <=  Nmax*(1-theta0)
        #     and    -Nmax*theta0 <= ClowerTwoSided <= ClowerOneSided
        #############################################################
        getab<-function(Clower,Cupper){
            a<-floor(round(N*theta0 + Clower,rdigit))
            b<-ceiling(round(N*theta0 + Cupper,rdigit))
            ### first remove a<0 and b>N
            a[a<0]<-NA
            b[b>N]<-NA
            ### now remove values that cannot get to, i.e., make boundary minimal
            # need to find values where a[j]=a[j+1] and set a[j+1] to NA
            diffa<-c(NA,diff(a))
            a[!is.na(diffa) & diffa==0]<-NA
            # need to find values where N[j]-b[j] = N[j+1]-b[j+1] and set b[j+1] to NA
            diffb<-c(NA,diff(N-b))
            b[!is.na(diffb) & diffb==0]<-NA
            # no need to keep values with NA for both a and b
            keep<-  !(is.na(a) & is.na(b))
            ab<-new("abparms",Nk=N[keep],a=a[keep],b=b[keep],binding=binding)
            ab
         }
         getUpperGivenLower<-function(clower){
         rootfuncUpper<-function(Cupper,Clower=ClowerOneSided){
             ab<-getab(Clower,Cupper)
             B<-abtoBound(ab)
             # for the upper side, use alternative="greater"
             pval<-pCalc(B@S,B@N,B@K,B@order,theta0=theta0,alternative="greater",ponly=TRUE)
             # we want the largest S at Nmax, say (Smax,Nmax) to reject, 
             # But we do not want (Smax-1,Nmax) to reject
             # so use unirootDiscrete so it just barely rejects at (Smax,Nmax)
             #
             # don't worry about the lower boundary for now
             nmax<-max(B@N)
             smax<-max(B@S[B@N==nmax])
             pval[B@S==smax & B@N==nmax] - TSalpha["alphaLower"]
          }
          # we want the rootfunc to give negatives values, i.e., pval[LastReject]<alpha
          # 
          # if CupperOneSided is already negative then no need to make it bigger
          if (rootfuncUpper(CupperOneSided)<0){
              CupperTwoSided<-CupperOneSided
          } else {
              CupperTwoSided<- unirootDiscrete(rootfuncUpper,
                lower=CupperOneSided,
                upper=Nmax*(1-theta0),
                tol=TOL,
                pos.side=FALSE,print.steps=FALSE)$root
          }
          CupperTwoSided 
          }
          getLowerGivenUpper<-function(cupper){
          rootfuncLower<-function(Clower, Cupper=cupper){
                  ab<-getab(Clower,Cupper)
                  B<-abtoBound(ab)
                  # for lower side use alternative=less
                  pval<-pCalc(B@S,B@N,B@K,B@order,theta0=theta0,alternative="less",ponly=TRUE)
                  # we want the smallest S at Nmax, say (Smin,Nmax) to reject, 
                  # but we do not want (Smin+1,Nmax) to reject
                  # so use unirootDiscrete so it just barely rejects at (Smin,Nmax)
                  nmax<-max(B@N)
                  smin<-min(B@S[B@N==nmax])
                  pval[B@S==smin & B@N==nmax] - TSalpha["alphaUpper"]
          }
          if (rootfuncLower(ClowerOneSided)<0){
              ClowerTwoSided<-ClowerOneSided
          } else {
              ClowerTwoSided<- unirootDiscrete(rootfuncLower,
                    lower=-Nmax*theta0,
                    upper=ClowerOneSided,
                    tol=TOL,
                    pos.side=FALSE,print.steps=FALSE)$root
          }
          ClowerTwoSided
          }
          # Now iterate
          niter<-0
          stopiter<-FALSE
          maxiter<-8
          ClowerTwoSided<-ClowerOneSided
          CupperTwoSided<-CupperOneSided
          while (niter<maxiter & !stopiter){
              cupperStart<-CupperTwoSided
              clowerStart<-ClowerTwoSided
              CupperTwoSided<-getUpperGivenLower(clowerStart)
              ClowerTwoSided<-getLowerGivenUpper(cupperStart)
              niter<-niter+1
              stopiter<-(identical(cupperStart,CupperTwoSided) && 
                         identical(clowerStart,ClowerTwoSided))
          }


          # one more time, get the final boundary 
          ab<-getab(ClowerTwoSided,CupperTwoSided)
          B<-abtoBound(ab)
    } 
    Best<-analyze(B,theta0=theta0,stats="all",tsalpha=TSalpha)
    Best
}

designAb<-function(Nk,a=NULL,b=NULL,theta0=NULL,tsalpha=NULL, alternative="two.sided", 
    conf.level=0.95,binding="both"){

    if (is.null(theta0)) stop("no value given for theta0")
    n<-length(Nk)
    ## start out with a and b with length 1 less than Nk
    if (length(a)>0 & length(a)+1!=n){
       stop("if a is not NULL it should be numeric of length one less than Nk")
    } else if (is.null(a)){
       a<-as.numeric(rep(NA,n))
    } else if (length(a)+1==n){
       a<-as.numeric(c(a,NA))
    }
    if (length(b)>0 & length(b)+1!=n){
       stop("if b is not NULL it should be numeric of length one less than Nk")
    } else if (is.null(b)){
       b<-as.numeric(rep(NA,n))
    } else if (length(b)+1==n){
       b<-as.numeric(c(b,NA))
    }


    ab<-new("abparms",Nk=Nk,a=a,b=b,binding=binding)
    #b<-as(ab,"bound")
    b<-abtoBound(ab)
    B<-analyze(b, theta0=theta0, tsalpha=getTSalpha(tsalpha,alternative,conf.level))
    B
}

designSimon<-function(theta0,theta1,alpha=.05,beta=.2,type=c("optimal","minimax")){
    # need clinfun R package
    # for notation and description 
    # see ?ph2simon after loading clinfun
    p0<-theta0
    p1<-theta1
    # MADE package depend on clinfun
    # so do not need the following commented lines
    #if (!require(clinfun)) stop("need to install 
    # the clinfun R package first")
    type<-match.arg(type)
    x<-ph2simon(p0,p1,alpha,beta)
    # see print.ph2simon
    xout <- x$out
    nmax <- x$nmax
    n <- nrow(xout)
    if (type=="optimal"){
        index <- ((1:n)[xout[, 5] == min(xout[, 5])])[1]
    } else index<-1
    out<-xout[index, ]
    B<-designAb(Nk=c(out[c("n1","n")]),a=c(out[c("r1")]),theta0=p0,
         conf.level=1-alpha,alternative="greater")
    B
}



#Bnew<-designAb(Nk=c(13,43),a=c(3,12),theta0=.2,conf.level=.95,alternative="greater")


##############################################################
#
#   METHODS
#
##############################################################
plotBsb<-function(x,rcol=c(orange="#E69F00",blue="#56B4E9",green="#009E73"), 
    rpch=c(openCircle=1,filledCircle=16,filledDiamond=18), 
    bplottype="NS",newplot=TRUE, dtext=NULL, grid=50, xlab=NULL,ylab=NULL,...){
    # default colors for rcol: 
    #              orange (fail to reject), put
    #              sky blue (reject theta1>theta0), 
    #              greenish blue (reject theta1<theta0)
    #  allow good differentiation for those with several of the 
    # most common type of colorblindness
    rcolNames<-names(rcol)
    rpchNames<-names(rpch)
    if (!is.logical(newplot)) stop("newplot should be logical")
    if (is.null(dtext)){
        # default dtext: if only 1 panel put text if newplot, otherwise do not
        if (all.equal(par()$mfrow,c(1,1))==TRUE){ 
            dtext<-newplot
        } else dtext<-FALSE
    }
    tsalpha<-x@tsalpha
    M<-max(x@N) + 2
    if (length(rpch)==1) rpch<-rep(rpch,3)
    if (length(rcol)==1) rcol<-rep(rcol,3)

    # the terminology is confusing
    # alphaLower is the error allowed on the lower side
    # so if alternative="less" then you want to reject if theta 
    # is smaller and you do not care about error on the lower side, 
    # so you set alphaLower=0, all the erro is in gt
    # also if alternative="less" more extreme values will be less
    # so you want to use the one-sided p-value, plower
    # but you reject when plower<gt
    RejectUpper<- x@pupper< tsalpha["alphaLower"]
    RejectLower<- x@plower< tsalpha["alphaUpper"]
    FailToReject<- !RejectUpper & !RejectLower


    if (bplottype=="NS"){
        if (is.null(xlab)){ 
            XLAB<- "Total (N)"
        } else { XLAB<-xlab }
        if (is.null(ylab)){ 
            YLAB<- "Sucesses (S)"
        } else { YLAB<-ylab }

        if (newplot){
            plot.default(x = 0:M, y = 0:M, xlab = XLAB, 
                ylab = YLAB, type = "n", 
                xaxs="i",yaxs="i",axes = FALSE,...)
            axis(1)
            axis(2,las=1)
            minX<- -10
            maxX<- M+10
            # gray out upper triangle
            polygon(c(minX,maxX,minX,minX),
              c(minX,maxX,maxX,minX),
              col=gray(.9),border=NA)
            if (max(x@N)<=grid){
                abline(h = 0, v= 0)
                segments(0:(M-1),rep(0, M), rep(M,M), M:1, lty = 3, col = "grey")
	          segments(1:M,1:M, rep(M,M),1:M, lty = 3, col = "grey")
             }
            # plot null hypothesis line
            lines(c(0,2*max(x@N)), c(0,2*max(x@N)*x@theta0), lwd=3, col="gray")
        }
	  points(x@N[FailToReject],x@S[FailToReject], pch=rpch[1], col = rcol[1])

        if (any(RejectUpper)) points(x@N[RejectUpper],x@S[RejectUpper], pch=rpch[2], col = rcol[2])
        if (any(RejectLower)) points(x@N[RejectLower],x@S[RejectLower], pch=rpch[3], col = rcol[3])
        if (dtext){            
            text(.1*max(x@N),.95*max(x@N),labels=c("Gray line is null hyothesis,"),pos=4)
            text(.1*max(x@N),.9*max(x@N),labels=bquote(theta[0] == .(x@theta0)),pos=4)
            if (is.null(rcolNames) & is.null(rpchNames)){
                text(.1*max(x@N),.85*max(x@N),labels=paste("Reject (top) if p[upper]<",x@tsalpha["alphaLower"],sep=""),pos=4)
                text(.1*max(x@N),.8*max(x@N),labels=paste("Reject (bottom) if p[lower]<",x@tsalpha["alphaUpper"],sep=""),pos=4)
            } else {
                text(.1*max(x@N),.85*max(x@N),labels=paste("Reject (",rcolNames[2]," ",rpchNames[2],") if p[upper]<",x@tsalpha["alphaLower"],sep=""),pos=4)
                text(.1*max(x@N),.8*max(x@N),labels=paste("Reject (",rcolNames[3]," ",rpchNames[3],") if p[lower]<",x@tsalpha["alphaUpper"],sep=""),pos=4)
            }
            if (x@binding=="both"){ btext<-" boundaries"
            } else btext<-" boundary"
            text(.1*max(x@N),.75*max(x@N),labels=paste("binding on ",x@binding,btext,sep=""),pos=4)
        }
    } else if (bplottype=="FS"){
        if (is.null(xlab)){ 
            XLAB<- "Failures (F)"
        } else { XLAB<-xlab }
        if (is.null(ylab)){ 
            YLAB<- "Sucesses (S)"
        } else { YLAB<-ylab }


        if (newplot){
            plot.default(x = 0:M, y = 0:M, xlab = XLAB, 
                ylab = YLAB, type = "n", 
                xaxs="i",yaxs="i",axes = FALSE,...)
            axis(1)
            axis(2,las=1)
            maxX<- max(x@N)
            
            # gray out upper triangle
            polygon(c(maxX,0,0,2*maxX,2*maxX,maxX),
              c(0,maxX,2*maxX,2*maxX,0,0),
              col=gray(.9),border=NA)

            if (max(x@N)<=grid){
                abline(h = 0, v= 0)
                segments(1:M,rep(0, M), 1:M, rep(2*M,M), lty = 3, col = gray(.9))
	          segments(rep(0,M),1:M,rep(2*M,M),1:M, lty = 3, col = gray(.9))
             }
            # plot null hypothesis line
            lines(c(0,2*max(x@N)*(1-x@theta0)), c(0,2*max(x@N)*x@theta0), lwd=3, col=gray(.9))
        }
	  points(x@N[FailToReject]-x@S[FailToReject],x@S[FailToReject], pch=rpch[1], col = rcol[1])

        if (any(RejectUpper)) points(x@N[RejectUpper]-x@S[RejectUpper],x@S[RejectUpper], pch=rpch[2], col = rcol[2])
        if (any(RejectLower)) points(x@N[RejectLower]-x@S[RejectLower],x@S[RejectLower], pch=rpch[3], col = rcol[3])
        if (dtext){
            text(.5*max(x@N),.95*max(x@N),labels=c("Gray line is null hyothesis,"),pos=4)
            text(.5*max(x@N),.9*max(x@N),labels=bquote(theta[0] == .(x@theta0)),pos=4)
            if (is.null(rcolNames) & is.null(rpchNames)){
                text(.5*max(x@N),.85*max(x@N),labels=paste("Reject (top) if p[upper]<",x@tsalpha["alphaLower"],sep=""),pos=4)
                text(.5*max(x@N),.8*max(x@N),labels=paste("Reject (bottom) if p[lower]<",x@tsalpha["alphaUpper"],sep=""),pos=4)
            } else {
                text(.5*max(x@N),.85*max(x@N),labels=paste("Reject (",rcolNames[2]," ",rpchNames[2],") if p[upper]<",x@tsalpha["alphaLower"],sep=""),pos=4)
                text(.5*max(x@N),.8*max(x@N),labels=paste("Reject (",rcolNames[3]," ",rpchNames[3],") if p[lower]<",x@tsalpha["alphaUpper"],sep=""),pos=4)
            }
            if (x@binding=="both"){ btext<-" boundaries"
            } else btext<-" boundary"
            text(.5*max(x@N),.75*max(x@N),labels=paste("binding on ",x@binding,btext,sep=""),pos=4)
        }
    } else if (bplottype=="NE"){
        if (is.null(xlab)){ 
            XLAB<- "Total (N)"
        } else { XLAB<-xlab }
        if (is.null(ylab)){ 
            YLAB<- "Probability of Sucess"
        } else { YLAB<-ylab }


        # plot N by estimates (and CI)
        # only plot estimate next to the continuation region 
        #
        # if more than 1 upper or more than 1 lower per N (before end) stop
        nup<-sort(x@N[x@UL=="upper"])
        nlo<-sort(x@N[x@UL=="lower"])
        if (any(diff(nup)==0) | any(diff(nlo)==0)) stop("bplotype='NE' does not work when more than one point per N on upper or lower boundary")

        uN<-sort(unique(x@N))
        nN<-length(uN)
        n<-length(x@N)
        keep<-rep("no",n)
        if (nN>1){
          for (i in 1:(nN-1)){
            IL<- x@N==uN[i] & x@UL=="lower"
            if (any(IL)){
                iL<- IL & x@S==max(x@S[IL])
                keep[iL]<-"lower"
            }
            IU<- x@N==uN[i] & x@UL=="upper"
            if (any(IU)){
                iU<- IU & x@S==min(x@S[IU])
                keep[iU]<-"upper"
            }
          }    
        } 
        estimate<-x@estimate
        lower<-x@lower
        upper<-x@upper
        N<-x@N
        I<- keep!="no"
        #YLIM<-range(c(lower[I],upper[I]))
        YLIM<-c(0,1)
        if (newplot){
            plot(c(0,max(N[I])),YLIM,type="n",xlab=XLAB,ylab=YLAB,
                xaxs="i",yaxs="i",axes=FALSE,...)
            axis(1)
            axis(2)
        }
        drawEstCI<-function(I,COL,LWD,PCH){
            points(N[I],estimate[I],col=COL,pch=PCH)
            segments(N[I],lower[I],N[I],upper[I],col=COL,lwd=LWD)
        }
        drawEstCI(keep=="upper",rcol[2],3,PCH=rpch[2])
        drawEstCI(keep=="lower",rcol[3],1,PCH=rpch[3])
    } else if (bplottype=="NZ" | bplottype=="NB"){
        Z<- (x@S/x@N - x@theta0)/sqrt( (x@theta0*(1-x@theta0))/x@N)
        # for graying out impossible areas, calculate Zhi and Zlo for 0:M
        Zhi<-(1- x@theta0)/sqrt( (x@theta0*(1-x@theta0))/0:M)
        Zlo<-(0- x@theta0)/sqrt( (x@theta0*(1-x@theta0))/0:M)
        if (is.null(xlab)){ 
            XLAB<- "Total (N)"
        } else { XLAB<-xlab }

        if (bplottype=="NZ"){
            # Y=Z score 
            Y<-Z
            Yhi<-Zhi
            Ylo<-Zlo
            if (is.null(ylab)){ 
                YLAB<-"Z-score"
            } else { YLAB<-ylab }
        } else if (bplottype=="NB"){
            # if 'NB'  change Y=Z-score to Y=B-value
            Y<- sqrt(x@N/max(x@N))*Z
            Yhi<-sqrt((0:M)/max(x@N))*Zhi
            Ylo<-sqrt((0:M)/max(x@N))*Zlo
            if (is.null(ylab)){ 
                YLAB<-"B-value"
            } else { YLAB<-ylab }
        }
        if (newplot){
            plot.default(x =c(0,M), y = 1.2*range(Y), xlab = XLAB, 
                ylab = YLAB, type = "n", 
                xaxs="i",yaxs="i",axes = FALSE,...)
            axis(1)
            axis(2,las=1)
            box()
            # gray out impossible values 
            polygon(c(0,0:M,M,0,0),
              c(0,Yhi,max(Yhi),max(Yhi),0),
              col=gray(.9),border=NA)
           polygon(c(0,0:M,M,0,0),
              c(0,Ylo,min(Ylo),min(Ylo),0),
              col=gray(.9),border=NA)


            # plot null hypothesis line
            lines(c(0,2*max(x@N)), c(0,0), lwd=3, col="gray")
        }
	  points(x@N[FailToReject],Y[FailToReject], pch=rpch[1], col = rcol[1])

        if (any(RejectUpper)) points(x@N[RejectUpper],Y[RejectUpper], pch=rpch[2], col = rcol[2])
        if (any(RejectLower)) points(x@N[RejectLower],Y[RejectLower], pch=rpch[3], col = rcol[3])
        #if (dtext){
        #    text(.1*max(x@N),.95*max(x@N),labels=c("Gray line is null hyothesis,"),pos=4)
        #    text(.1*max(x@N),.9*max(x@N),labels=bquote(theta[0] == .(x@theta0)),pos=4)
        #    text(.1*max(x@N),.85*max(x@N),labels=paste("Reject (",rcolNames[2],") if p[upper]<",x@tsalpha["alphaLower"],sep=""),pos=4)
        #    text(.1*max(x@N),.8*max(x@N),labels=paste("Reject (",rcolNames[3],") if p[lower]<",x@tsalpha["alphaUpper"],sep=""),pos=4)
        #    if (x@binding=="both"){ btext<-" boundaries"
        #    } else btext<-" boundary"
        #    text(.1*max(x@N),.75*max(x@N),labels=paste("binding on ",x@binding,btext,sep=""),pos=4)
        #}
    }
}

pointsBsb<-function(x,...){
    plotBsb(x,newplot=FALSE,...)
}

#define a plot method for the boundary class
setGeneric("plot")
setMethod("plot", signature(x="boundEst",y="missing"),plotBsb)

setGeneric("points")
setMethod("points", 
     signature(x="boundEst"),pointsBsb)

#pointsAbparms<-function(x,...){
#    plotAbparms(x,...,newplot=FALSE)
#}
#setMethod(f = "plot", signature(x="abparms",y="missing"),
#    definition = plotAbparms)

#setMethod(f = "points", signature(x="abparms",y="missing"),
#    definition = pointsAbparms)



missNAbparms<-function(ab,missN=NULL,...){
    if (!is.null(missN)){
        uMN<-sort(unique(missN))
        nMN<-length(uMN)
        Nk<-c(ab@Nk)
        a<-c(ab@a)
        b<-c(ab@b)
   
        if (any(uMN>max(Nk))) stop("missN at Nk value greater than Nmax")
        
          
        # create a function for missing only 1 Nk value
        # input abparms class value as a list
        # instead of as class abparms
        # to keep from doing extra abparms validation checks 
        # each iteration

        ablist<-list(a=a,b=b,Nk=Nk)

        modifyi<-function(ablist,missi){
            a<-ablist$a
            b<-ablist$b
            Nk<-ablist$Nk
            n<-length(Nk)
            keep<-rep(TRUE,n)
            I<- Nk==missi
            Nmax<-max(Nk)
            # if not planned to stop at missi, then no need to do it
            # do not call this program if not planned stop, give error if somehow call program
            if (!any(I)) stop("rewrite program, should not call modifyi if missi not in Nk")
            # if missi=Nk[I] then add 1 to get Nnew
            Nnew<- Nk[I]+1
            # for anew and bnew go forward 1 step, 
            # keeping bound closed
            # but not stopping more often then if did not miss
            # that means, anew=a and bnew=b+1
            anew<- a[I]
            bnew<- b[I]+1
            keep[I]<-FALSE
            addN<-adda<-addb<-NA
            if (Nnew==Nmax){
                # if Nnew=Nmax then do not need to do anything
                # since we stop for everything at Nmax anyway
            } else if (any(Nk==Nnew)){
                # already a stop at Nnew, so need to find if anew and bnew 
                # are part of the continuation region and should replace the old 
                # values, or are part of the stopping time already 
                iNew<-(1:n)[Nk==Nnew]
                # if old a=a[iNew] is NA, then replace with anew
                #   otherwise replace with max of both
                # if all NA then do not need to do anything
                if ( !all(is.na(c(anew,a[iNew]))) ) a[iNew]<- max(anew,a[iNew],na.rm=TRUE)
                # if old b=b[iNew] is NA, then replace with bnew
                #   otherwise replace with min of both
                # if all NA then do not need to do anything
                if ( !all(is.na(c(bnew,b[iNew]))) ) b[iNew]<- min(bnew,b[iNew],na.rm=TRUE)
            }  else if (!any(Nk==Nnew)){
                if (Nnew>Nmax){ Nmax<-Nnew
                } else { 
                    addN<- Nnew
                    adda<- anew
                    addb<- bnew
                }
            }         
            # remember adda and addb can be NA, so keep 
            # whenever !is.na(addN)
            aa<-c(a[keep],adda[!is.na(addN)])
            bb<-c(b[keep],addb[!is.na(addN)])
            NN<-c(Nk[keep],addN[!is.na(addN)]) 
            o<-order(NN)
            out<-list(a=aa[o],b=bb[o],Nk=NN[o])
        }
        # modify one at a time
        # need to do it this way in case the ith modify creates 
        # and N value that was not there before the modify
        for (i in 1:nMN){
           if (any(ablist$Nk %in% uMN[i])){
               ablist<-modifyi(ablist,missi=uMN[i])
           } else {
               warning(paste("missN=",uMN[i]," but no planned stop at that Nk"))
           }
        }

        out<-new("abparms",a=ablist$a,b=ablist$b,Nk=ablist$Nk,binding=ab@binding)
    } else {
        ## missN=NULL so return what you started with
        out<-ab
    }
    out
}
earlyCloseoutAbparms<-function(ab,earlyCloseout){
    Nk<-ab@Nk
    k<-length(Nk)
    pick<-(1:k)[Nk<=earlyCloseout]
    i<-length(pick)
    a<-ab@a[pick]
    # a and b at end are NA by definition
    a[i]<-NA
    b<-ab@b[pick]
    b[i]<-NA
    newab<-new("abparms",a=a,b=b,Nk=Nk[pick],binding=ab@binding)
    newab
}
#abnew<-missNAbparms(ab,missN=13:18)
#abnew
modify<-function(b,missN=NULL,theta0=NULL,tsalpha=NULL,conf.level=NULL,
    alternative=NULL,cipMatch=TRUE, binding=NULL, closeout=NULL,...){
    # Not the most efficient programming, since recalculates CI and p-values
    # even if they do not need to be recalculated. 
    # So there is room for efficiency gains.
    if (is.null(tsalpha) & is.null(conf.level) & is.null(alternative)){
        TSalpha<-b@tsalpha
    } else {
        TSalpha<-getTSalpha(tsalpha,alternative,conf.level)
    }
    if (is.null(theta0)){
        Theta0<-b@theta0
    } else {
        Theta0<-theta0
    }
    if (!is.null(missN) & !is.null(closeout)){
        if (max(missN)>=closeout) stop("missN values specified at or after closeout")
    }
    if (!is.null(missN)){
        ab<-missNAbparms(b,missN=missN)
        b<-abtoBound(ab)
    } 
    if (!is.null(closeout)){
        ab<-earlyCloseoutAbparms(b,earlyCloseout=closeout)
        b<-abtoBound(ab)
    }
    if (!is.null(binding)){
         b@binding<-binding
         # abtoBound can take 'bound' objects because they contain 'abparms' objects
         b<-abtoBound(b)
    }
    Best<-analyze(b,theta0=Theta0,stats="all",tsalpha=TSalpha, cipMatch=cipMatch)
    Best
}

#Bmod<-modify(B,conf.level=.95,alternative="two.sided")

powerBsb<-function(object,theta=.6,alternative=NULL){
    B<-object
    # if alternative=NULL, use alternative in object
    if (is.null(alternative)){
        alternative<-getAlternative(B@tsalpha)
    }

    poweri<-function(thetai,R=RejectUpper){
        sum(B@K[R]*(thetai)^B@S[R]*(1-thetai)^(B@N[R]-B@S[R]))
    }

    ntheta<-length(theta)
    pow<-rep(NA,ntheta)
        # the terminology is confusing
        # alphaLower is the error allowed on the lower side
        # so if alternative="less" then you want to reject if theta 
        # is smaller and you do not care about error on the lower side, 
        # so you set alphaLower=0, all the error is in alphaUpper
        # also if alternative="less" more extreme values will be less
        # so you want to use the one-sided p-value, plower
        # but you reject when plower<alphaUpper
    RejectUpper<-  B@pupper< B@tsalpha["alphaLower"]
    RejectLower<-  B@plower< B@tsalpha["alphaUpper"]

    if (alternative=="less"){
        for (i in 1:ntheta){
            pow[i]<-poweri(theta[i],R=RejectLower)
        }
    } else if (alternative=="greater"){
        for (i in 1:ntheta){
            pow[i]<-poweri(theta[i],R=RejectUpper)
        }
    } else if (alternative=="two.sided"){
        for (i in 1:ntheta){
            if (theta[i]<B@theta0){
                pow[i]<-poweri(theta[i],R=RejectLower)
            } else if (theta[i]>B@theta0){
                pow[i]<-poweri(theta[i],R=RejectUpper)
            } else if (theta[i]==B@theta0){
                if (B@tsalpha["alphaUpper"]<=B@tsalpha["alphaLower"]){
                pow[i]<-poweri(theta[i],R=RejectUpper)
                } else {
                pow[i]<-poweri(theta[i],R=RejectLower)
                }
            }
        }
    }
    pow
}

EN<-function(object,theta=.6){
    B<-object
    ENi<-function(thetai){
        sum(B@N*B@K*(thetai)^B@S*(1-thetai)^(B@N-B@S))
    }
    ntheta<-length(theta)
    ENout<-rep(NA,ntheta)
    for (i in 1:ntheta){
        ENout[i]<-ENi(theta[i])
    }
    ENout
}

prStop<-function(object,theta=NULL){
    if (is.null(theta)) theta<-object@theta0
    if (length(theta)>1) stop("length of theta cannot be bigger than 1")
    B<-object 
    iU<- B@UL=="upper"
    iL<- B@UL=="lower"

    Nupper<- sort(unique(B@N[iU]))
    Nlower<- sort(unique(B@N[iL]))
    Nend<- max(B@N)
    dStopUpper<-rep(NA,length(Nupper))
    for (i in 1:length(Nupper)){
        I<- B@N==Nupper[i] & B@UL=="upper"
        dStopUpper[i]<-  sum(B@K[I]*(theta)^B@S[I] * (1-theta)^(B@N[I]-B@S[I]))
    }
    dStopLower<-rep(NA,length(Nlower))
    for (i in 1:length(Nlower)){
        I<- B@N==Nlower[i] & B@UL=="lower"
        dStopLower[i]<-  sum(B@K[I]*(theta)^B@S[I] * (1-theta)^(B@N[I]-B@S[I]))
    }
    pStopUpper<-cumsum(dStopUpper)
    pStopLower<-cumsum(dStopLower)
    I<- B@N==Nend
    dStopEnd<- sum(B@K[I]*(theta)^B@S[I] * (1-theta)^(B@N[I]-B@S[I]))
    check<- sum(dStopUpper) + sum(dStopLower) + dStopEnd

    list(B=B,Nupper=Nupper,dStopUpper=dStopUpper,
        pStopUpper=pStopUpper, Nlower=Nlower, dStopLower=dStopLower,
        pStopLower=pStopLower,Nend=Nend,dStopEnd=dStopEnd,check=check)
}
#b<-designOBF(50,tsalpha=c(.05,.1))
#plotBsb(b)
#x<-prStop(b)
plotPrStop<-function(x,
    rcol=c(orange="#E69F00",blue="#56B4E9",green="#009E73"),
    rpch=c(openCircle=1,filledCircle=16,filledDiamond=18),total=FALSE){

   if (total){
       # plot total probability of stopping
       Ntotal<-sort(unique(c(x$Nlower,x$Nupper)))
       ptotal<-rep(NA,length(Ntotal))
       for (i in 1:length(Ntotal)){
           iu<- x$Nupper<=Ntotal[i]
           il<- x$Nlower<=Ntotal[i]
           ptotal[i]<- sum(x$dStopUpper[iu]) + sum(x$dStopLower[il])
       }
       YLIM<-range(c(0,ptotal,x$B@tsalpha))
   } else YLIM<- range(c(0,x$pStopUpper,x$pStopLower,x$B@tsalpha))
   plot(c(0,max(x$B@N)),YLIM,type="n",xlab="Total(N)",
     ylab="Probibility of Stopping by N")
   CEX<-c(.75,1,1.5)
   points(x$Nupper,x$pStopUpper,col=rcol[2],pch=rpch[2],cex=CEX[2])
   points(x$Nlower,x$pStopLower,col=rcol[3],pch=rpch[3],cex=CEX[3])
   if (total) points(Ntotal,ptotal,col=rcol[1],pch=rpch[1],cex=CEX[1])

   if (total){
       legend(0,YLIM[2],legend=c("Any stopping","upper","lower"),pch=rpch,col=rcol)
   } else {
       legend(0,YLIM[2],legend=c("upper","lower"),pch=rpch[-1],col=rcol[-1])
   }
}
#plotPrStop(x,total=TRUE)
stopTable<-function(object,Nrange=c(0,Inf),Srange=c(0,Inf),output="all",file="stopTableOutput.csv"){

    if (is.null(file)) file<- paste("tempfile.",paste(sample(letters,6,replace=TRUE),collapse=""),".csv",sep="")
    stab<-data.frame(
        N=object@N,
        S=object@S,
        estimate=object@estimate,
        lower=object@lower,
        upper=object@upper,
        pL=object@plower,
        pU=object@pupper,
        pts=pmin(1,2*object@plower,2*object@pupper),
        UL=object@UL)
        #rounded.pvalue=write.p(object@pval,digits[2]),


    N<-object@N
    S<-object@S

    if (length(Nrange)==1) Nrange<-c(Nrange,Nrange)
    if (length(Srange)==1) Srange<-c(Srange,Srange)

    # sort just in case order is backwards
    Nrange<-sort(Nrange)
    Srange<-sort(Srange)
    Nmin<-Nrange[1]
    Nmax<-Nrange[2]
    Smin<-Srange[1]
    Smax<-Srange[2]
    # pick out which rows to output
    i<-   (N>=Nmin & N<=Nmax & S>=Smin & S<=Smax)
    if (!any(i)) warning("no output: no stopping points meet selection criteria given by Nrange and/or Srange")

    SLIST<-list(conf.level=1-sum(object@tsalpha),
        lowerError=object@tsalpha[1],
        upperError=object@tsalpha[2],
        nullTheta=object@theta0,
        binding=object@binding,
        alternative=getAlternative(object@tsalpha),
        table=stab[i,])
    ## use S3 method, just need it for default printing 
    class(SLIST)<-"stopTable"
    STAB<-data.frame(stab[i,],
            conf.level=1-sum(object@tsalpha),
            lowerError=object@tsalpha[1],
            upperError=object@tsalpha[2],
            nullTheta=object@theta0,
            binding=object@binding,
            alternative=getAlternative(object@tsalpha)
            )   
    if (output=="csv"){   
        if (file!="") write.csv(STAB,file=file,row.names=FALSE)
    } else if (output=="data.frame"){
        out<-STAB
    } else {
        # outputs stopTable object unless output ='data.frame'
        out<-SLIST
    }
    out
}

print.stopTable<-function(x,digits=c(3,5),maxnprint=Inf,...){
         if (length(digits)==1) digits<-rep(digits,2)
    # write.p rounds and translates p-value vector to character vector, giving <.001 if needed
    write.p <-function(p,digits){
         p<-formatC(p,digits=digits,format='f')
         pout<-paste("",as.character(p),sep="")
         pout[as.numeric(p)==0]<- paste("<.",paste(rep("0",digits-1),collapse=""),"1",sep="")
         pout
    }
         cat(paste("H0: theta0=",x$nullTheta[1]),"\n")
         cat(paste("confidence level=",x$conf.level[1]),"\n")
         cat(paste("alternative=",x$alternative[1]),paste(": lowerError=",x$lowerError[1],
                 ", upperError=",x$upperError[1]),"\n")

         tempd<-data.frame(N=x$table$N,S=x$table$S,estimate=round(x$table$estimate,digits[1]),
                lower=round(x$table$lower,digits[1]),
                upper=round(x$table$upper,digits[1]),
                pLower=write.p(x$table$pL,digits[2]),
                pUpper=write.p(x$table$pU,digits[2]),
                pTwoSided=write.p(x$table$pts,digits[2]))


         cat(paste("Boundaries binding=",x$binding[1]),"\n")
         nrowprint<-nrow(tempd)
         if (nrowprint<maxnprint){ 
             print(tempd,row.names=FALSE)
         } else {
             print(tempd[1:maxnprint,],row.names=FALSE)
             cat(paste("...printed first",maxnprint,"rows only out of",nrowprint),"\n")
             cat("to see all rows, use maxnprint=Inf, or use ","\n")
             cat("   stopTable(...,output='data.frame')","\n")
             cat(" and look at that, or","\n")
             cat("see .csv file created using stopTable(...,output='csv')","\n")
         }
}

#stopTable(x,Nmin=40,Nmax=Inf)
#stopTable(x,S=10,N=43)

showBoundEst<-function(object){
   str(object)
   #stopTable(object,output="print")
   cat("Use stopTable to see confidence intervals and p-values","\n")
   cat("at each stopping point","\n")
}
#setGeneric("show")
setMethod("show", "boundEst",showBoundEst)

summaryBoundEst<-function(object){
   stopTable(object,output="print")
}
setGeneric("summary")
setMethod("summary", "boundEst",summaryBoundEst)




##################################################################
#   Fixed Sample Functions
#
#################################################################

powerUpper<-function(n,theta0,theta1,alpha){
    # reject if one-sided p-value<alpha
    # means reject if Cx>=qbinom(alpha,lower.tail=FALSE)
    Cx<- qbinom(alpha,n,theta0,lower.tail=FALSE)
    # power is Pr[reject | theta1]
    pbinom(Cx,n,theta1,lower.tail=FALSE)
}

powerLower<-function(n,theta0,theta1,alpha){
    # reject if one-sided p-value<alpha
    # means reject if Cx<=qbinom(alpha,lower.tail=TRUE)
    Cx<- qbinom(alpha,n,theta0,lower.tail=TRUE)
    # power is Pr[reject | theta1]
    pbinom(Cx-1,n,theta1,lower.tail=TRUE)
}


sampleSizeFixedBin<-function(theta0,theta1,alpha,power=.8,allNgreater=TRUE){
    if (theta1>theta0){ powerFunc<-powerUpper
    } else powerFunc<-powerLower

    # asymptotic normal theory ss
    ssAsy<- (qnorm(1-alpha) + qnorm(power))^2 * theta1*(1-theta1)/(theta1-theta0)^2
    ssAsy

    # get range
    ssMin<- max(2,floor(.5*ssAsy))
    # in case ssAsy<1 or something, make ssMax at least 20
    ssMax<- max(20,ceiling(2*ssAsy))
    ssRange<-ssMin:ssMax
    pow<-powerFunc(ssRange,theta0,theta1,alpha)
    if (allNgreater){
        # first pick all sample sizes with greater than nominal power
        ssGT<-ssRange[pow>power]
        diffs<-diff(c(0,ssGT))
        # usually power will cross back and forth across the nominal power 
        # as the sample size increases for a while, until all further
        # increased in sample size give pow greater than nominal power
        # since diffs[1]=ssGT[1] which is at least 2, we will always have 
        # a diffs>1 value 
        ss<-max(ssGT[diffs>1])
     } else {
    # just take smallest n that gives power greater than nominal
        ss<-min(ssRange[pow>power])
    }
    ss
}


designFixed<-function(Nmax, theta0=.5, tsalpha=NULL, alternative="two.sided", 
    conf.level=0.95){
    # binding not needed for Fixed designs, set to "both"
    B<-designAb(Nk=Nmax, theta0=theta0, tsalpha=tsalpha, alternative=alternative, conf.level=conf.level,
       binding="both")
    B
}

designFixedpower<-function(theta0=.5, theta1=.6, power=.8, maxNmax=Inf, tsalpha=NULL, alternative=NULL, 
    conf.level=0.95,allNgreater=FALSE){
    TSalpha<-getTSalpha(tsalpha,alternative,conf.level)
    if (theta1>theta0){
         onesided.alpha<-TSalpha[1]
         powerFunc<-powerUpper
    } else {
         onesided.alpha<-TSalpha[2]
         powerFunc<-powerLower
    }
    # if Nmax<Inf check that Nmax is large enough 
    Nmax<-maxNmax
    if (Nmax<Inf){
        powMax<-powerFunc(Nmax,theta0,theta1,onesided.alpha)
        if (powMax<power) stop("Nmax not large enough for specified power")
    }
    ss<-sampleSizeFixedBin(theta0,theta1,onesided.alpha,power=power,allNgreater=allNgreater)
    if (ss>Nmax) stop(paste("Sample Size of ",ss,"needed. Nmax is",Nmax))
    # binding does not mean anything for Fixed sample designs, so set to both
    B<-designAb(Nk=ss, theta0=theta0, tsalpha=TSalpha, binding="both")
    B
}
#designFixedpower(theta0=.5,theta1=.9)

#####################################################################################
#
#    Design using O-F method 
# 
#####################################################################################
designOBFpower<-function(
        theta0=.5,theta1=.6,
        k=Inf,
        power=.9,
        tsalpha=NULL,
        alternative="two.sided",
        conf.level=0.95,
        binding="both",
        allNgreater=FALSE,
        checkmax=10, maxNmax=2*ss){
    TSalpha<-getTSalpha(tsalpha,alternative,conf.level)
    alternative<-getAlternative(TSalpha)

    # want sample size from fixed design, that gives power=power
    # using a one-sided reject level.
    #  if theta1>theta0, then reject if thetahat is large
    #  so allow error on lower side, and onesided.alpha=TSalpha[1]
    if (theta1>theta0){
         onesided.alpha<-TSalpha[1]
         alt<-"greater"
    } else {
         onesided.alpha<-TSalpha[2]
         alt="less"
    }
    if (onesided.alpha==0){
        stop("tsalpha=0 on wrong side, e.g., tsalpha[1]=0 when theta1>theta0")
    }
 
    ## check that Nmax is large enough 
    ssFixed<-sampleSizeFixedBin(theta0,theta1,onesided.alpha,power,allNgreater)
    ss<-ssFixed 
    Nmax<-maxNmax
    pow<-0  
    if (Nmax<ss) stop("Nmax too small for this power")
    getPower<-function(ss){
        b<-designOBF(ss,k=k,theta0=theta0,tsalpha=TSalpha,binding=binding)
        # get one-sided p-value based on whether theta1> (or <) theta0, using alt
        pow<-powerBsb(b,theta=theta1)
        list(Bs=b,pow=pow)
    }

    iter<-0
    while (ss<=Nmax & pow<power){
        iter<-iter+1
        if (checkmax==iter){
            powmax<-getPower(Nmax)$pow
            if (powmax<power) stop("Nmax too small or alpha too large for this power")
        }
        powlist<-getPower(ss)
        pow<-powlist$pow
        Bs<-powlist$Bs
        ss<-ss+1
        #print(ss)
        #print(pow)
    }
    if (pow<power) warning("power less than nominal power, increase maxiter")
    #out<-list(B=Bs,power=pow,ssFixed=ssFixed,Nmax=max(Bs@N))
    #out
    Bs
}



