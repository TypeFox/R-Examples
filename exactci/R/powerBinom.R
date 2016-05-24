powerBinom<-function(n=NULL, p0=0.5, p1=NULL, sig.level=0.05, power=NULL, alternative=c("two.sided","one.sided"),
    strict=FALSE,type=c("standard","cilength","obs1ormore"),cilength=NULL,conf.level=0.95,direction=c("greater","less"),control=binomControl(),...){

    type<-match.arg(type)
    direction<-match.arg(direction)
    # not sure why this is better than list(...), but used in data.frame 
    # so I copied it
    dots<-as.list(substitute(list(...)))[-1L]

    if (type=="standard"){
        if (sum(sapply(list(n, p0, p1, power, sig.level), is.null)) != 
            1) stop("when type='standard' exactly one of 'n', 'p0', 'p1', 'power', and 'sig.level' must be NULL")
        alternative<-match.arg(alternative)
        if (alternative=="two.sided" & !strict){ 
            one.sided<-TRUE
            alpha<- sig.level/2
        } else if (alternative=="two.sided" & strict){
            one.sided<-FALSE
            alpha<-sig.level
        } else {
            one.sided<-TRUE
            alpha<-sig.level
        }
        if (is.null(power)){
            # in order to pass the results from ... (dots)  to the new call ..., we need to you do.call
            # power<-powerBinom.calcPower(n,p0,p1,alpha,one.sided,Control=control,...)$power
            power<-do.call("powerBinom.calcPower",c(list(n,p0,p1,alpha,one.sided,Control=control),dots))$power
        } 
        if (is.null(p1)){
            # in order to pass the results from ... (dots)  to the new call ..., we need to you do.call
            #p1<-powerBinom.calcP1(n,p0,alpha,one.sided,power,direction,...)
            p1<-do.call("powerBinom.calcP1",c(list(n,p0,alpha,one.sided,power,direction),dots))

        }
        if (is.null(n)){
            # in order to pass the results from ... (dots)  to the new call ..., we need to you do.call
            #nn<-powerBinom.calcN(p0,p1,alpha,one.sided,power,Control=control,...)
            nn<-do.call("powerBinom.calcN",c(list(p0,p1,alpha,one.sided,power,Control=control),dots))
            n<-nn$n
            power<-nn$power
        }
        method<-"power and sample size for single binomial response"
        if (strict & alternative=="two.sided"){ note<-"use rejections in opposite direction"
        } else { note<-"use rejections in correct direction only" } 
        out<-list(n=n,p0=p0,p1=p1,power=power,alternative=alternative,sig.level=sig.level,note=note,method=method)


    } else if (type=="cilength"){
          avgcilen<-function(n,theta=p1,...){
            dots<-as.list(substitute(list(...)))[-1L]
            if (theta==0){
                 X<-0
            } else if (theta==1){
                 X<-n
            } else {
                 X<-qbinom(control$pRange[1],n,theta):
                    qbinom(control$pRange[2],n,theta)
            }
            len<-rep(NA,length(X))

            for (i in 1:length(X)){
                # pass the dots  
                #ci<-binom.exact(X[i],n,conf.level=conf.level,...)$conf.int
                ci<-do.call("binom.exact",c(list(X[i],n,conf.level=conf.level),dots))$conf.int
                len[i]<-ci[2]-ci[1]
            }
            f<-dbinom(X,n,theta)
            f<-f/sum(f)
            sum(f*len)
          }
        if (is.null(p1)) p1<-0.5  
        if (is.null(n) & !is.null(cilength)){  
          cilength.minus.target<-function(n,theta=p1,...){
              dots<-as.list(substitute(list(...)))[-1L]
              do.call("avgcilen",c(list(n,theta),dots)) - cilength
          }      
        minn<-control$minn
        maxn<-control$maxn
        # fill in dots
        # cilength.minn<-cilength.minus.target(minn,theta=p1,...)
        cilength.minn<- do.call("cilength.minus.target",c(list(minn,theta=p1),dots))
        if (cilength.minn<0) stop("decrease control$minn or decrease 'cilength' ")
        # get a good starting step size
        # starting step size= 2^STEP.POWER
        # Note: 2*2*sqrt(p1(1-p1))/sqrt(n)=cilength   =>  n= (4*sqrt(p1*(1-p1))/cilength)^2   
        STEP.POWER<-ceiling(log2((4*sqrt(p1*(1-p1))/cilength)^2 ) )
        if (STEP.POWER>12){
            print("warning: calculation may take a long time")
            flush.console()
        }
        n<-try(uniroot.integer(cilength.minus.target,
            lower=control$minn,
            upper=control$maxn,
            step.up=TRUE,
            step.power=STEP.POWER,
            pos.side=TRUE,print.steps=control$PRINT.STEPS)$root)
        if (class(n)=="try-error") stop("increase control$maxn or 
               increase 'cilength' ")


        }  else if (is.null(cilength)){
            if (is.null(n)) stop("when type='cilength' must supply either 'cilength' or 'n'")
            # cilength<-avgcilen(n,theta=p1,...)
            cilength<-do.call("avgcilen",c(list(n,theta=p1),dots))
        }
        method<-"Minimum sample size so that the expected CI length is cilength, 
               when the true proportion is p1"
        note<-"sample size is maximized when p1=0.5"
        out<-list(n=n,p1=p1,cilength=cilength,conf.level=conf.level,method=method,note=note)
    } else if (type=="obs1ormore"){
        # probability to observe at least one event
        # use equation:   power = 1- (1-p1)^n
        # solve for unknown
        if (sum(sapply(list(n, p1, power), is.null)) != 
            1) stop("when type='obs1ormore' exactly one of 'n', 'p1', and 'power' must be NULL")
        if (is.null(n)){
            n<- ceiling( log(1-power)/log(1-p1) )
        } else if (is.null(power)){
            power<- 1 - (1-p1)^n
        } else if (is.null(p1)){
            p1<- 1-exp( log(1-power)/n )
        }
        method<-"power (i.e., probability) to observe at least one event from binomial(n,p)"
        note<-""
        out<-list(n=n,p=p1,power=power,method=method)
    }
    ## if midp is passed through ... then add to note
    if (length(dots)>0){
        if (any(names(dots)=="midp")){
            if (dots[names(dots)=="midp"]$midp){
                out$note<-paste(note,"
                    midp=TRUE",collapse="")
            }
        }
    }
    class(out)<-"power.htest"
    out
}
#library(ssanv)
#library(exactci)
#powerBinom(power=.9,p1=.01,type="obs1ormore")
#powerBinom(cilength=0.1,p1=.5,type="cilength")


powerBinom.calcPower<-function(n,theta0=.5,theta1=0.7,
     alpha=0.025,one.sided=TRUE,fast=FALSE,Control,...){
    # use 'Control' here to allow a second 'control' to be passed 
    # to binom.exact through '...' if needed

    dots<-as.list(substitute(list(...)))[-1L]
    reject<-rep(FALSE,n+1)
    if (fast){
           if (theta1==0){
                 X<-0
            } else if (theta1==1){
                 X<-n
            } else {
                 X<-qbinom(Control$pRange[1],n,theta1):
                    qbinom(Control$pRange[2],n,theta1)
            }
    } else { X<-0:n }
    if (one.sided){
       for (i in X){
          alt<-ifelse(theta1>theta0,"greater","less")
          # put results in reject[i+1] because goes from 0:n
          # use binom.exact so that midp option can be passed
          # use do.call to pass ...
          # pvalue<-binom.exact(i,n,p=theta0,alternative=alt,...)$p.value
          pvalue<-do.call("binom.exact",c(list(i,n,p=theta0,alternative=alt),dots))$p.value
          if (pvalue<=alpha) reject[i+1]<-TRUE
       }
    } else {
       for (i in X){
          # use do.call to pass ...
          # pvalue<-binom.exact(i,n,p=theta0,...)$p.value
          pvalue<-do.call("binom.exact",c(list(i,n,p=theta0),dots))$p.value
          if (pvalue<=alpha) reject[i+1]<-TRUE
       }
    }
    ## if fast only count rejections 
    ## from qbinom(pRange[1],n,theta1):qbinom(pRange[2],n,theta1)
    list(power=sum(dbinom((0:n)[reject],n,theta1)),reject=reject)
}

powerBinom.calcP1<-function(n,theta0=.5,alpha=0.025,one.sided=TRUE,power=.80,direction="greater",...){
    dots<-as.list(substitute(list(...)))[-1L]
    # use powerBinom.calcPower to get reject region, use theta1=.5 but we will not use power from this calculation
    theta1<-ifelse(direction=="greater",1,0)
    # need to use fast=FALSE so we get the full range of reject, since 
    # theta1 will be changing during uniroot call below
    # Note: do not need 'Control' since fast=FALSE 
    # Note: pass ... using do.call
    #r<-powerBinom.calcPower(n,theta0,theta1,alpha,one.sided,fast=FALSE,Control=NULL,...)$reject
    r<-do.call("powerBinom.calcPower",c(list(n,theta0,theta1,alpha,one.sided,fast=FALSE,Control=NULL),dots))$reject
    rootfunc<-function(p,reject=r,nominalPower=power){
         n<-length(reject)-1
         sum(dbinom((0:n)[reject],n,p)) -nominalPower
    }
    if (direction=="greater"){
        p1<-uniroot(rootfunc,c(theta0,1))$root
    } else {
        p1<-uniroot(rootfunc,c(0,theta0))$root
    } 
    p1
}

powerBinom.calcN<-function(theta0,theta1,alpha,one.sided,
        power,Control,...){
    dots<-as.list(substitute(list(...)))[-1L]
    minn<-Control$minn
    maxn<-Control$maxn
    N<-minn
    # Note: pass ... using do.call
    #pow0<-powerBinom.calcPower(N,theta0,theta1,alpha,one.sided,fast=TRUE,Control=Control,...)$power
    pow0<-do.call("powerBinom.calcPower",c(list(N,theta0,theta1,alpha,one.sided,fast=TRUE,Control=Control),dots))$power
    if (pow0<power){
        for (i in N:maxn){
            N<-i
            #pow0<-powerBinom.calcPower(N,theta0,theta1,alpha,one.sided,fast=TRUE,Control=Control,...)$power
            pow0<-do.call("powerBinom.calcPower",c(list(N,theta0,theta1,alpha,one.sided,fast=TRUE,Control=Control),dots))$power
            if (pow0>=power) break()
        }
    }
    if (N==maxn & pow0<power) warning("N equals maxn and power is less than nominal power")
    list(n=N,power=pow0,minn=minn,maxn=maxn)
}
#powerBinom.calcPower(2,.7,.5,.025,TRUE)
#powerBinom.calcN(.7,.5,.025,TRUE,.8)

#t0<-proc.time()
#source("binomControl.R")
#powerBinom(n=NULL, p0=0.5, p1=.6, sig.level=0.05, power=.80, alternative=c("two.sided","one.sided"),
#    strict=FALSE,type="standard",cilength=NULL,conf.level=0.95,control=binomControl())
#t1<-proc.time()
#t1-t0
#powerBinom(n=90, p0=0.5, p1=NULL, sig.level=0.05, power=.8, direction="less")

#powerBinom(type="cilength",cilength=.3)