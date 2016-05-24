
## n <- 1e4; ACE <- c(1,1,1)/3
## logscale <- -4.5; logshape <- .7
## p2 <- .065
## a2 <- -10; b2 <- 0.15 
## a1 <- 85; b1 <- 0.1

## pmvn(upper=c(q2,q2),sigma=diag(2)*(1-2/3)+2/3)
## pmvn(c(q2,q2),sigma=diag(2)*(1-2/3)+2/3)

## pnorm(q2,sd=1) ## Marginal / Perfect dependence
## pmvn(upper=c(q2,q2),sigma=diag(2)*(1-2/3)+2/3) ## Concordance
## pnorm(q2,sd=1)^2 ## Independence
## (lambdaR <- pmvn(upper=c(q2,q2),sigma=diag(2)*(1-2/3)+2/3)/pnorm(q2,sd=1)^2)

bicomprisksim <- function(n=1e4,
                          ACE=c(1/3,1/3,1/3),
                          logscale=-4.5,logshape=.7,
                          a1=85,b1=0.1,
                          a2=-10,b2=0.15,
                          p2=.065,
                          tt,
                          ...) {

    ACE <- ACE/sum(ACE)
    ialpha <- function(v,a,b) -(log(-v)+a)/b    
    q2 <- qnorm(p2) ## p2: Prostata cancer prevalence 
    p1 <- 1-p2; q1 <- qnorm(p1) ## Death without cancer
    alpha <- function(t,a,b) -exp(-(b*t+a))
    F2s <- function(t) pnorm(alpha(t,b=b2,a=a2)+q2) ## Marginal Cumulative Incidence (cancer)
### Random effects
    R <- diag(2)*.5+.5
    J <- matrix(1,ncol=2,nrow=2)
    I <- diag(2)
    zyg <- rep(c(0,1),each=n)
    A <- rbind(rmvn(n,sigma=R*ACE[1]),rmvn(n,sigma=J*ACE[1]))
    C <- rmvn(2*n,sigma=J*ACE[2])
### Random effects 'death'
    eta1 <- C
### Random effects 'cancer'
    eta2 <- A+C;
### Subject-specific probability of lifetime cancer 
    probcanc <- pnorm(q2+eta2,sd=ACE[3]^.5)   
### Cancer/Death without cancer realizations
    cancertrue <- (runif(length(probcanc))<probcanc)*1+1
### Event times given failure and random effect
    ## inverse P(T<t| eta,cause=1)
    iF1 <- function(u,eta,sd,a,b) (qnorm(u,sd=sd)-eta)/b+a 
    iF2 <- function(u,eta,pr,sd,q,a,b) ialpha(qnorm(u*pr,sd=sd)-eta-q,a=a,b=b)
    u <- runif(length(probcanc))
    t2 <- iF2(u,eta2,probcanc,sd=sqrt(1-sum(ACE[1:2])),q=q2,a=a2,b=b2)
    t1 <- iF1(u,eta1,sd=sqrt(1-sum(ACE[2])),a=a1,b=b1)   
    t0 <- t1; t0[which(cancertrue==2)] <- t2[which(cancertrue==2)] 
### Censoring distributon
    cens <- matrix(rep(rweibull(length(probcanc)/2,shape=exp(logshape),scale=1/exp(logscale)),2),
                   ncol=2)
    ##hist(cens,xlab=c(0,400),200)    
    cause <- cancertrue
    cause[t0>cens] <- 0
    t <- pmin(t0,cens)
    (censtab <- table(cause)/length(cause))    
### Data.frame
    d <- data.frame(t,cause,zyg,cancertrue-1,cens[,1],t0[,1],t0[,2]);
    names(d) <- c("time1","time2","cause1","cause2","zyg","cancertrue1","cancertrue2","cens.time","T01","T02")
### Long format
    dd <- fast.reshape(d)
    dd$cancer <- (dd$cause==2)*1

    if (missing(tt)) tt <- seq(0,max(dd$time))    
    Smz <- J*ACE[1]+J*ACE[2]+I*ACE[3]
    Sdz <- R*ACE[1]+J*ACE[2]+I*ACE[3]
    rr <- alpha(tt,b=b2,a=a2)+q2
    Cmz <- pmvn(upper=cbind(rr,rr),mu=matrix(0,ncol=2,nrow=length(rr)),sigma=Smz)
    Cdz <- pmvn(upper=cbind(rr,rr),mu=matrix(0,ncol=2,nrow=length(rr)),sigma=Sdz)
    
    true <- list(p2=p2,
                 p22mz=pmvn(lower=c(-Inf,-Inf),upper=c(q2,q2),sigma=Smz),
                 p12mz=pmvn(lower=c(q2,-Inf),upper=c(Inf,q2),sigma=Smz),
                 p12mz=pmvn(lower=c(q2,-Inf),upper=c(Inf,q2),sigma=Smz),
                 p11mz=pmvn(lower=c(q2,q2),upper=c(Inf,Inf),sigma=Smz),
                 p22dz=pmvn(lower=c(-Inf,-Inf),upper=c(q2,q2),sigma=Sdz),
                 p12dz=pmvn(lower=c(q2,-Inf),upper=c(Inf,q2),sigma=Sdz),
                 p12dz=pmvn(lower=c(q2,-Inf),upper=c(Inf,q2),sigma=Sdz),
                 p11dz=pmvn(lower=c(q2,q2),upper=c(Inf,Inf),sigma=Sdz),
                 time=tt,
                 F2=F2s(tt),
                 Cmz=Cmz,
                 Cdz=Cdz
                 )
    attributes(dd) <- c(attributes(dd),true)                  
    return(dd)
}




##################################################
##################################################
##################################################
##################################################

##################################################
### simulation for gamma distributed cif model
##################################################
###
##### {{{ 
###lap<-function(theta,t) {
###    return( (1+t/theta)^(-theta))
###}
###ilap<-function(theta,t) {
###    itheta<-1/theta; return((t^(-itheta)-1)/(itheta))
###}
###
###F1clust<-function(t,rtheta=1,theta=1,lam0=0.5,beta=0.3,x=0) {
###    return(1-exp(-rtheta*ilap(theta,exp(-t*lam0-t*x*beta))))
###}
###
######################################################################
###F1<-function(t,lam0=0.5,beta=0.3,x=0) { # additive version
###    return( 1 - exp(-t*lam0-t*x*beta))
###}
###
######F1<-function(t,lam0=0.5,beta=0.3,x=0) # proportional version 
######{ return( 1 - exp(-(t*lam0)*exp(x*beta))) }
###
###sim.F1<-function(n,theta=1,lam0=0.5,beta=0.3,crate=2) {
###    x<-runif(n); tt<-seq(0,1,length=100)
###    F11x<-F1(1,x=x,beta=beta,lam0=lam0)
###    cause1<-rbinom(n,1,F11x)
###
###    stime<-rep(100,n); 
###    for (i in 1:n)
###        {
###            if (cause1[i]==1) {
###                myhazx<-F1(tt,x=x[i],beta=beta,lam0=lam0)/F11x[i]
###                stime[i]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
###            } 
###        }
###    ctime<-runif(n)*crate
###    time<-apply(cbind(ctime,stime),1,min)
###    status<-(stime<ctime); cause1[status==0]<-0; 
###    data<-data.frame(time=time,status=status,X=x,cause=cause1)
###    return(data)
###} 
###
###sim.F1clust<-function(n,theta=1,lam0=0.5,beta=0.3,alpha=0,crate=2,
###                      same.cens=FALSE,fix.cens=FALSE) {
###    k<-n/2; tt<-seq(0,1,length=100)
###    rtheta<-rgamma(k,theta,scale=1/theta)
###    stime<-c();cause1<-c();id<-c();vtheta<-c(); X<-c()
###    cause1 <- c()
###
###    for (i in 1:k) { 
###        x<-runif(2); 
###        X<-c(X,x); 
###        F11x<-F1clust(1,rtheta=rtheta[i],theta=theta,x=x,beta=beta,lam0=lam0) 
###        cause<-rbinom(2,1,F11x); 
###        ##cause1<-c(cause1,cause); 
###        id<-c(id,rep(i,2)); vtheta<-c(vtheta,rep(rtheta[i],2))
###        
###        for (j in 1:2) {
###            if (cause[j]==1) {
###                myhazx<-F1clust(tt,x=x[j],rtheta=rtheta[i],theta=theta, beta=beta,lam0=lam0)/F11x[j]
###                stime<-c(stime,Cpred(cbind(myhazx,tt),runif(1))[1,2]+ runif(1,0,0.001))
###                cause1 <- c(cause1,1);
###            } else { stime<-c(stime,runif(1)); cause1 <- c(cause1,2);}
###        }
###    }
###    
###    if (same.cens) ctime <- rep(ctime<-runif(k)*crate,each=2) else ctime<-runif(n)*crate;
###    if (fix.cens) ctime <- rep(ctime<-rep(crate,each=2))
###
###    time<-apply(cbind(ctime,stime),1,min)
###    status<-(stime<ctime); 
###    cause1[status==0]<-0; 
###    data<-data.frame(time=time,X=X,cause=cause1,stime=stime,ctime=ctime,
###                     id=id,theta=vtheta)
###    return(data)
###} 
###
##### }}} 
###
######################################################
### test ting  #######################################
######################################################

## beta=0; n=200; theta=0.1; 
## simdata<-sim.F1clust(n=200,theta=0.1,beta=0,crate=3.0,lam0=0.5)
## table(simdata$cause)
## simdata$tv <- rep(1:2,n/2)

## tt <- seq(0,1,by=0.01)
## F1(tt)
 
## p11<- 1 - (1-F1(tt)) - (1-F1(tt)) + lap(1/theta, 2*ilap(1/theta,1-F1(tt)) )
## plot(tt,F1(tt),type="l")
## lines(tt,p11)
## lines(tt,(F1(tt)^2),col=2)

## beta=0; n=1000; theta=0.1; 
## gem <- c()
## for (i in 1:10) {
## simdata<-sim.F1clust(n,theta=theta,beta=beta,crate=3.0,lam0=0.5,same.cens=TRUE)
## table(simdata$cause)
## simdata$tv <- rep(1:2,n/2)
## simd <- reshape(simdata,direction="wide",idvar="id",timevar="tv")
## simd$cancer <- (simd$cause.1==1 ) + (simd$cause.2==1 )  
## tabd <- table(simd$cancer)
## casewise <- tabd["2"]/ ( tabd["2"] + tabd["1"]*0.5)  
## gem <- c(gem,casewise)
## }
## summary(gem)

## ### plot

## casewiset <- p11/F1(tt)
## plot(tt,casewiset,type="l")

## F1t <- F1(tt)


###################################################################
###################################################################

### increase censoring by crate lower, dependence theta lower  

## theta=1.0;  
## p11t<- 1 - (1-F1(tt)) - (1-F1(tt)) + lap(theta, 2*ilap(theta,1-F1(tt)) )
## F1t <- F1(tt)
## Gc <- 1 - 2*tt
## Gc[Gc<0] <- 0
## mn11t <- cumsum(diff(p11t)* Gc[-1])
## mndt <- cumsum(diff(F1t)* Gc[-1])
## casewisem <- mn11t[100]/mndt[100]
## casewiset <- p11t/F1(tt)

## gem <- c()
## for (i in 1:10) {
## print(i)
## beta=0; n=10000; 
## theta=theta;  
## simdata<-sim.F1clust(n,theta=theta,beta=beta,crate=0.5,lam0=0.5,same.cens=TRUE,fix.cens=FALSE)
## table(simdata$cause)
## simdata$tv <- rep(1:2,n/2)
## simd <- reshape(simdata,direction="wide",idvar="id",timevar="tv")
## simd$cancer <- (simd$cause.1==1 ) + (simd$cause.2==1 )  
## tabd <- table(simd$cancer)
## casewise <- tabd["2"]/ ( tabd["2"] + tabd["1"]*0.5)  
## ###
## outm <- prodlim(Event(time,cause)~+1,data=simdata)
## par(mfrow=c(1,3))
## plot(outm)
## lines(tt,F1t,col=3,lwd=3)
## legend("topleft",c("cif","sand"),lty=1,col=c(1,3))
## ###          
## cc <- bicomprisk(Event(time,cause)~+1+id(id),data=simdata,cause=c(1,1),prodlim=TRUE)
## p11t<- 1 - (1-F1(tt)) - (1-F1(tt)) + lap(theta, 2*ilap(theta,1-F1(tt)) )
## plot(cc)
## lines(tt,p11t,col=3,lwd=3)
## lines(tt,F1t^2,col=4,lwd=1)
## legend("topleft",c("conc","sand","indep"),lty=1,col=c(1,3,4))
## ###
## case <- casewise(cc,outm,cause.marg=1)
## ###
## casewiset <- p11t/F1t
## plot(case$casewise[,1],case$casewise[,2],type="l")
## lines(tt,casewiset,col=3,lwd=3)
## abline(h=casewise,col=2)
## legend("topleft",c("est","sand","constant"),lty=1,col=c(1,3,2))
## ###
## Gc <- 1 - 2*tt
## Gc[Gc<0] <- 0
## ###
## abline(h=casewisem,col=4)
## gem <- c(gem,casewise)
## }

## abline(h=casewisem,col=4)
## casewisem
## mean(gem)

## pdf("caswise-vs-surv.pdf")
## par(mfrow=c(1,1))

## dev.off()

############################################################################
############################################################################
############################################################################

## library(doMC)
## registerDoMC()
## library(mets)

## onesim <- function(i,n,theta=2,beta=0) {
## 	theta=2; beta=0; 
## 	n=400

##         simdata<-sim.F1clust(n,theta=theta,beta=beta,crate=0.5,lam0=0.5,same.cens=TRUE,fix.cens=FALSE)
##         simdata$tv <- rep(1:2,n/2)
##         simd <- reshape(simdata,direction="wide",idvar="id",timevar="tv")
##         simd$cancer <- (simd$cause.1==1 ) + (simd$cause.2==1 )  
##         tabd <- table(simd$cancer)
##         casewise <- tabd["2"]/ ( tabd["2"] + tabd["1"]*0.5)  
## ###
##         outm <- prodlim(Hist(time,cause)~+1,data=simdata)
##         par(mfrow=c(1,3))
##         plot(outm)
##         lines(tt,F1t,col=3,lwd=3)
##         legend("topleft",c("cif","sand"),lty=1,col=c(1,3))
## ###          
##         cc <- bicomprisk(Event(time,cause)~+1+id(id),data=simdata,cause=c(1,1),prodlim=TRUE)
##         p11t<- 1 - (1-F1(tt)) - (1-F1(tt)) + lap(theta, 2*ilap(theta,1-F1(tt)) )
##         case <- casewise(cc,outm,cause.marg=1)
##         case$concordance
##         case$casewise
##         simdata$time
        
##         times <- seq(0,0.4,by=0.01)
##         outm <-comp.risk(Surv(time,cause==0)~+1+cluster(id),data=simdata,simdata$cause,
## 		                         causeS=1,times=times,max.clust=200,n.sim=0)
##         cifmz <-predict(outm,X=1,uniform=0,resample.iid=1)
##         cc <-bicomprisk(Event(time,status)~+1+id(id),data=simdata,
##                         cause=c(1,1),se.clusters=outm$clusters)
##         cdz <- cc$model$"DZ"
##         cmz <- cc$model$"MZ"
        
##         cdz <- casewise.test(cdz,cifmz,test="case") ## test based on casewise
##         cmz <- casewise.test(cmz,cifmz,test="conc") ## based on concordance:w
        
##         out <- comp.risk(Event(time,cause)~+1,data=simdata,causeS=1)
##         cc <- bicomprisk(Event(time,cause)~+1+id(id),data=simdata,cause=c(1,1))
        
        
##     }


