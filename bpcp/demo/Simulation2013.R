### Simulation code
###
###     Increase parameters in initialization section for more replications, 
###     and more methods, see comments    
###
###     Simulation for unconstrained Bootstrap is in a separate file, since it was requested by a reviewer after the initial submission  
#################################################################
#   Initialization of parameters
#
# NSIM=number of replications
#NSIM<-10^5
NSIM<-10^2
#
#
# NMonteCarlo is the number of Monte Carlo replications for the BPCP MC method
# NMonteCarlo=0, just means do not do the BPCP  MC method
#NMonteCarlo<-10^5
NMonteCarlo<-0
#
# NBoot is the number of bootstrap replications, NBoot=0 means do not do bootstrap method
#NBoot<-1000
NBoot<-0
###################################################################
#   Part 1
#
#   Define functions for later use 
#
#
###################################################################
# most functions within bpcp package, load it
library(bpcp)
# time how long simulation takes
t0<-proc.time()


createData<-function(n=30,muX=10,urange=c(0,5)){
    x<-rexp(n,1/muX)
    cens<-runif(n,min=urange[1],max=urange[2])
    time<-pmin(x,cens)
    status<-rep(0,n)
    status[time==x]<-1
    list(time=time,status=status,xi=x,ci=cens)
}

createData2<-function(n=30,muX=10){
    x<-rexp(n,1/muX)
    cens<-rep(Inf,n)
    time<-pmin(x,cens)
    status<-rep(0,n)
    status[time==x]<-1
    list(time=time,status=status,xi=x,ci=cens)
}
binTest<-function(xi,ci,ttrue,Strue){
    I<- ci>ttrue
    x<-xi[I]
    Y<- length(x[x>ttrue])
    N<-length(x)
    if (N==0){ bci<-c(0,1) 
    } else bci<-binom.test(Y,N)$conf.int
    reject<-rejectFromInt(Strue,bci,thetaParm=TRUE)
    reject
}
### need survival Package for modTest function
library(survival)
modTest<-function(time,status,ttrue,Strue){
    mod<-survfit(Surv(time,status)~1,conf.lower='modified')
    ## if last observations is a failure, then gives NA for 
    ## last lower and upper
    ## replace lower with 0, upper with previous upper
    if (any(is.na(mod$lower))){
         k<-length(mod$lower)
         imax<-max((1:k)[!is.na(mod$upper)])
         mod$lower[is.na(mod$lower)]<-0
         mod$upper[is.na(mod$upper)]<-mod$upper[imax]
    }
    sci<-StCI(mod,ttrue,afterMax="zero")
    reject<-rejectFromInt(Strue,sci[3:4],thetaParm=TRUE)
    reject
}
normlogTest<-function(time,status,ttrue,Strue){
    mod<-survfit(Surv(time,status)~1)
    sci<-StCI(mod,ttrue,afterMax="zero")
    reject<-rejectFromInt(Strue,sci[3:4],thetaParm=TRUE)
    reject
}


#binomTest<-binTest(d$xi,d$ci,3,SMix(3))
#modTest(x$time,x$status,5,.64)

##################################################
# For the simulations for the paper, I did the simulations with 
# the bootstrap and BPBP MC first and most of the methods first. 
# This took over a week on the biowulf network of computers
#
# Later I went back and did the modified, and binomTest
# to save time when I went back, I created the rJUNK output
# and just substituted it in each replication of the simulation
# instead of redoing the kmtestALL again 
#rJUNK<-kmtestALL(d$time,d$status,3,SMix(3),M=0,NMC=0)
#rbind(binomTest,r)

onestep<-function(d,ttrue,Strue,M=NBoot,NMC=NMonteCarlo){
    r<-kmtestALL(d$time,d$status,ttrue,Strue,M=M,NMC=NMC)
    #r<-rJUNK
    binomTest<-binTest(d$xi,d$ci,ttrue,Strue)
    normLog<-normlogTest(d$time,d$status,ttrue,Strue)
    modified<-modTest(d$time,d$status,ttrue,Strue)
    #check01<-function(x){ any(is.na(x)) | !all(x==0 | x==1) }
    #if (check01(modified)) browser()
    reject<-rbind(binomTest,normLog,modified,r)

    b<-bpcp(d$time,d$status,stype="mue")
    bsci<-StCI(b,ttrue)
    mue<-bsci[,"survival"]
    bpci<-as.vector(as.matrix(bsci[,3:4]))
    fit<-survfit(Surv(d$time,d$status)~1)
    km<-StCI(fit,ttrue)[,"survival"]
    kmInsideBPCP<-0
    eps<-.Machine$double.eps^0.5
    if ((bpci[1]<=km+eps & km-eps<=bpci[2])) kmInsideBPCP<-1

    eventsBeforeT<- length(d$xi[d$xi<ttrue & d$xi<d$ci])
    atRiskatT<- length(d$xi[d$ci>ttrue & d$xi>ttrue])
    output<-c(
        events=eventsBeforeT,
        atrisk=atRiskatT,
        mue=mue,
        km=km,
        kminside=kmInsideBPCP,
        bplo=bpci[1],
        bphi=bpci[2])
     out<-list(reject=reject,output=output)

    out
}

## run everything once to get number of tests
d<-createData()
s<-onestep(d,3,.5)
Ntests<-nrow(s$reject)

MSE<-function(x,xtrue){
   mean( (x-xtrue)^2 )
}

vectWithNames<-function(mat,byrow=TRUE){
    ### write matrix as vector but name each element with rname:cname
    rnames<-dimnames(mat)[[1]]
    cnames<-dimnames(mat)[[2]]
    nr<-nrow(mat)
    nc<-ncol(mat)
    if (byrow) mat<-t(mat)
    outNames<-paste(rep(rnames,each=nc),rep(cnames,nr),sep=":")
    out<-as.vector(mat)
    names(out)<-outNames
    out
}


onerow<-function(nsim,createD,ttrue,Strue,seed=12323511){
    set.seed(seed)
    R<-matrix(0,Ntests,3)
    OUT<-matrix(NA,nsim,7)
    numNoDeath<-0
    for (i in 1:nsim){
        #if (i>12) browser()

        d<-createD()
        while (all(d$status==0)){
             d<-createD()
             numNoDeath<-numNoDeath+1
             if (numNoDeath>100) break()
        }
        s<-try(onestep(d,ttrue,Strue))
        if (class(s)=="try-error") browser()
        R<-R+s$reject
        OUT[i,]<-s$output
        if (s$output["kminside"]==0) browser()
    }
    R<-100*R/nsim
    rowvalues<-c(
        TestTime=ttrue,
        MeanEvents=mean(OUT[,1]),
        MeanAtRisk=mean(OUT[,2]),
        MSEmue=1000*MSE(OUT[,3],Strue),
        MSEkm=1000*MSE(OUT[,4],Strue),
        KMInside=100*sum(OUT[,5])/nsim,
        vectWithNames(R[,-3]),
        numNoDeath=numNoDeath)
    rowvalues
}
####################################################################
#   Part 2
#
#   run simulation
#
#####################################################################
## do one replication just to get row names, etc
out<-onerow(1,function(){createData(30)},1,1-pexp(1,1/10))
TABLE2<-matrix(NA,12,length(out),dimnames=
    list(c(rep("Exp(10) with U(0,5),n=30",12)),
           names(out)))

TABLE2[1,]<-onerow(NSIM,function(){createData(30)},1,1-pexp(1,1/10))
TABLE2[2,]<-onerow(NSIM,function(){createData(30)},2,1-pexp(2,1/10))
TABLE2[3,]<-onerow(NSIM,function(){createData(30)},3,1-pexp(3,1/10))

TABLE2[5,]<-onerow(NSIM,function(){createData(30)},1,1-pexp(1,1/2.5))
TABLE2[6,]<-onerow(NSIM,function(){createData(30)},2,1-pexp(2,1/2.5))
TABLE2[7,]<-onerow(NSIM,function(){createData(30)},3,1-pexp(3,1/2.5))

TABLE2[9,]<-onerow(NSIM,function(){createData(30)},1,1-pexp(1,1/100))
TABLE2[10,]<-onerow(NSIM,function(){createData(30)},2,1-pexp(2,1/100))
TABLE2[11,]<-onerow(NSIM,function(){createData(30)},3,1-pexp(3,1/100))


TABLE2[4,]<-onerow(NSIM,function(){createData(30)},4,1-pexp(4,1/10))
TABLE2[8,]<-onerow(NSIM,function(){createData(30)},4,1-pexp(4,1/2.5))
TABLE2[12,]<-onerow(NSIM,function(){createData(30)},4,1-pexp(4,1/100))

t1<-proc.time()
t1-t0

write.csv(TABLE2,file="SimulationTable.csv",row.names=TRUE)


