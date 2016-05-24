### Simulation code for unconstrained Bootstrap only
#################################################################
#   Initialization of parameters
#
# NSIM=number of replications
#NSIM<-10^5
NSIM<-10^4
#
#
# NBoot is the number of bootstrap replications
NBoot<-1000
#NBoot<-100
###################################################################
#   Part 1
#
#   Define functions for later use 
#
#
###################################################################
# most functions within bpcp package, load it
library(bpcp)
library(survival)

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


onestep<-function(d,ttrue,Strue,M=NBoot,NMC=NMonteCarlo){
    #r<-kmtestALL(d$time,d$status,ttrue,Strue,M=M,NMC=NMC)
    #r<-rJUNK
    #binomTest<-binTest(d$xi,d$ci,ttrue,Strue)
    #normLog<-normlogTest(d$time,d$status,ttrue,Strue)
    #modified<-modTest(d$time,d$status,ttrue,Strue)
    #check01<-function(x){ any(is.na(x)) | !all(x==0 | x==1) }
    #if (check01(modified)) browser()
    #reject<-rbind(binomTest,normLog,modified,r)

    #b<-bpcp(d$time,d$status,stype="mue")
    #bsci<-StCI(b,ttrue)
    #mue<-bsci[,"survival"]
    #bpci<-as.vector(as.matrix(bsci[,3:4]))
    #fit<-survfit(Surv(d$time,d$status)~1)
    #km<-StCI(fit,ttrue)[,"survival"]
    #kmInsideBPCP<-0
    #eps<-.Machine$double.eps^0.5
    #if ((bpci[1]<=km+eps & km-eps<=bpci[2])) kmInsideBPCP<-1

    eventsBeforeT<- length(d$xi[d$xi<ttrue & d$xi<d$ci])
    atRiskatT<- length(d$xi[d$ci>ttrue & d$xi>ttrue])
    output<-c(
        events=eventsBeforeT,
        atrisk=atRiskatT
        )
     reject<-kmtestBoot(d$time,d$status,ttrue,Strue,M=M)
     out<-list(reject=reject,output=output)

    out
}

## run everything once to get number of tests
d<-createData()
s<-onestep(d,3,.5)
Ntests<-1


onerow<-function(nsim,createD,ttrue,Strue,seed=12323511){
    set.seed(seed)
    R<-matrix(0,Ntests,3)
    OUT<-matrix(NA,nsim,2)
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
    }
    R<-100*R/nsim
    rowvalues<-c(
        TestTime=ttrue,
        Strue=Strue,
        MeanEvents=mean(OUT[,1]),
        MeanAtRisk=mean(OUT[,2]),
        low=R[,1],
        high=R[,2],
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
TABLE2<-matrix(NA,4,length(out),dimnames=
    list(c(rep("Exp(10) with U(0,5),n=30",4)),
           names(out)))

TABLE2[1,]<-onerow(NSIM,function(){createData(30)},1,1-pexp(1,1/10))
TABLE2[2,]<-onerow(NSIM,function(){createData(30)},2,1-pexp(2,1/10))
TABLE2[3,]<-onerow(NSIM,function(){createData(30)},3,1-pexp(3,1/10))

#### other rows give power
#TABLE2[5,]<-onerow(NSIM,function(){createData(30)},1,1-pexp(1,1/2.5))
#TABLE2[6,]<-onerow(NSIM,function(){createData(30)},2,1-pexp(2,1/2.5))
#TABLE2[7,]<-onerow(NSIM,function(){createData(30)},3,1-pexp(3,1/2.5))

#TABLE2[9,]<-onerow(NSIM,function(){createData(30)},1,1-pexp(1,1/100))
#TABLE2[10,]<-onerow(NSIM,function(){createData(30)},2,1-pexp(2,1/100))
#TABLE2[11,]<-onerow(NSIM,function(){createData(30)},3,1-pexp(3,1/100))


TABLE2[4,]<-onerow(NSIM,function(){createData(30)},4,1-pexp(4,1/10))
#TABLE2[8,]<-onerow(NSIM,function(){createData(30)},4,1-pexp(4,1/2.5))
#TABLE2[12,]<-onerow(NSIM,function(){createData(30)},4,1-pexp(4,1/100))

t1<-proc.time()
t1-t0

write.csv(TABLE2,file="SimulationTableBootOnly.csv",row.names=TRUE)


