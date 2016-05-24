#Script for Empirical Variance of H
#This script requires a Beowulf cluster computer with Rmpi library
nB<-10^5 #Number of bootstrap iterations
ns<-c(100, 200, 500, 1000, 2000)
HS<-c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
#
oneSimulation<-function(ns, HS){
    nns<-length(ns)
    nHS<-length(HS)
    H<-numeric(nns*nHS)
    ij<-0
    for (i in 1:length(ns)){
    	n<-ns[i]
    	for (j in 1:length(HS)){
        	ij<-ij+1
        	H0<-HS[j]
        	z<-SimulateFGN(n, H0)
        	H[ij]<-GetFitFGN(z)$H
    	}
    }
    H
}

#start Rmpi and setup R slaves
library(Rmpi)
#setup slaves
mpi.spawn.Rslaves()
#send data and function to all slaves
mpi.bcast.Robj2slave(ns)
mpi.bcast.Robj2slave(HS)
mpi.bcast.Robj2slave(oneSimulation)

#tell slaves to load FGN
mpi.bcast.cmd(library(FGN))

#setup parallel RNG. seed can be specified.
mpi.setup.rngstream(19480813)

startTime<-proc.time()
#start parallel bootstrap
out <- mpi.parReplicate(nB, oneSimulation(ns=ns, HS=HS))
averVarH <- apply(out, 1, var)
endTime<-proc.time()
totalTime<-endTime-startTime
totalTime<-(totalTime/3600)[3]
#
#tabulate result
tb<-rep(ns, rep(length(HS),length(ns)))*averVarH
tb<-matrix(tb, ncol=length(ns))
dimnames(tb)<-list(HS,ns)
#
#output and close R
save(tb, out, totalTime, file="tbVarH.R")  #
print(paste("Completed Var HHat Simulations. Elapsed Time = ", totalTime , "Hours"))
print(paste("nB = ",nB))
print(tb)

#close all slaves
mpi.close.Rslaves()
mpi.quit()

