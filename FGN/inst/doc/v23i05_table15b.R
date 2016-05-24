#Script for comparing FGN/ARMA forecast performance
#This script requires a Beowulf cluster computer with Rmpi library
#attach library
library(FGN)
#input data 
data(NileMin)
y<-NileMin
nB<-10^4 #Number of bootstrap iterations
n<-length(y) #length of series
K<-100 #number of out-of-sample data values
n1<-n-K #length of training series
sdy<-sd(y) #sd of original series
#syd is the long-run prediction sd
outy<-FitFGN(y) #fix model for bootstrap
#
#combine all into together except outy,
inputdata=c(n,K)

onebootFGNARMA=function(outy,inputdata){
    n <- inputdata[1]
    K <- inputdata[2]
    n1=n-K
#
#FGN fit to z1 and forecast using z2. 
#FGN and ARMA use independent bootstraps
    z<-Boot(outy)
    z1<-z[1:n1] #training data
    z2<-z[-(1:n1)] #testing data
    H<-outy$H  #same as H=0.831 (original estimate for NileMin complete series)
    mu<-outz1$muHat
    maxLead<-3
    rFGN<-var(z1)*FGNAcf(0:(n+maxLead-1), H)
    ans<-TrenchForecast(c(z1,z2), rFGN, mu, n1, maxLead=maxLead)
    F<-ans$Forecasts
    nF<-nrow(F)
    err1<-z2-F[,1][-nF]
    err2<-z2[-1]-F[,2][-c(nF,(nF-1))]
    err3<-z2[-c(1,2)]-F[,3][-c(nF,(nF-1),(nF-2))]
    rmse1<-sqrt(mean(err1^2))
    rmse2<-sqrt(mean(err2^2))
    rmse3<-sqrt(mean(err3^2))
    FGNrmse<-c(rmse1,rmse2,rmse3)
#
#Average mse
    FGNrmse
}

####################################################
#start Rmpi and setup R slaves
library(Rmpi)
#setup slaves
mpi.spawn.Rslaves()
#send data and function to all slaves
mpi.bcast.Robj2slave(inputdata)
mpi.bcast.Robj2slave(outy)
mpi.bcast.Robj2slave(onebootFGNARMA)

#tell slaves to load FGN
mpi.bcast.cmd(library(FGN))

#setup parallel RNG. seed can be specified.
mpi.setup.rngstream(19480813)

startTime<-proc.time()
#start parallel bootstrap
out <- mpi.parReplicate(nB, onebootFGNARMA(outy=outy,inputdata=inputdata))
averRMSE <- apply(out, 1, mean)
endTime<-proc.time()
totalTime<-endTime-startTime
totalTime<-(totalTime/3600)[3]
#
#tabulate result
tb<-averRMSE
names(tb)<-c("lead1","lead2","lead3")
print(tb)
#
#output and close R
save(tb, out, totalTime, file="tbParameterCertainty.R")  #
write(round(tb,3), file="tbParUn.txt", ncolumns=2)
print(paste("Elapsed Time = ", totalTime , "Hours"))
print(paste("nB = ",nB))

#close all slaves
mpi.close.Rslaves()
mpi.quit()
