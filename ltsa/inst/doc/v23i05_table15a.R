#Script for comparing FGN/ARMA forecast performance
#This script requires a Beowulf cluster computer with Rmpi library
#attach library
library(FGN)
#input data 
data(NileMin)
y<-NileMin
#ARMA(p,q) model order
p<-2
q<-1
#
nB<-10^4 #Number of bootstrap iterations
n<-length(y) #length of series
K<-100 #number of out-of-sample data values
n1<-n-K #length of training series
sdy<-sd(y) #sd of original series
#syd is the long-run prediction sd
outy<-FitFGN(y) #fix model for bootstrap
#
#combine all into together except outy,
inputdata=c(p,q,n,K)

onebootFGNARMA=function(outy,inputdata){
    MAXIT <- 10
    p <- inputdata[1]
    q <- inputdata[2]
    n <- inputdata[3]
    K <- inputdata[4]
    n1=n-K
#
#FGN fit to z1 and forecast using z2. 
#FGN and ARMA use independent bootstraps
    NotOK <- TRUE
    ITER<-0
    while (NotOK && ITER<MAXIT){
	z<-Boot(outy)
 	z1<-z[1:n1] #training data
   	z2<-z[-(1:n1)] #testing data
    	outz1<-tryCatch(FitFGN(z1), error = function(e) FALSE)
    	if (!(is.logical(outz1)))  NotOK <- FALSE
    	else ITER<-ITER+1
    }
    H<-outz1$H
    mu<-outz1$muHat
    rFGN<-var(z1)*FGNAcf(0:(n+3-1), H)
    F<-TrenchForecast(c(z1,z2), rFGN, mu, n1, maxLead=3)$Forecasts
    nF<-nrow(F)
    err1<-z2-F[,1][-nF]
    err2<-z2[-1]-F[,2][-c(nF,(nF-1))]
    err3<-z2[-c(1,2)]-F[,3][-c(nF,(nF-1),(nF-2))]
    rmse1<-sqrt(mean(err1^2))
    rmse2<-sqrt(mean(err2^2))
    rmse3<-sqrt(mean(err3^2))
    FGNrmse<-c(ITER,rmse1,rmse2,rmse3)
#
#ARMA(p,q) fit to z1 and forecast using z2
    NotOK <- TRUE
    ITER<-0
    while (NotOK && ITER<MAXIT){
	z<-Boot(outy)
 	z1<-z[1:n1] #training data
   	z2<-z[-(1:n1)] #testing data
    	outz1<-tryCatch(arima(z1,c(p,0,q)), error = function(e) FALSE)
    	if (!(is.logical(outz1)))  NotOK <- FALSE
    	else ITER<-ITER+1
    }
    z1<-z[1:n1] #training data
    z2<-z[-(1:n1)] #testing data
    phi<-theta<-numeric(0)
    if (p>0) phi<-coef(ansz1)[1:p]
    if (q>0) theta<-coef(ansz1)[(p+1):(p+q)]
    zm<-coef(ansz1)[p+q+1]
    sigma2<-ansz1$sigma2
    r<-var(z)*ARMAacf(ar=phi, ma=theta, lag.max=n+3-1)
    F<-TrenchForecast(c(z1,z2), r, zm, n1, maxLead=3)$Forecasts
    err1<-z2-F[,1][-nF]
    err2<-z2[-1]-F[,2][-c(nF,(nF-1))]
    err3<-z2[-c(1,2)]-F[,3][-c(nF,(nF-1),(nF-2))]
    rmse1<-sqrt(mean(err1^2))
    rmse2<-sqrt(mean(err2^2))
    rmse3<-sqrt(mean(err3^2))
    ARMArmse<-c(ITER, rmse1, rmse2, rmse3)
#Average mse
    c(FGNrmse, ARMArmse)
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
mpi.setup.rngstream(19000713)

startTime<-proc.time()
#start parallel bootstrap
out <- mpi.parReplicate(nB, onebootFGNARMA(outy=outy,inputdata=inputdata))
failed<-apply(out[c(1,5),], 1, sum)
averRMSE=apply(out[-c(1,5),], 1, mean)
endTime<-proc.time()
totalTime<-endTime-startTime
totalTime<-(totalTime/3600)[3]
#
#tabulate result
tb<-matrix(c(averRMSE[1:3],sdy,failed[1],averRMSE[4:6],sdy,failed[2]),ncol=2)
dimnames(tb)<-list(c("lead1","lead2","lead3","infty","failed"),c("FGN","ARMA"))
#
#output and close R
save(tb, out, totalTime, file="tbBoot.Rdata")  #
write(totalTime, "totalTimeMPI.txt")
write(t(round(tb,3)), file="tbMPI.txt", ncolumns=2)
print(paste("Elapsed Time = ", totalTime , "Hours"))
print(paste("nB = ",nB))
print(tb)

#close all slaves and quit
mpi.close.Rslaves()
mpi.quit()
