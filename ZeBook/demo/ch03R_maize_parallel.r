################################################################################
# ZeBook. Working with dynamic models for agriculture
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2010-08-09
# R introduction. Improving performance with parallel code : example with maize.model
library(ZeBook)
library(parallel)
# example with a simulation design with random generation with an uniform distribution
n <- 100
mat.param <- sapply(1:7,function(k)runif(n,min=maize.define.param()["binf",k],
   max=maize.define.param()["bsup",k]))
colnames(mat.param) <- colnames(maize.define.param())
mat.param=cbind(id=1:nrow(mat.param),mat.param)
# for one single weather, sdate and ldate
weather <- maize.weather(working.year=2010, working.site=1)
sdate <- 100
ldate <- 250

# 1) define the function that is run in each parallel PU
# id is used to construct arguments and the input variables required
# sdate,ldate, ncore are arguments to be passed to the parallel process
maize_cluster <-function (id, sdate,ldate,ncore)
{
Bmax=max(maize.model2(mat.param[id,],weather,sdate=sdate,ldate=ldate)$B,na.rm=TRUE)
# return output variables for process id
return (c(id=id,Bmax=Bmax))
}

# test of maize_cluster for id=1
maize_cluster(1,sdate,ldate,ncore)

# 2) creating the Cluster: definition of computational units mobilized.
# Example for the definition of four units.
ncore <- 4
cl <-makeCluster (ncore)

# 3) define a calculation sequence.
sequence <- mat.param[,"id"]

# 4) dxport global variables to parallel PU.
# For each calculation unit has good access to global variables
# including function, data defined in the main process data
clusterExport (cl, list ("maize.model2","mat.param", "weather"))

# 5) Allocation of the calculation sequence calculation units.
# The same R function maize_cluster is applied to the cluster "c1"
# where id is picked from sequence
# supplementary arguments (sdate,ldate,ncore) are passed to maize_cluster too
time_with_parallel = system.time(list.result <- clusterApply (cl, sequence, maize_cluster, sdate,ldate,ncore))[3]

# 6) Closing cluster
stopCluster(cl)

# 7) Aggregate results from the parallel computations (a list) to the appropriate structure
mat.result_paral = matrix(unlist(list.result),byrow=TRUE,ncol=2)

# compare to a classical version
time_without_parallel = system.time(mat.result_nonparal <- maize.simule(mat.param[,2:8], weather, sdate, ldate, all=FALSE))[3]
all(mat.result_nonparal==mat.result_paral[,2])

barplot(c(time_with_parallel, time_without_parallel), names.arg = c("with parallel", "without_parallel"), ylab="time (s)")

################################################################################
# performance test run on FBRUN PC (2012-08-23)
# test on a HP Z400 with Intel Xeon CPU W3550 3.07GHz and 4.00Go of Ram
#system.time(mat.result_nonparal <- maize.simule(mat.param[,2:8], weather, sdate, ldate, all=FALSE))["elapsed"]
tps_cal=as.data.frame(rbind(
c(1,10,0.11 ,0.74,1.35 ,1.91,2.39,NA, 3.63 ,NA,NA),
c(1,100,1.05,1.59,1.79,2.16,2.59 ,NA,3.73,NA,NA),
c(1,1000,8.44,7.82,5.05,NA,5.27,NA,6.71,NA,NA),
c(1,5000,44.79,31.36, 17.33,NA, 15.59,NA,14.4,NA,NA),
c(1,10000,102.38,94.24,33.96 ,NA,30.63 ,NA,22.71,NA,NA),
c(1,15000,133.35,116.94 ,45.6 ,NA,37.85,NA,30.69,NA,NA),
c(1,20000, 189.87 , 178.44 ,61.64,46.31,49.58 ,47.45,38.23,34.88,41.04)
))
names(tps_cal) <- c("ncore","n","nonP","P1","P2","P3","P4","P5","P6","P8","P10")

plot(tps_cal$n,tps_cal$nonP ,type="l",col="black",xlab="number of simulation (n)",
ylab="duration of computation (s)",lty=2, main="performance (Intel CPUW3550)")
text(18000,175,"no parallelization", cex=0.9)
lines(tps_cal$n,tps_cal$P2,col="black",lwd=2)
text(18000,64,"ncore=2", cex=0.9)
lines(tps_cal$n,tps_cal$P4,col="black",lwd=3)
text(18000,52,"ncore=4", cex=0.9)
lines(tps_cal$n,tps_cal$P6,col="black",lwd=4)
text(18000,42,"ncore=6", cex=0.9)
lines(tps_cal$n,tps_cal$P8,col="black",lwd=4)
text(18000,20,"ncore=8", cex=0.9)

# end of file