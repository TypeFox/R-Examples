
rm(list=ls())
library (msm)
# The matrix of transition rates (M-matrix)
  M <- cbind(c(  0.15,-0.10,-0.05),        
             c( -0.07, 0.10,-0.03),
             c( -0.02,-0.05, 0.07))
  diag(M) <- diag(M) + 0.00000001                      
  namstates <- c("A","B","C")
  dimnames(M) <- list(destination=namstates,origin=namstates)
# Sample size
 nsample <- 200
# Subject identification numbers and ages at birth
 ID <- 1:nsample
 born <- rep(0,nsample)
# Life segment to be simulated
 Age.min <- rep(20,nsample)
 Age.max <- rep (40,nsample)
# State occupied at onset of simulation
 st_entry <- sample (1:length(namstates),nsample,replace=TRUE)
# Covariate (identical individuals)
 cov1 <- rep("X",nsample)
# Simulation: transitions: dates, states occupied and state sequence (path)
 dates <- array (as.numeric(NA),dim=c(nsample,30),dimnames=list (ID=ID,trans=paste ("Tr",1:30,sep="")))
 states <- array (as.character(NA), dim=c(nsample,30),dimnames=list (ID=ID,trans=paste ("Tr",1:30,sep="")))
 path <- vector (mode="character",length=nsample)
 for (i in 1:nsample)
  { # Simulate life history of an individual
  	bio <- sim.msm (-t(M),mintime=Age.min[i],maxtime=Age.max[i],start=st_entry[i])
  	if (length (bio$times)>2) dates[i,1:(length (bio$times)-2)] <- bio$times[-c(1,length(bio$times))]
  	path[i] <- namstates[bio$states][1]
  	if (length(bio$times)>2) for (j in 2:(length(bio$states)-1)) path[i] <- paste (path[i],namstates[bio$states][j],sep="")
  	states[i,1:length (bio$states)] <- bio$states
 }
# Create Biograph object  
ndatesmax <- max(apply (dates,1,function(x) length(na.omit(x))))
D.sim <- data.frame (ID=ID,born=born,start=Age.min,end=Age.max,cov=cov1,path=path,round(dates[,1:ndatesmax],2),stringsAsFactors = FALSE)
namcov <- c(paste("cov",1:1,sep=""))
colnames(D.sim) <- c("ID","born","start","end",namcov,"path",paste("Tr",1:ndatesmax,sep=""))
attr(D.sim,"format.date") <- "year"
attr(D.sim,"format.born") <- "year"
attr(D.sim,"param") <- Parameters (D.sim)



