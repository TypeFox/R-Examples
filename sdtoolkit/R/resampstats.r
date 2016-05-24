

#function to compute distribution on statistics of a box given resampling of dataset

#source("ofun.r")

resampstats <- function(ntimes=100, resampsize=.5,dset, y=NULL, abox, boxtype="OLD"){

  statmat <- matrix(nrow=ntimes,ncol=2)
  dsetsize <- nrow(dset)
  
  if (resampsize <= 1){
    rnum <- resampsize*nrow(dset)
  } else{rnum <- resampsize}  
  
  if(!is.null(y)){dset <- cbind(dset,y)}
   
  for (i in 1:ntimes){
  
    tempind <- sample(dsetsize,rnum)
  
    tempdset <- dset[tempind,]
  
    if (boxtype=="OLD"){
      svect <- beval(tempdset,list(abox))
    } else {stop("Currently only set up to use 'old' box formulation")}
    
    statmat[i,1] <- svect[2] #coverage
    statmat[i,2] <- svect[1] #density
    
  }
 
  covv <- as.numeric(statmat[,1])
  denv <- as.numeric(statmat[,2])
  
  statmat <- cbind(covv,denv)
 
  descrip <- apply(statmat,2,summary)
  
  sdboth <- apply(statmat,2,sd)
  
  descrip <- rbind(descrip,sdboth,(as.numeric(beval(dset,list(abox))[2:1])))
  
  
  gmean <- mean(dset[,ncol(dset)])
  mdiff <-  descrip[4,2] - gmean   
  pseudot <- mdiff/sdboth[2]
  
  rownames(descrip) <- c(rownames(descrip)[1:6],"stdev","full")
  colnames(descrip) <- c("coverage","density")
  attr(descrip,"globalmean") <- gmean
  attr(descrip,"pseudot") <- pseudot
  
  return(descrip)
  
}


  