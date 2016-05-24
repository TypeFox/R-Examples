
permutation <- function(pheno,gt,n.permu,process=TRUE,fileName="",t.thres=1.5)  {
#this version permute some Part samples (e.g. 5%) at a time
#when permutation on cluster, the usage of tMatFunc needs to use different file names #dec 18 13

  gt0         <- gt
  N           <- ncol(gt)
  wls.score.permu <- NULL
  for(i.permu in 1:n.permu){
    if(process){
    	cat( paste("Permutation ",i.permu,"/",n.permu," \n",sep="") )
    	flush.console()
    }
    startTime <- proc.time()[3];
   
  fileName.i        <- paste(fileName, i.permu, sep="")
	#manually select 1% samples
	nMisSample        <- ceiling(0.01 * N)
	if(	nMisSample < 2) nMisSample  <- 2 
	ind.wls0          <- sample(ncol(gt),nMisSample)
	 
	#add wls to genotype
	gt 				        <-  gt0  #add here Dec 4 2012
  elements          <-  unique(c(gt))[!is.na(unique(c(gt)))]  #revsie here May 27 2013
	gt.permu          <-  matrix(sample(elements,nMisSample*nrow(gt),replace =TRUE),nrow=nrow(gt))
	gt[,ind.wls0]     <-  gt.permu
   
    
    t.mat0  <- tMatFunction(pheno,gt,fileName.i) 
    #for permutaion : only use permuted gt as input here(compute wls score for onlyy permuted samples
    #NOTE: t.mat0 must input here , it must be computed based on all genotype data
    delta.t.mat.allmk.list  <- .deltaTMatAllmk_list(gt ,pheno,t.mat0,indSample=ind.wls0,t.thres=t.thres,fileName.i)#a list with length=nMK, each element is a matrix nSample by nSigGene ,      
    deltaT.sampleByall <- NULL   #row:samples, all: all markers for sig genes at each marker
    for(k in 1:nMisSample){
      deltaT.sampleByall <-  rbind(deltaT.sampleByall ,unlist(lapply(delta.t.mat.allmk.list,function(x) x[,k])) )
    }     
    #calculat ethe area of deltaT>0
    wls.score      <- apply( deltaT.sampleByall,1,.area.my)                   #only take the score of swapped samples
    wls.score.permu <- rbind(wls.score.permu,wls.score)
    if(process){
    	newTime <- proc.time()[3];
	cat("   done in ", newTime-startTime, " seconds\n");
	flush.console()
    }
  }    
  return(wls.score.permu)
}
