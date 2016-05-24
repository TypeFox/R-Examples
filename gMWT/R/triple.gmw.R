# Changes:
# 30-06-2013, Added the keepPM option, DF

triple.gmw <- function(X,g,goi,type,nper,alternative,mc,PARAMETERS,output,alg, keepPM,order){

 res <- list()
 pm <- NULL
 if(keepPM) pm <- matrix(NA,ncol=ncol(X),nrow=nper)
 diffTests <- getComb(goi,"triple",order=order)

 METHOD <- c("********* Triple based Test *********")
 DNAME <- PARAMETERS[[1]]
 TEST  <- PARAMETERS[[2]]
 TYPE  <- PARAMETERS[[3]]
 ALTERNATIVE <- PARAMETERS[[4]]
 STATISTIC   <- PARAMETERS[[5]]
 PVAL        <- PARAMETERS[[6]]

 dimX      <- PARAMETERS[[7]]
 XisVector <- PARAMETERS[[8]]

## Case: X is vector
    if(XisVector){
##---------------------------------------------------------------------------------------------------------------------------------------
       if(alternative=="two.sided"){
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, two sided, X is vector
	    for(testRun in 1:nrow(diffTests))
	    {
	      obsValue1 <- as.numeric(getP.Cnaive(X[g==diffTests[testRun,1]],X[g==diffTests[testRun,2]],X[g==diffTests[testRun,3]]))
	      obsValue2 <- as.numeric(getP.Cnaive(X[g==diffTests[testRun,3]],X[g==diffTests[testRun,2]],X[g==diffTests[testRun,1]]))
	      nullDist1 <- perm.triple(X[g==diffTests[testRun,1]],X[g==diffTests[testRun,2]],X[g==diffTests[testRun,3]],nper,algorithm=alg)
	      PVAL <- min(sum(nullDist1>=obsValue1)/nper,sum(nullDist1>=obsValue2)/nper)
	
	      names(PVAL) <- "p.value"
	      STATISTIC <- max(obsValue1,obsValue2)
	      names(STATISTIC) <- "obs.value"
	      ALTERNATIVE <- "two.sided"
	      resTemp<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
	      class(resTemp)<-"htest"
	      
              res[[testRun]] <- resTemp
	      names(res)[testRun] <- paste("H1: P",diffTests[testRun,1],diffTests[testRun,2],diffTests[testRun,3]," > 1/6 or P",diffTests[testRun,3],diffTests[testRun,2],diffTests[testRun,1]," > 1/6",sep="")
	    }
	    if(output=="min")
	    {
	      resMin <- matrix(NA,ncol=1,nrow=length(res))
	      colnames(resMin) <- "pValues"
	      rownames(resMin) <- names(res)
	      for(i in 1:length(res))
	      {
		resMin[i,1] <- res[[i]]$p.value
	      }
	      res <- resMin
	    }
	  } else if(type=="asymptotic"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: asymptotic, two sided, X is vector

	    res <- c()
            stop("We do not have a asymptotic two-sided version for the triple test, sorry!!!")
          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, two sided, X is vector
	    res <- c()
	    stop("We do not have this kind of type for the triple test!")
	  }
##---------------------------------------------------------------------------------------------------------------------------------------
       } else if(alternative=="greater"){
	    if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, greater, X is vector
	    for(testRun in 1:nrow(diffTests))
	    {
	      obsValue <- as.numeric(getP.Cnaive(X[g==diffTests[testRun,1]],X[g==diffTests[testRun,2]],X[g==diffTests[testRun,3]]))
	      nullDist <- perm.triple(X[g==diffTests[testRun,1]],X[g==diffTests[testRun,2]],X[g==diffTests[testRun,3]],nper,algorithm=alg)
	      PVAL <- sum(nullDist>=obsValue)/nper
	
	      names(PVAL) <- "p.value"
	      STATISTIC <- obsValue
	      names(STATISTIC) <- "obs.value"
	      ALTERNATIVE <- "greater"
	      resTemp<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
	      class(resTemp)<-"htest"
	      
              res[[testRun]] <- resTemp
	      names(res)[testRun] <- paste("H1: P",diffTests[testRun,1],diffTests[testRun,2],diffTests[testRun,3]," > 1/6",sep="")
	    }
	    if(output=="min")
	    {
	      resMin <- matrix(NA,ncol=1,nrow=length(res))
	      colnames(resMin) <- "pValues"
	      rownames(resMin) <- names(res)
	      for(i in 1:length(res))
	      {
		resMin[i,1] <- res[[i]]$p.value
	      }
	      res <- resMin
	    }

	  } else if(type=="asymptotic"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: asymptotic, greater, X is vector
	    res <- c()
            stop("We do not have a asymptotic greater version for the triple test, sorry!!!")
          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, greater, X is vector
	    res <- c()
	    stop("We do not have this kind of type for the triple test!")
	  }
       } else if(alternative=="smaller"){
##---------------------------------------------------------------------------------------------------------------------------------------
	   if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, smaller, X is vector
	    for(testRun in 1:nrow(diffTests))
	    {
	      obsValue <- as.numeric(getP.Cnaive(X[g==diffTests[testRun,1]],X[g==diffTests[testRun,2]],X[g==diffTests[testRun,3]]))
	      nullDist <- perm.triple(X[g==diffTests[testRun,1]],X[g==diffTests[testRun,2]],X[g==diffTests[testRun,3]],nper,algorithm=alg)
	      PVAL <- sum(nullDist<obsValue)/nper
	
	      names(PVAL) <- "p.value"
	      STATISTIC <- obsValue
	      names(STATISTIC) <- "obs.value"
	      ALTERNATIVE <- "smaller"
	      resTemp<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
	      class(resTemp)<-"htest"
	      
              res[[testRun]] <- resTemp
	      names(res)[testRun] <- paste("H1: P",diffTests[testRun,1],diffTests[testRun,2],diffTests[testRun,3]," < 1/6",sep="")
	    }
	    if(output=="min")
	    {
	      resMin <- matrix(NA,ncol=1,nrow=length(res))
	      colnames(resMin) <- "pValues"
	      rownames(resMin) <- names(res)
	      for(i in 1:length(res))
	      {
		resMin[i,1] <- res[[i]]$p.value
	      }
	      res <- resMin
	    }

	  } else if(type=="asymptotic"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: asymptotic, greater, X is vector
	    res <- c()
            stop("We do not have a two-sided version for the triple test, sorry!!!")
          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, one sided, X is vector
	    res <- c()
	    stop("We do not have this kind of type for the triple test!")
	  }
       } else {
	    res <- c()
	    stop("There is no other option than small, greater or two-sided...")
       }
## Case: X is a matrix
    } else{
##----------------------------------------------------------------------------------------------------------------------------------------
#Preparational things for the case that X is a matrix
    # First, restrict the cores to maximum of possible tests
    if(mc>detectCores()){
	mc <- detectCores()
	warning("You do not have so many cores on this machine! I automatically reduced it to your maximum number ",mc)
    }
    mc <- min(dimX[2],mc)

    if(alternative=="two.sided"){
       	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, two sided, X is matrix
	    innerLoop <- function(i,testRun){
             nullDist1 <- perm.triple(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i],nper,algorithm=alg)
             obsValue1 <- as.numeric(getP.Cnaive(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i]))
             obsValue2 <- as.numeric(getP.Cnaive(X[g==diffTests[testRun,3],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,1],i]))
             pValue <- min(sum(nullDist1>=obsValue1)/nper,sum(nullDist1>=obsValue2)/nper)
             return(list(pValue=pValue,obsValue=max(obsValue1,obsValue2)))
            }

	    innerLoopPM <- function(i,testRun){
             nullDist1 <- perm.triple(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i],nper,algorithm=alg)
             obsValue1 <- as.numeric(getP.Cnaive(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i]))
             obsValue2 <- as.numeric(getP.Cnaive(X[g==diffTests[testRun,3],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,1],i]))
             pValue <- min(sum(nullDist1>=obsValue1)/nper,sum(nullDist1>=obsValue2)/nper)
             return(list(pValue=pValue,obsValue=max(obsValue1,obsValue2), nullDist=nullDist1))
            }

	    if(keepPM){
	        nullDistRES <- list()
		STATISTIC <- list()
		for(i in 1:nrow(diffTests)){
		  nullDistRES[[i]] <- matrix(0, ncol=dimX[2],nrow=nper)
		  STATISTIC[[i]] <- c(rep(-1,dimX[2]))
		}
	    }	 

	    for(testRun in 1:nrow(diffTests))
	    { 
	      resTemp <- list()

	      if(keepPM==TRUE){
   	        resInner <-  unlist(mclapply(c(1:dimX[2]),innerLoopPM,testRun=testRun,mc.cores=mc))
		#nullDistRES <- matrix(0, ncol=dimX[2],nrow=nper)
              } else {
   	        resInner <- unlist(mclapply(c(1:dimX[2]),innerLoop,testRun,mc.cores=mc))
              }

	      for(i in 1:dimX[2])
	      {
		if(keepPM==TRUE){
                  PVAL <- resInner[nper*(i-1) + 2*(i) - 1]
                  STATISTIC[[testRun]][i] <- resInner[nper*(i-1) + 2*i]
                  nullDistRES[[testRun]][,i] <- resInner[(nper*(i-1) + 2*i + 1):(nper*i + 2*i)]
                } else {
		  PVAL <- resInner[2*i-1]
                  STATISTIC <- resInner[2*i]
		}
		obsValue <- STATISTIC
		names(PVAL) <- "p.value"
		ALTERNATIVE <- "two.sided"
		names(STATISTIC) <- "obs.value"
		resTemp[[i]]<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
		class(resTemp[[i]])<-"htest"	    
	      }
	     res[[testRun]] <- resTemp
	     names(res)[testRun] <- paste("H1: P",diffTests[testRun,1],diffTests[testRun,2],diffTests[testRun,3]," > 1/6 or P",diffTests[testRun,3],diffTests[testRun,2],diffTests[testRun,1]," > 1/6",sep="")
	    }
	    if(output=="min")
	    {
	      resMin <- matrix(NA,ncol=dimX[2],nrow=length(res))
	      colnames(resMin) <- colnames(X)
	      rownames(resMin) <- names(res)
	      for(i in 1:length(res))
	      {
		for(j in 1:dimX[2])
		{
		  resMin[i,j] <- res[[i]][[j]]$p.value
		}
	      }
	      res <- resMin
	    }
	  } else if(type=="asymptotic"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: asymptotic, two sided, X is matrix
	    res <- c()
	    stop("We do not have a two-sided version for the triple test, sorry!!!,A,T,M")
          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, two sided, X is vector
	    res <- c()
	    stop("We do not have this kind of type for the triple test!,O,T,M")
	  }
    } else if(alternative=="greater"){
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, greater, X is matrix
	      # Define the function, that is performed for column i (important for parallelization)
	   innerLoop <- function(i,testRun){
             nullDist <- perm.triple(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i],nper,algorithm=alg)
             obsValue <- as.numeric(getP.Cnaive(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i]))
             pValue <- sum(nullDist>=obsValue)/nper
	     return(list(pValue=pValue,obsValue=obsValue))
            }

	   innerLoopPM <- function(i,testRun){
             nullDist <- perm.triple(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i],nper,algorithm=alg)
             obsValue <- as.numeric(getP.Cnaive(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i]))
             pValue <- sum(nullDist>=obsValue)/nper
	     return(list(pValue=pValue,obsValue=obsValue, nullDist=nullDist))
            }

	    if(keepPM){
	        nullDistRES <- list()
		STATISTIC <- list()
		for(i in 1:nrow(diffTests)){
		  nullDistRES[[i]] <- matrix(0, ncol=dimX[2],nrow=nper)
		  STATISTIC[[i]] <- c(rep(-1,dimX[2]))
		}
	    }

	    for(testRun in 1:nrow(diffTests))
	    { 
	      resTemp <- list()

	      if(keepPM==TRUE){
   	        resInner <-  unlist(mclapply(c(1:dimX[2]),innerLoopPM,testRun=testRun,mc.cores=mc))
		#nullDistRES <- matrix(0, ncol=dimX[2],nrow=nper)
              } else {
   	        resInner <-  unlist(mclapply(c(1:dimX[2]),innerLoop,testRun=testRun,mc.cores=mc))
              }

	      for(i in 1:dimX[2])
	      {

		if(keepPM==TRUE){
                  PVAL <- resInner[nper*(i-1) + 2*(i) - 1]
                  STATISTIC[[testRun]][i] <- resInner[nper*(i-1) + 2*i]
                  nullDistRES[[testRun]][,i] <- resInner[(nper*(i-1) + 2*i + 1):(nper*i + 2*i)]
                } else {
		  PVAL <- resInner[2*i-1]
                  STATISTIC <- resInner[2*i]
		}
		obsValue <- STATISTIC
		names(PVAL) <- "p.value"
		ALTERNATIVE <- "greater"
		names(STATISTIC) <- "obs.value"
		resTemp[[i]]<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
		class(resTemp[[i]])<-"htest"	    
	      }
	     res[[testRun]] <- resTemp
	     names(res)[testRun] <- paste("H1: P",diffTests[testRun,1],diffTests[testRun,2],diffTests[testRun,3]," > 1/6",sep="")
	    }
	    if(output=="min")
	    {
	      resMin <- matrix(NA,ncol=dimX[2],nrow=length(res))
	      colnames(resMin) <- colnames(X)
	      rownames(resMin) <- names(res)
	      for(i in 1:length(res))
	      {
		for(j in 1:dimX[2])
		{
		  resMin[i,j] <- res[[i]][[j]]$p.value
		}
	      }
	      res <- resMin
	    }
	  } else if(type=="asymptotic"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: asymptotic, greater, X is matrix
	    res <- c()
            stop("We do not have a two-sided version for the triple test, sorry!!!A,2S,V")
	  } else {
	    res <- c()
	    stop("We do not have this kind of type for the UIT!,O,G,M")
	  }
    } else if(alternative=="smaller"){
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, smaller, X is matrix
	    innerLoop <- function(i,testRun){
             nullDist <- perm.triple(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i],nper,algorithm=alg)
             obsValue <- as.numeric(getP.Cnaive(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i]))
             pValue <- sum(nullDist<obsValue)/nper
	     return(list(pValue=pValue,obsValue=obsValue))
            }

	    innerLoopPM <- function(i,testRun){
             nullDist <- perm.triple(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i],nper,algorithm=alg)
             obsValue <- as.numeric(getP.Cnaive(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i]))
             pValue <- sum(nullDist<obsValue)/nper
	     return(list(pValue=pValue,obsValue=obsValue, nullDist=nullDist))
            }

	    if(keepPM){
	        nullDistRES <- list()
		STATISTIC <- list()
		for(i in 1:nrow(diffTests)){
		  nullDistRES[[i]] <- matrix(0, ncol=dimX[2],nrow=nper)
		  STATISTIC[[i]] <- c(rep(-1,dimX[2]))
		}
	    }

	    for(testRun in 1:nrow(diffTests))
	    { 
	      resTemp <- list()

	      if(keepPM==TRUE){
   	        resInner <-  unlist(mclapply(c(1:dimX[2]),innerLoopPM,testRun=testRun,mc.cores=mc))
		#nullDistRES <- matrix(0, ncol=dimX[2],nrow=nper)
              } else {
   	        resInner <-  unlist(mclapply(c(1:dimX[2]),innerLoop,testRun=testRun,mc.cores=mc))
              }

	      for(i in 1:dimX[2])
	      {

		if(keepPM==TRUE){
                  PVAL <- resInner[nper*(i-1) + 2*(i) - 1]
                  STATISTIC[[testRun]][i] <- resInner[nper*(i-1) + 2*i]
                  nullDistRES[[testRun]][,i] <- resInner[(nper*(i-1) + 2*i + 1):(nper*i + 2*i)]
                } else {
		  PVAL <- resInner[2*i-1]
                  STATISTIC <- resInner[2*i]
		}
		obsValue <- STATISTIC
		names(PVAL) <- "p.value"
		ALTERNATIVE <- "smaller"
		names(STATISTIC) <- "obs.value"
		resTemp[[i]]<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
		class(resTemp[[i]])<-"htest"	    
	      }
	     res[[testRun]] <- resTemp
	     names(res)[testRun] <- paste("H1: P",diffTests[testRun,1],diffTests[testRun,2],diffTests[testRun,3]," < 1/6",sep="")
	    }
	    if(output=="min")
	    {
	      resMin <- matrix(NA,ncol=dimX[2],nrow=length(res))
	      colnames(resMin) <- colnames(X)
	      rownames(resMin) <- names(res)
	      for(i in 1:length(res))
	      {
		for(j in 1:dimX[2])
		{
		  resMin[i,j] <- res[[i]][[j]]$p.value
		}
	      }
	      res <- resMin
	    }
	  } else if(type=="asymptotic"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: asymptotic, smaller, X is matrix
	    res <- c()
	    stop("We do not have this kind of type for the triple test!,A,S,M")
          } else {
	    res <- c()
	    stop("We do not have this kind of type for the triple test!,O,S,M")
	  }
    } else {
	    res <- c()
	    stop("There are no other alternatives possible, sorry! All other....")
     }
  }
  if(type=="permutation"){
    ifelse(keepPM,res <- list(p.values=res, nullDist=nullDistRES, obsValue=obsValue), res <- list(p.values=res))
  } else {
    res <- list(p.values=res)
  }
  res
}