# Version: 03-07-2013

# Changes:
# 28-06-2013: Added the option keepPM, DF
# 30-06-2013: Finished the keepPM option, DF
# 03-07-2013: Adjusted the warning/stop messages, DF

uit.gmw <- function(X,g,goi,type,nper,alternative,mc,PARAMETERS,output, keepPM){

 res <- list()
 diffTests <- getComb(goi, "triple", order=T)

 METHOD <- c("********* Union-Intersection Test *********")
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

## Alternative: two-sided
##---------------------------------------------------------------------------------------------------------------------------------------
       if(alternative=="two.sided"){
	  stop("There is no 2-sided alternative available for UIT, please choose 'greater'!")
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, two sided, X is vector
	    res <- c()
	    stop("There is no permutation type test implemented for the two-sided UIT!")
	  } else if(type=="asymptotic"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: asymptotic, two sided, X is vector
	    res <- c()
	    stop("There is no asymptotic test implemented for the two-sided UIT!")
          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, two sided, X is vector
	    res <- c()
	    stop("There is no such test type implemented for the two-sided UIT!")
	  }
## Alternative: greater
##---------------------------------------------------------------------------------------------------------------------------------------
       } else if(alternative=="greater"){
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, greater, X is vector
	    for(testRun in 1:nrow(diffTests))
	    {
	      obsValue <- uit.C(X[g==diffTests[testRun,1]],X[g==diffTests[testRun,2]],X[g==diffTests[testRun,3]])
	      nullDist <- uitPTest(X[g==diffTests[testRun,1]],X[g==diffTests[testRun,2]],X[g==diffTests[testRun,3]],nper)
	      PVAL <- sum(obsValue<=nullDist)/nper
	
	      DNAME <- paste("Data:",deparse(substitute(X)),", Groups:",deparse(substitute(g)),", Order: max(P",diffTests[testRun,1],diffTests[testRun,3],",P",diffTests[testRun,2],diffTests[testRun,3],")",sep="")
	      names(PVAL) <- "p.value"
	      STATISTIC <- obsValue
	      names(STATISTIC) <- "obs.value"
	      ALTERNATIVE <- "greater"
	      resTemp<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
	      class(resTemp)<-"htest"
	      res[[testRun]] <- resTemp
	      names(res)[testRun] <- paste("H1: Max(P",diffTests[testRun,1],diffTests[testRun,3],",P",diffTests[testRun,2],diffTests[testRun,3],") > 0.5",sep="")
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
	    for(testRun in 1:nrow(diffTests))
	    {
	      obsValue <- uit.C(X[g==diffTests[testRun,1]],X[g==diffTests[testRun,2]],X[g==diffTests[testRun,3]])
	      PVAL <- uitATest(X[g==diffTests[testRun,1]],X[g==diffTests[testRun,2]],X[g==diffTests[testRun,3]],obsValue)
	
	      DNAME <- paste("Data:",deparse(substitute(X)),", Groups:",deparse(substitute(g)),", Order: max(P",diffTests[testRun,1],diffTests[testRun,3],",P",diffTests[testRun,2],diffTests[testRun,3],")",sep="")
	      names(PVAL) <- "p.value"
	      STATISTIC <- obsValue
	      names(STATISTIC) <- "obs.value"
	      ALTERNATIVE <- "greater"
	      resTemp<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
	      class(resTemp)<-"htest"
	      res[[testRun]] <- resTemp
	      names(res)[testRun] <- paste("H1: Max(P",diffTests[testRun,1],diffTests[testRun,3],",P",diffTests[testRun,2],diffTests[testRun,3],") > 0.5",sep="")
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
          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, greater, X is vector
	    res <- c()
	    stop("We do not have this kind of type for the UIT!")
	  }
       } else if(alternative=="smaller"){
##---------------------------------------------------------------------------------------------------------------------------------------
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, smaller, X is vector
	    res <- c()
	    stop("There is no smaller alternative available, please choose greater!")
	  } else if(type=="asymptotic"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: asymptotic, smaller, X is vector
	    res <- c()
            stop("There is no smaller alternative available, please choose greater!")
          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, one sided, X is vector
	    res <- c()
	    stop("There is no smaller alternative available, please choose greater!")
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

	    res <- c()
	    stop("There is no 2-sided alternative available, please choose greater!")

	  } else if(type=="asymptotic"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: asymptotic, two sided, X is matrix
	    res <- c()
            stop("There is no 2-sided alternative available, please choose greater!")
          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, two sided, X is vector
	    res <- c()
	    stop("There is no 2-sided alternative available, please choose greater!")
	  }
    } else if(alternative=="greater"){
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, greater, X is matrix
	      # Define the function, that is performed for column i (important for parallelization)
	    innerLoop <- function(i,testRun){
	      nullDist <- uitPTest(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i],nper)
	      obsValue <- uit.C(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i])
	      pValue <- sum(obsValue<=nullDist)/nper
	      return(list(pValue=pValue,obsValue=obsValue))
            }

	    innerLoopPM <- function(i,testRun){
	      nullDist <- uitPTest(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i],nper)
	      obsValue <- uit.C(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i])
	      pValue <- sum(obsValue<=nullDist)/nper
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
		ALTERNATIVE <- "greater"
		DNAME <- paste("Data:",deparse(substitute(X)),", Groups:",deparse(substitute(g)),", Order: max(P",diffTests[testRun,1],diffTests[testRun,3],",P",diffTests[testRun,2],diffTests[testRun,3],")",sep="")
		names(STATISTIC) <- "obs.value"
		resTemp[[i]]<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
		class(resTemp[[i]])<-"htest"	    
	      }
	     res[[testRun]] <- resTemp
	     names(res)[testRun] <- paste("H1: Max(P",diffTests[testRun,1],diffTests[testRun,3],",P",diffTests[testRun,2],diffTests[testRun,3],") > 0.5",sep="")
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
	    innerLoop <- function(i,testRun){
	      obsValue <- uit.C(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i])
	      pValue <- uitATest(X[g==diffTests[testRun,1],i],X[g==diffTests[testRun,2],i],X[g==diffTests[testRun,3],i],obsValue)
	      return(list(pValue=pValue,obsValue=obsValue))
            }
	    
	    for(testRun in 1:nrow(diffTests))
	    {
	      resTemp <- list()
	      resInner <- unlist(mclapply(c(1:dimX[2]),innerLoop,testRun,mc.cores=mc))
	      for(i in 1:dimX[2])
	      {
		PVAL <- resInner[2*i-1]
		STATISTIC <- resInner[2*i]
		names(PVAL) <- "p.value"
		ALTERNATIVE <- "greater"
		DNAME <- paste("Data:",deparse(substitute(X)),", Groups:",deparse(substitute(g)),", Order: max(P",diffTests[testRun,1],diffTests[testRun,3],",P",diffTests[testRun,2],diffTests[testRun,3],")",sep="")
		names(STATISTIC) <- "obs.value"
		resTemp[[i]]<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
		class(resTemp[[i]])<-"htest"	    
	      }
	     res[[testRun]] <- resTemp
	     names(res)[testRun] <- paste("H1: Max(P",diffTests[testRun,1],diffTests[testRun,3],",P",diffTests[testRun,2],diffTests[testRun,3],") > 0.5",sep="")
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
          } else {
	    res <- c()
	    warning("We do not have this kind of type for the UIT!")
	  }
    } else if(alternative=="smaller"){
     #  res <- do.call(rbind,mclapply(c(1:dimX[2]),innerLoop,mc.cores=mc))
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, smaller, X is matrix
	    res <- c()
	    stop("There is no smaller alternative available, please choose greater!")
	  } else if(type=="asymptotic"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: asymptotic, smaller, X is matrix

	    res <- c()
	    stop("We do not have this kind of type for the UIT!")
          } else {
	    res <- c()
	    stop("We do not have this kind of type for the UIT!")
	  }
    } else {
	    res <- c()
	    stop("There are no other alternatives possible, sorry!")
     }
  }
  if(type=="permutation"){
    ifelse(keepPM,res <- list(p.values=res, nullDist=nullDistRES, obsValue=obsValue), res <- list(p.values=res))
  } else {
    res <- list(p.values=res)
  }
  res
}