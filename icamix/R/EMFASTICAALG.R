###########################
## EMFASTICAALG.R
##
## Written by Xiaotian Zhu
###########################

##
## accessory functions
##

save_rng <- function(savefile="tempfile1"){# function for saving radom seed info
  oldseed <- oldRNGkind <- NULL
  rm(oldseed)
  rm(oldRNGkind)
  if (exists(".Random.seed"))  {
    oldseed <- get(".Random.seed", .GlobalEnv)
  } else stop("don't know how to save before set.seed() or r*** call")
  oldRNGkind <- RNGkind()
  save("oldseed","oldRNGkind",file=savefile)
  invisible(savefile)
}

restore_rng <- function(savefile) {# function for loading radom seed info
  oldseed <- oldRNGkind <- NULL
  rm(oldseed)
  rm(oldRNGkind)
  load(savefile)
  do.call("RNGkind",as.list(oldRNGkind))  
  assign(".Random.seed", oldseed, .GlobalEnv)
}

ESTIMATEDMEMBER <-
  function(rstICAMIX){# function calculates estimated class membership from an EMFASTICAALG object
    postmembership <- rstICAMIX$MembershipProbs == apply(rstICAMIX$MembershipProbs, 1, max) # estimated membership matrix
    ncluster <- ncol(rstICAMIX$MembershipProbs)
    estimatedmember <- factor(rep((1:ncluster),ncol(rstICAMIX$InputData))[as.vector(t(postmembership))]) # estimated membership info
    estimatedmember
  }

CLASSDIFFRATE <-
  function(factor1, factor2){# function calculates classification difference rate between two factors
    length(which(as.numeric(levels(factor1))[factor1] != as.numeric(levels(factor2))[factor2]))/length(factor1)
  }

ATRANSDENSITY <- function(grid, A, f){
  # function evaluates density of linearly transformed random vector
  # on a given grid of which the columns store the grid points
  n <- ncol(grid)
  grid_ <- solve(A) %*% grid
  answer <- f(grid_)
  answer <- answer / abs(det(A))
  answer
}


##
## main function
##
## an R wrapper for carrying out NSMM-ICA
## on nonparametric multivariate ICA mixture data
## 
##
EMFASTICAALG <-
  function (DataMatrix, numCluster, h = 0, 
            maxiter = 300, icaiter = 150, tol = 1e-6, 
            verb = TRUE, combine = TRUE, seednum = 82196, ...){
    RNGstore <- save_rng() # get current random seed info
    
    t0 <- proc.time() # get current time
    
    n <- nrow(DataMatrix) # number of observation
    r <- ncol(DataMatrix) # data dimension
    
    #if(h == 0){
    #  h = 0.9 * (n*r)^{-0.2} * min( sd(DataMatrix), IQR(DataMatrix)/1.34 )
    #}
    
    MemberProbs <- matrix(0, n, numCluster)
    set.seed(seednum)
    preCluster <- kmeans(DataMatrix, numCluster)
    for(i in 1:numCluster)
    {
      MemberProbs[preCluster$cluster==i, i] <- 1
    }
    res <- .Call('icamix_EMInterwovenFastICA', PACKAGE = 'icamix', t(DataMatrix), MemberProbs, h, maxiter+1, icaiter, tol, verb, combine)
    for (i in 1:numCluster){
      res$WMtrs[[i]] <- matrix(res$WMtrs[[i]],r)
      res$WUnmixZ[[i]] <- matrix(res$WUnmixZ[[i]],r)
      res$OriginalSignals[[i]] <- matrix(res$OriginalSignals[[i]],r)
    }
    #for (i in 1:r){
    #  tempresult$Densities[[i]] <- matrix(tempresult$Densities[[i]], numCluster)
    #}
    t1 <- proc.time()
    computingtime <- (t1 - t0)[3]
    cat("\n computing time:", computingtime, " s\n")
    restore_rng(RNGstore)
    file.remove(RNGstore)
    res$call <- match.call()
    res$time <- computingtime
    class(res) <- "EMFASTICAALG"
    res
  }


##
## methods functions
##
print.EMFASTICAALG <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nNumber of Observations: ")
  cat(ncol(x$InputData))
  cat("\nAttributes per Observation: ")
  cat(nrow(x$InputData))
  cat("\nNumber of Mixture Components: ")
  cat(ncol(x$Lambdas))
  cat("\nNumber of Iterations: ")
  cat(nrow(x$Lambdas)-1)
  cat("\nComputing Time Elapsed (sec): ")
  cat(as.numeric(x$time))
  cat("\nData Likelihood from Last Interation:\n")
  cat(x$ObjValue[nrow(x$ObjValue),1])
  cat("\nEstimated Lambdas:\n")
  cat(x$Lambdas[nrow(x$Lambdas),])
  invisible(x)
}

summary.EMFASTICAALG <- function(object, ...)
{
  # extract simple info from EMFASTICAALG result
  n <- ncol(object$InputData)
  r <- nrow(object$InputData)
  m <- ncol(object$Lambdas)
  t <- nrow(object$Lambdas)-1
  obj <- object$ObjValue[nrow(object$ObjValue),1]
  lambdas <- object$Lambdas[nrow(object$Lambdas),]
  estimatedclass <- ESTIMATEDMEMBER(object) # estimated class info
  datalikelihood <- object$ObjValue[-1,]
  
  # mixture component means and variances
  compmeans <- list()
  compvars <- list()
  for(j in 1:m){
    compmeans[[j]] <- object$InputData %*% object$MembershipProbs[,j] / sum(object$MembershipProbs[,j])
    centeredData <- object$InputData - compmeans[[j]] %*% rep(1,n)
    compvars[[j]] <- centeredData %*% diag(object$MembershipProbs[,j]) %*% t(centeredData) / ( (sum(object$MembershipProbs[,j]))^2 - sum((object$MembershipProbs[,j])^2) ) * sum(object$MembershipProbs[,j])
  }
  
  # compile the summary list and return
  res <- list(inputData=object$InputData, originSig=object$OriginalSignals, call=object$call, time=as.numeric(object$time), numIter=t, lastObj=obj, compMeans=compmeans, compVars=compvars, numObs=n, numAtr=r, numCls=m, estWts=lambdas, estCls=estimatedclass, dataLklhd=datalikelihood)
  class(res) <- "summary.EMFASTICAALG"
  res
}

print.summary.EMFASTICAALG <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat(x$numObs)
  cat(" Observations, ")
  cat(x$numAtr)
  cat(" Attributes, and ")
  cat(x$numCls)
  cat(" Classes.\n")
  cat("Computation Lasts for ")
  cat(x$time)
  cat(" Seconds and ")
  cat(x$numIter)
  cat(" Iterations.")
  cat("\nData Likelihood from Last Interation:\n")
  cat(x$lastObj)
  cat("\nEstimated Lambdas:\n")
  cat(x$estWts)
  cat("\nHard Classification: ")
  print(table(x$estCls))
  cat("\nMean and Covariance of Mixture Component:\n")
  for(j in 1:(x$numCls)){
    cat(paste("Component #",j,"\n Mean:\n"))
    print(x$compMeans[[j]])
    cat(paste(" Covariance:\n"))
    print(x$compVars[[j]])
  }
  invisible(x)
}

plot.summary.EMFASTICAALG <- function(x, vec1=c(1:2), vec2=c(1:2), ...){
  # plot data likelihood vs interation
  if(length(x$dataLklhd) > 1){
    plot.ts(x$dataLklhd,xlab="iteration",ylab="data loglikelihood",main="Objective Function Values")
  }
  
  # plot data with color according to NSMM-ICA classification
  plot((t(x$inputData))[,vec1], asp = 1, col=as.numeric(levels(x$estCls))[x$estCls]+1, pch = as.numeric(levels(x$estCls))[x$estCls], main = "NSMM-ICA Classification",cex.lab=1.1)
  
  # plot recovered independent signals
  par(mfrow = c(1,x$numCls))
  for(clsind in 1:(x$numCls)){
    keyind <- rank(levels(x$estCls))[clsind]
    plot(t(x$originSig[[keyind]])[x$estCls==clsind,vec2], asp = 1, col = clsind + 1, pch = clsind, xlab = 'X', ylab = 'Y', main = paste("Comp #", clsind," IndSig"))
  }
  par(mfrow=c(1,1)) # reset par
  invisible(x)
}

plot.EMFASTICAALG <- function(x, vec1=c(1:2), vec2=c(1:2), ...){
  plot(summary(x), vec1, vec2)
  invisible(x)
}









