DIME.classify <-
function(data, obj, obj.cutoff = 0.1, obj.sigma.diff.cutoff = NULL,
  obj.mu.diff.cutoff=NULL)
{
  if(!is.null(obj$best)){
    obj$best <- DIME.classify(data, obj$best, obj.cutoff = obj.cutoff, 
      obj.sigma.diff.cutoff = obj.sigma.diff.cutoff, 
      obj.mu.diff.cutoff = obj.mu.diff.cutoff);
    obj$inudge <- DIME.classify(data, obj$inudge, obj.cutoff = obj.cutoff,
      obj.sigma.diff.cutoff = obj.sigma.diff.cutoff,
      obj.mu.diff.cutoff = obj.mu.diff.cutoff);
    obj$gng <- DIME.classify(data, obj$gng, obj.cutoff = obj.cutoff,
      obj.sigma.diff.cutoff = obj.sigma.diff.cutoff,
      obj.mu.diff.cutoff = obj.mu.diff.cutoff);
    obj$nudge <- DIME.classify(data, obj$nudge, obj.cutoff = obj.cutoff,
      obj.sigma.diff.cutoff = obj.sigma.diff.cutoff,
      obj.mu.diff.cutoff= obj.mu.diff.cutoff);
  }
  else{
    data <- unlist(data);
    fsigma.diff.in <- TRUE;
    fmu.diff.in <- TRUE;
    # calculate interquartile region
    iqr <- abs((quantile(data,3/4) - quantile(data,1/4)));
    
    if(is.null(obj.sigma.diff.cutoff)){
      obj.sigma.diff.cutoff <- NULL;
      fsigma.diff.in <- FALSE;
      }
    if(is.null(obj.mu.diff.cutoff)){ 
      obj.mu.diff.cutoff <- iqr/0.6745;
      fmu.diff.in <- FALSE;
      }   
    
    #####################################
    ## Classification for GNG
    if(obj$name=="GNG"){
    idxRem <- NULL;
    nonDiffIdx <- NULL;
    diffPiIdx <- NULL;
    # set all the normal as non-differential initially
    nonDiffIdx <- (1:(obj$K));
    
    ## Categorize component based on mean ##
    for (i in 1:length(nonDiffIdx)){
      if(!fmu.diff.in){
        mco <- obj.mu.diff.cutoff;
      }else if(i <= length(obj.mu.diff.cutoff)){
        mco <- obj.mu.diff.cutoff[i];
      }else if(i > length(obj.mu.diff.cutoff)){
        mco <- obj.mu.diff.cutoff[length(obj.mu.diff.cutoff)];
      }
    	if (abs(obj$mu[nonDiffIdx[i]]) > mco)
      {
    	# save index of object to be removed from differential group
    		if (!is.null(idxRem)){
    			idxRem <- c(idxRem,i);
    		}else{
    			idxRem <-c(i);
    		}
    	}
    }
    # remove index from differential group
    if (!is.null(idxRem) && length(nonDiffIdx) >= length(idxRem)){
    	nonDiffIdx <- nonDiffIdx[-idxRem];
      # find diffIdx and nonDiffIdx in terms of pi index
      nonDiffPiIdx <- (nonDiffIdx+1);
    }
    if(length(nonDiffIdx) == 0){
      cat("Warning: no non-differential component. \n Automatically set",
        " the component with smallest absolute mean as non-differential component\n");
      nonDiffIdx <- which(abs(obj$mu)==min(abs(obj$mu)));
    }
    ## Categorize components based on standard deviation ##
    idxRem <- NULL;
    for (i in 1:length(nonDiffIdx)){
      if(!fsigma.diff.in){ # sigma.diff.cutoff is NOT inputed by user
        sco <- (1.5*abs(iqr)-abs(obj$mu[nonDiffIdx[i]]))/2;
        obj.sigma.diff.cutoff <- c(obj.sigma.diff.cutoff,sco);
      }else if(i <= length(obj.sigma.diff.cutoff)){
        # sigma.diff.cutoff input is less than number of non-diff components
        sco <- obj.sigma.diff.cutoff[i];
      }else if(i > length(obj.sigma.diff.cutoff)){
        sco <- obj.sigma.diff.cutoff[length(obj.sigma.diff.cutoff)];
      }
    	if(abs(obj$sigma[nonDiffIdx[i]]) > sco)
      {
    	# save index of object to be removed from Non-differential group
    		if (!is.null(idxRem)){
    			idxRem <- c(idxRem,i);
    		}else{
    			idxRem <-c(i);
    		}
    	}
    }
    # remove index from Non-differential group
    if (!is.null(idxRem) && length(nonDiffIdx) >= length(idxRem)){
    	nonDiffIdx <- nonDiffIdx[-idxRem];
    }
    if(length(nonDiffIdx) == 0){
      cat("Warning: no non-differential component. \n Automatically set",
        " the component with smallest sigma as non-differential component\n");
      nonDiffIdx <- which(obj$sigma==min(obj$sigma));
    }
    # find diffIdx and nonDiffIdx in terms of pi index
    nonDiffPiIdx <- (nonDiffIdx+1);
    diffIdx <- setdiff(1:obj$K,nonDiffIdx);
    
    # added the exponential components into the differential Pi index
    diffPiIdx <- c((diffIdx+1),1,(obj$K+2));
    obj$diffPiIdx <- diffPiIdx;
    
    # calculate component density of non differentials
    numNonDiff <- length(nonDiffIdx);
    # calculating the estimated model for all probes
    f <- rowSums(obj$phi[,1:(obj$K+2)]);
    if (numNonDiff == 1){
    	f0_psi <- obj$phi[,nonDiffPiIdx];
    } else{
    	f0_psi <- rowSums(obj$phi[,nonDiffPiIdx[1:numNonDiff]]);
    }
    # calculating the local fdr
    if (length(nonDiffPiIdx) > 1) { # no rowsums is needed diffPiIdx <=1
    	obj$fdr <- f0_psi/(f*sum(obj$pi[nonDiffPiIdx]));
    } else {
    	obj$fdr <- f0_psi/(f*obj$pi[nonDiffPiIdx]);
    }
  
    n <- length(obj$fdr);
    obj$class <- rep(0,n);
    obj$class[obj$fdr <= obj.cutoff] <- 1;
  }
  
  ###################################################
  ## Classification for iNudge
  if(obj$name=="iNUDGE"){
    idxRem <- NULL;
    nonDiffIdx <- NULL;
    diffPiIdx <- NULL;
    # set all the normal as non-differential initially
    nonDiffIdx <- (1:(obj$K));
    
    ## Categorize component based on mean ##
    for (i in 1:length(nonDiffIdx)){
      if(!fmu.diff.in){
        mco <- obj.mu.diff.cutoff;
      }else if(i<=length(obj.mu.diff.cutoff)){
        mco <- obj.mu.diff.cutoff[i];
      }else if(i>length(obj.mu.diff.cutoff)){
        mco <- obj.mu.diff.cutoff[length(obj.mu.diff.cutoff)];
      }
    	if (abs(obj$mu[nonDiffIdx[i]]) > mco)
      {
    	# save index of object to be removed from non-differential group
    		if (!is.null(idxRem)){
    			idxRem <- c(idxRem,i);
    		}else{
    			idxRem <-c(i);
    		}
    	}
    }
    # remove index from non-differential group
    if (!is.null(idxRem) && length(nonDiffIdx) >= length(idxRem)){
    	nonDiffIdx <- nonDiffIdx[-idxRem];
    }
    if(length(nonDiffIdx) == 0){
      cat("Warning: no non-differential component. \n Automatically set",
        " the component with smallest absolute mean as non-differential component\n");
      nonDiffIdx <- which(abs(obj$mu)==min(abs(obj$mu)));
    }
    
    ## Categorize component based on standard deviation ##
    idxRem <- NULL;
    for (i in 1:length(nonDiffIdx)){
      if(!fsigma.diff.in){ # sigma.diff.cutoff is NOT inputed by user
        sco <- (1.5*abs(iqr)-abs(obj$mu[nonDiffIdx[i]]))/2;
        obj.sigma.diff.cutoff <- c(obj.sigma.diff.cutoff,sco);
      }else if(i <= length(obj.sigma.diff.cutoff)){
        # sigma.diff.cutoff input is less than number of non-differential normal
        sco <- obj.sigma.diff.cutoff[i];
      }else if(i > length(obj.sigma.diff.cutoff)){
        sco <- obj.sigma.diff.cutoff[length(obj.sigma.diff.cutoff)];
      }
#      cat(fsigma.diff.in," ",sco,"\n");
    	if(abs(obj$sigma[nonDiffIdx[i]]) > sco)
      {
    	# save index of object to be removed from Non-differential group
    		if (!is.null(idxRem)){
    			idxRem <- c(idxRem,i);
    		}else{
    			idxRem <-c(i);
    		}
    	}
    }
    # remove index from Non-differential group
    if (!is.null(idxRem) && length(nonDiffIdx) >= length(idxRem)){
    	nonDiffIdx <- nonDiffIdx[-idxRem];
    }
    if(length(nonDiffIdx) == 0){
      cat("Warning: no non-differential component. \n Automatically set",
        " the component with smallest sigma as non-differential component\n");
      nonDiffIdx <- which(obj$sigma==min(obj$sigma));
    }
    # also add the uniform component
    diffIdx <- setdiff(1:(obj$K),nonDiffIdx)
    obj$diffPiIdx <- c(1,diffIdx+1);
    
    # calculate component density of non differentials
    numNonDiff <- length(nonDiffIdx);
    # calculating the estimated model for all probes
    f <- rowSums(obj$phi[,1:(obj$K+1)]);
    if (numNonDiff == 1){
    	f0_psi <- obj$phi[,nonDiffIdx];
    } else{
    	f0_psi <- rowSums(obj$phi[,nonDiffIdx[1:numNonDiff]]);
    }
    # calculating the local fdr
    if (length(nonDiffIdx) > 1) { # no rowsums is needed diffPiIdx <=1
    	obj$fdr <- f0_psi/(f*sum(obj$pi[nonDiffIdx+1]));
    } else {
      obj$fdr <- f0_psi/(f*obj$pi[nonDiffIdx+1]);
    }
    n <- length(obj$fdr);
    obj$class <- rep(0,n);
    obj$class[obj$fdr <= obj.cutoff]<-1;
    }
    
    ######################################
    ###### Classification for NUDGE #####
    if(obj$name =="NUDGE"){
      n <- length(obj$fdr);
      obj$class <- rep(0,n);
      thresh <- (1-obj.cutoff);
      obj$diffPiIdx <- c(1);
      obj$class <- rep(0,n,1);
      obj$class[obj$pdiff >= thresh] <- 1;
      obj$class[obj$pdiff < thresh] <- 0;
      obj$sigma.diff.cutoff <- NULL;
      obj$mu.diff.cutoff <- NULL;
    }
    
    # save cutoff
    obj$mu.diff.cutoff <- obj.mu.diff.cutoff;
    obj$sigma.diff.cutoff <- obj.sigma.diff.cutoff;
    if(obj$name !="NUDGE"){
      tmp <- which(obj$sigma.diff.cutoff != Inf);
      id <- setdiff(1:obj$K,tmp);
      obj$sigma.diff.cutoff[id] <- Inf;
    }
  }
  return (obj);
}