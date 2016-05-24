predict.bayesQR <- function(object,X,burnin=0,...){

  # Error handling
  pandterm = function(message) {
    stop(message, call. = FALSE)
  }

  if (!all(sapply(object, "[[", "method") %in% c("QRb","QRb.AL"))) {
	  pandterm("object is not based on 'bayesQR' method with binary dependent variable")
	}

	nqr <- length(object)
  if (nqr < 9){
    pandterm("object should contain the results of at least 9 estimated quantiles") 
  }

  if (!identical(sapply(object,"[[","quantile"),sort(sapply(object,"[[","quantile")))){
	  pandterm("The elements of object are not sorted in increasing order of 'quantile'")
	}

  if(ncol(object[[1]]$betadraw) != ncol(X)){
    pandterm("object has different number of estimated parameters than X has predictors") 
  }

  # Make prediction of latent utility based on the Bayes estimate 
	X <- as.matrix(X)
  n <- nrow(X)
  nvar <- ncol(X)
	outsum <- summary(object=object, burnin=burnin)
  bayesest <- sapply(outsum,"[[","betadraw")[1:nvar,]
	preds <- X%*%bayesest
	
  # Find interval that contains zero
  preds <- cbind(0,preds)
  preds <- t(apply(preds,FUN=sort,MARGIN=1))
  preds <- t(apply(preds,FUN="==",MARGIN=1,0))
  preds <- apply(preds,FUN=which,MARGIN=1)
  
  # Link utilities to probabilities
	pvec <- sapply(object,"[[","quantile")
  pvec <- (c(0,pvec)+c(pvec,1))/2

  # Pr(y=1) = 1 - p
  return(1-pvec[preds])
}
