summary.bayesQR.single <- function(QRobj, burnin=0, credint=c(.025,.975)){

  # Dimensions of QRobj
  k <- ncol(QRobj$betadraw)
  ndraw <- nrow(QRobj$betadraw)

  # Error handling
  pandterm <- function(message) {
    stop(message, call. = FALSE)
  }

  # Check burnin
  burnin <- burnin+1
  if(burnin>ndraw){
    pandterm("Burnin cannot be larger than number of mcmc draws R") 
  }

  # Check specified credible interval 
  if ((credint[1]>credint[2])|(credint[1]<=0)|(credint[2]>=1)){
    pandterm("Incorrect credible interval")
  }

  # Model information
	out <- list()
  out$method <- QRobj$method
  out$quantile <- QRobj$quantile
  out$names <- QRobj$names
  out$burnin <- burnin-1
  out$retained <- ndraw-out$burnin
  out$credint <- credint

  # Summarize betadraw
  out$betadraw <- cbind(apply(QRobj$betadraw[burnin:ndraw,],FUN="mean",MARGIN=2),
                        t(apply(QRobj$betadraw[burnin:ndraw,],FUN="quantile",MARGIN=2,prob=credint)))
  colnames(out$betadraw) <- c("Bayes Estimate", "lower", "upper") 
  rownames(out$betadraw) <- QRobj$names

  # If present, summarize sigmadraw
  if (QRobj$method %in% c("QRc","QRc.AL")){
    out$sigmadraw <- matrix(c(mean(QRobj$sigmadraw[burnin:ndraw]),
                              quantile(QRobj$sigmadraw[burnin:ndraw],prob=credint)),ncol=3)
    colnames(out$sigmadraw) <- c("Bayes Estimate", "lower", "upper") 
    rownames(out$sigmadraw) <- "sigma" 
  }

  return(out)
}
