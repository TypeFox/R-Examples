print.bayesQR.summary <- function(x, digits=3, ...){

  # Number of quantile regression summaries to print
	nqr <- length(x)

	# Loop trough every quantile regression 
	for (i in 1:nqr){
	QRsub <- x[[i]]
    cat("\n")
    if (QRsub$method=="QRc"){
      cat("Type of dependent variable: continuous\n")
      cat("Lasso variable selection: no\n")
    } else if (QRsub$method=="QRc.AL"){
      cat("Type of dependent variable: continuous\n")
      cat("Lasso variable selection: yes\n")
    } else if (QRsub$method=="QRb"){
      cat("Type of dependent variable: binary\n")
      cat("Lasso variable selection: no\n")
    } else if (QRsub$method=="QRb.AL"){
      cat("Type of dependent variable: binary\n")
      cat("Lasso variable selection: yes\n")
    }
    cat(paste("Estimated quantile: ",QRsub$quantile),"\n")
    cat(paste("Lower credible bound: ",QRsub$credint[1]),"\n")
    cat(paste("Upper credible bound: ",QRsub$credint[2]),"\n")
    cat(paste("Number of burnin draws: ",QRsub$burnin),"\n")
    cat(paste("Number of retained draws: ",QRsub$retained),"\n")
    cat("\n")
    cat("\n")
    cat("Summary of the estimated beta:\n")
    cat("\n")
    print(QRsub$betadraw,digits=digits)
    cat("\n")
    if (QRsub$method %in% c("QRc","QRc.AL")){
      cat("\n")
      cat("Summary of the estimated sigma:\n")
      cat("\n")
      print(QRsub$sigmadraw,digits=digits)
      cat("\n")
    }
  	if ((nqr>1)&(i<nqr)) cat("*****************************************\n")
	}
}
