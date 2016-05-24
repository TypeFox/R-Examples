print.fitsaemodel <-
function (x, digits=6, ...){
   saemodel <- attr(x, "saemodel")
   # check whether the model converged
   if (x$converged == 1){
      yname <- attr(saemodel, "yname")
      xnames <- attr(saemodel, "xnames")
      areadef <- attr(saemodel, "areadef")
      # retrieve the estimating method
      method <- attr(x, "method") 
      cat("ESTIMATES OF SAE-MODEL (model type B) \n")
      cat("Method: ", method$type, "\n")
      # branch: robust vs. non-robust methods
      if (length(method) > 1){
	 tuning <- method$tuning$k
	 if (length(tuning) == 1){
	    cat(paste("Robustness tuning constant: k = ", tuning, "\n"))
	 }else{
	     cat(paste("Robustness tuning constants: k_beta = ", tuning[1], ", k_v = ", tuning[2], ", k_d = ", tuning[3],"\n", sep=""))
	 }
      }
      cat("---\n")
      cat("Fixed effects\n")
      cat(paste("Model: ", yname, " ~ ", paste(xnames, collapse=" + "), sep=""), "\n")
      cat("  Coefficients: \n")
      beta <- x$beta
      names(beta) <- xnames
      if (saemodel$intercept == 1){
	 names(beta)[1] <- "(Intercept)"
      }
      print.default(format(beta, digits = digits), print.gap = 2, quote = FALSE)
      cat("--- \n")
      cat("Random effects \n")
      cat(paste("  Model: ~1| ", areadef, sep=""), "\n")
      # warn if the raneff variance is almost zero
      if (x$theta[2] <= .Machine$double.eps^(1/4)){
	 cat("---\n")
	 cat("NOTE THAT THE VARIANCE OF THE AREA-LEVEL RANDOM \nEFFECT IS ALMOST ZERO! DO YOU REALLY NEED THE \nRANDOM EFFECT? IF SO, GO AHEAD. HOWEVER, YOU \nSHOULD CONSIDER FITTING A (ROBUST) GLS MODEL.\n")
	 cat("---\n")
      }
      theta <- sqrt(x$theta)
      # change the order of theta
      theta <- c(theta[2], theta[1])
      names(theta) <- c("(Intercept)", "Residual")
      theta <- as.matrix(theta)
      colnames(theta) <- "Std. Dev."
      print.default(format(t(theta), digits = digits), print.gap = 2, quote = FALSE)
      cat("--- \n")
      cat("Number of Observations: ", saemodel$n, "\n")
      cat("Number of Areas: ", saemodel$g, "\n")
      cat("\n")
   }else{
      #not converged
      cat("THE METHOD DID NOT CONVERGE!\n")
      cat("---\n")
      cat("  1) use convergence() of your fitted model to learn more \n")
      cat("  2) study the documentation using the command ?fitsaemodel \n")
      cat("  3) you may call fitsaemodel with 'init' equal to (either) 'lts'\n")
      cat("     or 's' (this works also for ML, though it may no be very efficient)\n")
      cat("  4) if it still does not converge, the last resort is to modify\n")
      cat("     'acc' and/or 'niter' (and hope and pray)\n")
   }
   invisible(x)
}

