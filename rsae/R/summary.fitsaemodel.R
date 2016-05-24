summary.fitsaemodel <-
function (object, digits=6, ...){
   saemodel <- attr(object, "saemodel")
   n <- saemodel$n
   g <- saemodel$g
   p <- saemodel$p
   fixeffnames <- attr(saemodel, "xnames")
   # retrieve the estimating method
   method <- attr(object, "method") 
   # check whether the model converged
   converged <- object$converged
   if (converged != 1){
      cat("ALGORITHM FAILED TO CONVERGE! use convergence() to learn more \n")
   }else{
      cat("ESTIMATION SUMMARY \n")
      cat("Method: ", method$type, "\n")
      # branch: robust vs. non-robust methods
      if (length(method) > 1){
	 tuning <- method$tuning
	 if (length(tuning) == 1){
	    cat(paste("Robustness tuning constant: ",names(tuning), " = ", as.numeric(tuning), "\n"))
	 }else{
	    for (i in 1:length(tuning)){
	       cat(tuning[i], "\n")  
	    }
	 }
      }
      #----------------------
      # fixed effects table
	 df <- n - g - p + 1
	 fixed <- object$beta
	 stdfixed <- sqrt(diag(object$vcovbeta))
	 tTable <- cbind(fixed, stdfixed, fixed/stdfixed, df, fixed)
	 colnames(tTable) <- c("Value", "Std.Error", "t-value", "df", "p-value")
	 rownames(tTable) <- fixeffnames 
	 tTable[, 5] <- 2 * pt(-abs(tTable[, 3]), tTable[, 4])
	 cat("---\n")
	 cat("Fixed effects\n")
	 printCoefmat(tTable, digits=digits, P.values=TRUE, has.Pvalue=TRUE)
   
      #----------------------
      # robustness properties
      robustness <- attr(object, "robustness") 
      # branch robust vs non-robust
      if (!is.null(robustness)){
	 wgt <- t(robustness$wgt) / n
	 colnames(wgt) <- c("fixeff", "residual var", "area raneff var")
	 rownames(wgt) <- "sum(wgt)/n"
	 cat("---\n")
	 cat("Degree of downweighting/winsorization:\n")
	 cat("\n")
	 print.default(format(t(wgt), digits = digits), print.gap = 2, quote = FALSE)
      }
   }
}

