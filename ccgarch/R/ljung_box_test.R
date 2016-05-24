ljung.box.test <- function(x){
  if(!is.vector(x)){
      cat("Error: The argument must be a vector.\n")
      # stop("\n")
   } else {
      nobs <- length(x); nlag <- seq(1,50,1)
      y <- x-mean(x)
      denom <- sum(y^2)
      rho <- numeric(length(nlag))
         for(i in 1:length(nlag)){
            rho[i] <- sum(y[(nlag[i]+1):nobs]*y[1:(nobs-nlag[i])])/denom
         }
      rho <- nobs*(nobs+2)*(rho^2/(nobs-nlag))
      Q <- cumsum(rho)
      p.val <- pchisq(Q, df=nlag, lower.tail=FALSE)
      ans <- cbind(Q, p.val); colnames(ans) <- c("test stat", "p-value"); rownames(ans) <- paste("Lag", nlag)
      ans[seq(5,50,5),]
   }
}
