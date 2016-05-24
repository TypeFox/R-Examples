dparetotol.int <- function(x, m = NULL, alpha=0.05, P=0.99, side = 1, ...){
  if (side != 1 && side != 2) {
    stop(paste("Must specify a one-sided or two-sided procedure!", 
               "\n"))
  }
  if (side == 2) {
    alpha <- alpha/2
    P <- (P + 1)/2
  }
  n <- length(x)
  if(is.null(m)) m <- n
  out <- dpareto.ll(x, ...)
  theta <- as.numeric(stats4::coef(out))
  CI <- theta+c(-1,1)*qnorm(1-alpha)*sqrt(abs(stats4::vcov(out)[1]))*sqrt(n/m)
  CI[1] <- max(CI[1],1e-14)
  CI[2] <- min(CI[2],1)
  lower <- ifelse(CI[1]==0,0,max(qdpareto(1 - P, theta = CI[1]),0))
  upper <- ifelse(CI[2]==1,Inf,qdpareto(P, theta = CI[2]))
  if (side == 2) {
    alpha <- 2 * alpha
    P <- (2 * P) - 1
  }
  temp <- data.frame(cbind(alpha, P, theta, lower, upper))
  if (side == 2) {
    colnames(temp) <- c("alpha", "P", "theta", "2-sided.lower", "2-sided.upper")
  }
  else {
    colnames(temp) <- c("alpha", "P", "theta", "1-sided.lower", "1-sided.upper")
  }
  temp	
}
