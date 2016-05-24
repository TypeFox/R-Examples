mctr.Bayes <- function(tickers, weights = rep(1,length(tickers)),
                       start, end, data, sim.size = 1000){
  p <- length(tickers)
  if(length(weights)!=p)
    stop("unequal lengths of the arguments 'tickers' & 'weights'", call. = FALSE)
  
  if(p>1)
  {
    dat <- access(tickers, start, end, data)
    # Check if sufficient data (at least 2 values needed) is available for all the tickers to compute variance-covariance matrix
    chk <- colSums(is.na(dat)) < (nrow(dat)-1)
    # If not, delete columns with insufficient data
    dat <- dat[,chk]
    if(sum(chk)!=p)
      warning("insufficient data available for the ticker(s) ", paste(tickers[!chk], collapse=", "), "; at least 2 values for each of the tickers should be available for the given time period", call. = FALSE)
    tickers <- tickers[chk]
    p <- length(tickers)
  }
  
  n <- nrow(dat)
  S <- cov(dat, use="pairwise.complete.obs")
  
  if(p >= n){
    stop("insufficient data; large p, small n")
    #     n0 <- (p-n)+3
    #     df <- n0+(n-1)
    #     Psi <-  diag(1,nrow=nrow(S))*max(diag(S))
  } else{
    df <- n
    Psi <- S
  }
  
  w <- weights/sum(weights)
  mctr.sim <- matrix(NA, nrow = sim.size, ncol = p)
  for(i in 1:sim.size){
    Sigma <- riwish(v = df, S = Psi)
    Sigma <- Sigma*(df-p-1)
    sigma_p <- sqrt((t(w) %*% Sigma %*% w)[1,1])
    mctr.sim[i,] <- ((Sigma %*% w)/sigma_p)[,1]
  }
  
  colnames(mctr.sim) <- tickers
  return(100*mctr.sim)
}