gpr <- function (h,rec3.sd,rec3.delta,ssd,sdelta,responses) {
  #*********************************************************************
  #  This function uses GPR to find next slope in the flow field to be
  #  used in forecasting.  Historical covariance is determined by using
  # the squared exponential covariance function.
  #
  #  Input: skeleton - Matrix containing the data skeleton
  #         hist.sd - History space of SDCs from psr function
  #         hist.delta - History space of forward response derivatives
  #                      from psr function
  #         rec3.sd - Most recent 3 SDCs
  #         rec3.delta - Most recent 3 forward response derivatives
  #
  #  Output: kmat and kvec from GPR to aid in forecasting
  #
  #  References: 1. C. E. Rasmussen and C. K. I. Williams, Gaussian Processes
  #                 for Machine Learning, Cambridge, MA, MIT Press, 2006
  #
  #              2. Frey, MR and Caudle, KA “Flow field forecasting for 
  #                 univariate time series,” Statistical Analysis and Data 
  #                 Mining, 2013    
  #
  #**********************************************************************
  IDR <- quantile(h[,3],0.9) - quantile(h[,3],0.1)
  corr.l <- (0.5)*IDR/length(h[,1])^(1/3)
  #corr.l <- 0.5
  
  ls <- corr.l*ssd
  ld <- corr.l*sdelta
  
  
  tau2 <- var(responses)
  
  kmat <- matrix(data=0,nrow=length(h[,1]),ncol=length(h[,1]))
  kvec <- matrix(data=0,nrow=length(h[,1]),ncol=1)
  
  for (i in 1:length(h[,1])){
    for (j in 1:length(h[,1])){
      rs <- (h[i,3]-h[j,3])^2+(h[i,2]-h[j,2])^2+(h[i,1]-h[j,1])^2
      rd <- (h[i,6]-h[j,6])^2+(h[i,5]-h[j,5])^2+(h[i,4]-h[j,4])^2;
      kmat[i,j] <- tau2*exp(-rs/(2*ls*ls))*exp(-rd/(2*ld*ld))
    }
  }
  for (i in 1:length(h[,1])){
    rs <- (rec3.sd[3]-h[i,3])^2 + (rec3.sd[2]-h[i,2])^2 + (rec3.sd[1]-h[i,1])^2
    rd <- (rec3.delta[3]-h[i,6])^2 + (rec3.delta[2]-h[i,5])^2 + (rec3.delta[1]-h[i,4])^2
    kvec[i,1] <- tau2*exp(-rs/(2*ls^2))*exp(-rd/(2*ld^2))
  }  
  return(list(first=kmat, second=kvec))
}
