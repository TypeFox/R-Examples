epi.prcc <- function(dat, sided.test = 2){
   # Calculate mu and number of parameters:
   N <- dim(dat)[1]
   K <- dim(dat)[2] - 1
   # Return an error message if the number of parameters is greater than the number of model replications:
   if(K > N) 
     stop("Error: the number of replications of the model must be greater than the number of model parameters")
   
   mu <- (1 + N) / 2
   # Compute ranks:
   for(i in 1:(K + 1)){
      dat[,i] <- rank(dat[,i])
   }

   # Create a K + 1 by K + 1 matrix:
   C <- matrix(rep(0, times = (K + 1)^2), nrow = (K + 1))

   # Correlations for each parameter pair:
   for(i in 1:(K + 1)){
    for(j in 1:(K + 1)){
      r.it <- dat[,i]
      r.jt <- dat[,j]
      r.js <- r.jt
      c.ij <- sum((r.it - mu) * (r.jt - mu)) / sqrt(sum((r.it - mu)^2) * sum((r.js - mu)^2))
      C[i,j] <- c.ij
    }
   }

  # Fudge to deal with response variables that are all the same:
  if(is.na(C[K + 1,K + 1]))
     {gamma.iy <- rep(0, times = K)
      t.iy <- gamma.iy * sqrt((N - 2) / (1 - gamma.iy^2))
      p <- rep(1, times = K)
      df <- rep((N - 2), times = K)
      # Results:
      rval <- as.data.frame(cbind(gamma = gamma.iy, test.statistic = t.iy, df = df, p.value = p)) 
      return(rval)
      }
   
   else { 
      # Matrix B is defined as the inverse of c.ij:
      B <- solve(C)
      # PRCC (gamma.iy) between the ith input parameter and the yth outcome is defined by Kendall and Stewart (1979) as follows:
      gamma.iy <- c()
      for(i in 1:K){
         num <- -B[i,(K + 1)]
         den <- sqrt(B[i,i] * B[(K + 1),(K + 1)])
         gamma.iy <- c(gamma.iy, num/den)
      }

      # Email Andrew Hill (mailto:fyu7@cdc.gov) 14 August 2009. 
      # I think there may be an error in the epi.prcc function. Looking at the example in the package documentation, I note that the p-values all close to 0 yet if we switch the sign of y to force negative correlation with the x's we get p-values near 1. Backtracking, I think the problem is in the definition of the test statistic. There is a typo in the Blower-Dowlatabadi paper. I believe they misstate the test statistic at the end of the Appendix. It should be (dropping their subscripts): t = gamma * sqrt((N-2) / (1 - gamma^2)).
      # Equivalently, the square of the PRCC gamma^2 is asymptotically Beta(1/2,(N-2)/2). (ref. Muirhead's book 'Aspects of Multivariate Statistical Theory').

      # Blower and Dowlatabadi appendix: 
      # t.iy <- gamma.iy * sqrt((N - 2) / (1 - gamma.iy))
      # Andrew Hill's correction:
      t.iy <- gamma.iy * sqrt((N - 2) / (1 - gamma.iy^2))
      df <- rep((N - 2), times = K)      

      # Blower and Dowlatabadi appendix: 
      # p <- 1 - pt(q = t.iy, df = (N - 2))
      
      if(sided.test == 2){
         # Andrew Hill's correction:
         p <- 2 * pt(abs(t.iy), df = N - 2, lower.tail = FALSE)
         # p <- pbeta(t.iy^2, shape1 = 1/2, shape2 = (N - 2)/2, lower.tail = FALSE)
         }

      if(sided.test == 1){
         # Andrew Hill's correction:
         p <- pt(abs(t.iy), df = N - 2,lower.tail = FALSE)
         }
         
      # Results:
      rval <- as.data.frame(cbind(gamma = gamma.iy, test.statistic = t.iy, df = df, p.value = p)) 
      return(rval)
      }
}