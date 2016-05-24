
alphaf <- function(nsbp, cc) {
  
  cc <- as.vector(cc)
  times <- nsbp - 1
  chiM <- matrix(1 - pchisq(rep(cc, each = times), 1 : times), ncol = times, nrow = length(cc), 
                 byrow = TRUE)
  as.vector(chiM %*% dbinom(1 : times, times, 0.5))
  
}

alphaV <- c(0.71, 1.73, 2.59, 3.39, 4.14, 4.86, 5.56, 6.75, 6.92, 8.24, 9.53, 10.79,
            12.03, 13.26, 16.28, 19.25, 1.64, 2.95, 4.01, 4.96, 5.84, 6.67, 7.48, 8.26,
            9.02, 10.50, 11.93, 13.33, 14.70, 16.04, 19.34, 22.55, 2.71, 4.23, 5.43,
            6.50, 7.48, 8.41, 9.29, 10.15, 10.99, 12.60, 14.15, 15.66, 17.13, 18.57,
            22.09, 25.50, 3.84, 5.54, 6.86, 8.02, 9.09, 10.09, 11.05, 11.98, 12.87,
            14.57, 16.25, 17.85, 19.41, 20.94, 24.65, 28.23, 9.55, 11.76, 13.47, 14.95,
            16.32, 17.58, 18.78, 19.93, 21.03, 23.13, 25.17, 27.10, 28.96, 30.79,
            35.15, 39.33)
alphaM <- matrix(alphaV, 16, 5)
sbp <- c(2 : 10, 12, 14, 16, 18, 20, 25, 30)
sig <- c(0.20, 0.10, 0.05, 0.025, 0.001)

cvaluef <- function(nsbp, alpha) {
  
  if(any(sbp == nsbp) & any(sig == alpha)) {
    cvalue <- as.numeric(alphaM[which(sbp == nsbp), which(sig == alpha)])
  } else {
    cinit <- as.numeric(alphaM[which.min(abs(sbp - nsbp)), which.min(abs(sig - alpha))])
    if(alphaf(nsbp, cinit) > alpha) {
      repeat{
        cinit <- cinit + 0.01
        if(alphaf(nsbp, cinit) < alpha) {
          cvalue <- cinit
          break
        }
      }
    } else if(alphaf(nsbp, cinit) < alpha) {
      repeat{
        cinit <- cinit - 0.01
        if(alphaf(nsbp, cinit) > alpha) {
          cvalue <- cinit + 0.01
          break
        }
      }
    } 
  }
  
  cvalue
  
}
