###################################################
### code chunk number 1: chunklibraries
###################################################

#Libraries required
require(aod)

###################################################
### code chunk number 2: chunkfunction
###################################################

contrast_wald <- function(est,
                          t=quantile(est$time, 0.75) # A time of interest for comparison, optional, default is the 75th percentile of time
) {
  
  if (t >= max(est$time)) stop("Time beyond the scope of the data")
  if (t <= min(est$time)) stop("Time beyond the scope of the data")
  
  #Find the index for estiamtes
  index <- which(abs(est$time-t)==min(abs(est$time-t)[which(est$time<=t)]))
    
  #Combine survival estimates as column vector
  SURV <- matrix(c(est$SURV11[index], est$SURV12[index], est$SURV21[index], 
                   est$SURV22[index]), nrow=4, ncol=1)
  #Combine variance estimates as matrix
  VAR <- matrix(c(est$SE11[index]^2, est$COV1112[index], 0, 0, est$COV1112[index], 
                  est$SE12[index]^2, 0, 0,
                  0, 0, est$SE21[index]^2, est$COV2122[index], 0, 0, 
                  est$COV2122[index], est$SE22[index]^2), nrow=4, ncol=4)
  
  #Test for H0: S11=S12=S21=S22  
  test_overall <- wald.test(b = SURV, Sigma = VAR, L=matrix(c(1,1,1,-1,0,0,0,-1,0,0,0,-1), nrow=3, ncol=4))
  
  #Test for H0: S11=S12
  test_1112 <- wald.test(b = SURV, Sigma = VAR, L=matrix(c(1,-1,0,0), nrow=1, ncol=4))
  
  #Test for H0: S11=S21
  test_1121 <- wald.test(b = SURV, Sigma = VAR, L=matrix(c(1,0,-1,0), nrow=1, ncol=4))

  #Test for H0: S11=S22
  test_1122 <- wald.test(b = SURV, Sigma = VAR, L=matrix(c(1,0,0,-1), nrow=1, ncol=4))
  
  #Test for H0: S12=S21
  test_1221 <- wald.test(b = SURV, Sigma = VAR, L=matrix(c(0,1,-1,0), nrow=1, ncol=4))
  
  #Test for H0: S12=S22
  test_1222 <- wald.test(b = SURV, Sigma = VAR, L=matrix(c(0,1,0,-1), nrow=1, ncol=4))
  
  #Test for H0: S21=S22
  test_2122 <- wald.test(b = SURV, Sigma = VAR, L=matrix(c(0,0,1,-1), nrow=1, ncol=4))

  #Combine results
  results <- t(data.frame(test_overall$result, test_1112$result, test_1121$result, 
                          test_1122$result, test_1221$result, test_1222$result, test_2122$result))
  rownames(results) <- NULL
  #Return restuls
  TEST <- data.frame(c("A1B1=A1B2=A2B1=A2B2", "A1B1=A1B2", "A1B1=A2B1", "A1B1=A2B2", 
                       "A1B2=A2B1", "A1B2=A2B2", "A2B1=A2B2"),
                     results)
  names(TEST) <- c(paste("H0 (t=", round(t,2), ")", sep=""), "test statistic", "df", "p")
  return(TEST)    
}  