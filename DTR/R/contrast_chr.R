###################################################
### code chunk number 1: chunklibraries
###################################################

#Libraries required
require(aod)

###################################################
### code chunk number 2: chunkfunction
###################################################

contrast_chr <- function(est,
                          t=quantile(est$time, 0.75) # A time of interest for comparison, optional, default is the 75th percentile of time
) {
  
  if (t >= max(est$time)) stop("Time beyond the scope of the data")
  if (t <= min(est$time)) stop("Time beyond the scope of the data")
  
  #Round the time
  t <- round(t,2)
  
  #Find the index for estiamtes
  index <- which(abs(est$time-t)==min(abs(est$time-t)[which(est$time<=t)]))
    
  #Combine survival estimates as column vector
  EST <- matrix(c(est$CHR1211.LOG[index], est$CHR2111.LOG[index], 
                   est$CHR2211.LOG[index], est$CHR2112.LOG[index],
                   est$CHR2212.LOG[index], est$CHR2221.LOG[index]), nrow=6, ncol=1)
  #Combine variance estimates as matrix
  VAR <- matrix(c(est$SE1211.LOG[index]^2, est$COV1211_2111.LOG[index], 
                  est$COV1211_2211.LOG[index], est$COV1211_2112.LOG[index],
                  est$COV1211_2212.LOG[index], est$COV1211_2221.LOG[index],
                  est$COV1211_2111.LOG[index], est$SE2111.LOG[index]^2,
                  est$COV2111_2211.LOG[index], est$COV2111_2112.LOG[index], 
                  est$COV2111_2212.LOG[index], est$COV2111_2221.LOG[index],
                  est$COV1211_2211.LOG[index], est$COV2111_2211.LOG[index], 
                  est$SE2211.LOG[index]^2, est$COV2211_2112.LOG[index],
                  est$COV2211_2212.LOG[index], est$COV2211_2221.LOG[index],
                  est$COV1211_2112.LOG[index], est$COV2111_2112.LOG[index],
                  est$COV2211_2112.LOG[index], est$SE2112.LOG[index]^2,
                  est$COV2112_2212.LOG[index], est$COV2112_2221.LOG[index],
                  est$COV1211_2212.LOG[index], est$COV2111_2212.LOG[index],
                  est$COV2211_2212.LOG[index], est$COV2112_2212.LOG[index],
                  est$SE2212.LOG[index]^2, est$COV2212_2221.LOG[index],
                  est$COV1211_2221.LOG[index], est$COV2111_2221.LOG[index],
                  est$COV2211_2221.LOG[index], est$COV2112_2221.LOG[index],
                  est$COV2212_2221.LOG[index], est$SE2221.LOG[index]^2), nrow=6, ncol=6)

  #Test for H0: A1B1=A1B2=A2B1=A2B2
  #log(theta1211)=log(theta2111)=log(theta2211)=0
  test_overall <- wald.test(b = EST, Sigma = VAR, 
                            L=matrix(c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0), nrow=3, ncol=6, byrow=T))
  
  #Test for H0: A1B1=A1B2
  #log(theta1211)=0
  test_1112 <- wald.test(b = EST, Sigma = VAR, 
                         L=matrix(c(1,0,0,0,0,0), nrow=1, ncol=6, byrow=T))
  
  #Test for H0: A1B1=A2B1
  #log(theta2111)=0
  test_1121 <- wald.test(b = EST, Sigma = VAR, 
                         L=matrix(c(0,1,0,0,0,0), nrow=1, ncol=6, byrow=T))
  
  #Test for H0: A1B1=A2B2
  #log(theta2211)=0
  test_1122 <- wald.test(b = EST, Sigma = VAR, 
                         L=matrix(c(0,0,1,0,0,0), nrow=1, ncol=6, byrow=T))
  
  #Test for H0: A1B2=A2B1
  #log(theta2112)=0
  test_1221 <- wald.test(b = EST, Sigma = VAR, 
                         L=matrix(c(0,0,0,1,0,0), nrow=1, ncol=6, byrow=T))
  
  #Test for H0: A1B2=A2B2
  #log(theta2212)=0
  test_1222 <- wald.test(b = EST, Sigma = VAR, 
                         L=matrix(c(0,0,0,0,1,0), nrow=1, ncol=6, byrow=T))
  
  #Test for H0: A2B1=A2B2
  #log(theta2221)=0
  test_2122 <- wald.test(b = EST, Sigma = VAR, 
                         L=matrix(c(0,0,0,0,0,1), nrow=1, ncol=6, byrow=T))

  #Combine results
  results <- t(data.frame(test_overall$result, test_1112$result, test_1121$result, 
                          test_1122$result, test_1221$result, test_1222$result, test_2122$result))
  rownames(results) <- NULL

  #Return results
  TEST <- data.frame(c("A1B1=A1B2=A2B1=A2B2", "A1B1=A1B2", "A1B1=A2B1", "A1B1=A2B2", 
                       "A1B2=A2B1", "A1B2=A2B2", "A2B1=A2B2"),
                     results)
  names(TEST) <- c(paste("H0 (t=", t, ")", sep=""), "test statistic", "df", "p")
  return(TEST)    

}  