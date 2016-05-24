DeltaC <-
function (Data,order) {
lmResults <- lm(y ~ x1 + x2 + x1*x2,data=Data) # fit lm() 
# extract coefficients
coefficient <- summary(lmResults)$coefficients 
# extract variance covariance matrix
covariance <- vcov(lmResults)
B2 <- coefficient[3,1] # estimation of B2
B3 <- coefficient[4,1] # estimation of B3
COV22 <- covariance[3,3] # variance of B2
COV33 <- covariance[4,4] # variance of B3
k <- 1.96 # 95% percentile for N(0,1)
if (order==1) {
 SE_delta <- sqrt((((COV22)/(B3^2))+(((B2^2)*COV33)/(B3^4))))
}
else if (order==2) {
SE_delta <- sqrt((((COV22)/(B3^2))+(((B2^2)*COV33)/(B3^4)))
      + (COV22*COV33+2*(B2^2)*(COV33)^2)/(B3^4))
}
C_Hat <- (-1)*B2/B3
LowCI <- C_Hat - k*SE_delta
UpperCI <- C_Hat + k*SE_delta
results <- list(LowCI = LowCI, UpperCI = UpperCI)
return(results)
}
