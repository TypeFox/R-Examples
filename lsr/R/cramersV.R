# file:    cramersV 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 29 June 2013

# cramersV() calculates the Cramer's V measure of effect size for chi-square tests.
# I haven't thought about this function in a while: it might not be the right way
# to go about this one.
cramersV <- function (...) {
  
  test <- chisq.test(...)
  chi2 <- test$statistic
  N <- sum(test$observed)
  
  if (test$method =="Chi-squared test for given probabilities"){
    # for GOF test, calculate max chi-square value 
    ind <- which.min(test$expected)
    max.dev <- test$expected
    max.dev[ind] <- N-max.dev[ind]
    max.chi2 <- sum( max.dev ^2 / test$expected ) 
    V <- sqrt( chi2 / max.chi2 )
  } 
  else { 
    # for test of association, use analytic expression
    k <- min(dim(test$observed))
    V <- sqrt( chi2 / (N*(k-1)) )
  }
  names(V) <- NULL
  return(V)  
  
}
