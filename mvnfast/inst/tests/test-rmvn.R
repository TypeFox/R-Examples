context("rmvn()")

test_that("Checking rmvn()", {
  
  ##########
  ###### d = 1, n = 1 case
  ##########
  N <- c(1e4, 1e4, 1e4)
  d <- c(1,     2,   3) 
  
  message("Testing rmvn()")
  for(ii in 1:length(N))
  {
    mu <- 1:d[ii]
    X <- t(t(matrix(rnorm(N[ii]*d[ii]), N[ii], d[ii])) + mu)
    tmp <- matrix(rnorm(d[ii]^2), d[ii], d[ii])
    mcov <- tcrossprod(tmp, tmp)
    myChol <- chol(mcov)
    
    tolMu <- 0.1 * sum(abs(mu))
    tolCov <- 0.1 * sum(sum(abs(mcov)))
        
    ####### Sequential
    # Using covariance
    X <- rmvn(N[ii], mu, mcov)
    expect_less_than(sum(sum(abs(colMeans(X) - mu))), tolMu)
    expect_less_than(sum(sum(abs(cov(X) - mcov))), tolCov)
    # Using Cholesky
    X <- rmvn(N[ii], mu, myChol, isChol = TRUE)
    expect_less_than(sum(sum(abs(colMeans(X) - mu))), tolMu)
    expect_less_than(sum(sum(abs(cov(X) - mcov))), tolCov)
    
    ####### Parallel
    # Using covariance
    X <- rmvn(N[ii], mu, mcov, ncores = 2)
    expect_less_than(sum(sum(abs(colMeans(X) - mu))), tolMu)
    expect_less_than(sum(sum(abs(cov(X) - mcov))), tolCov)
    # Using Cholesky
    X <- rmvn(N[ii], mu, myChol, ncores = 2, isChol = TRUE)
    expect_less_than(sum(sum(abs(colMeans(X) - mu))), tolMu)
    expect_less_than(sum(sum(abs(cov(X) - mcov))), tolCov)
 
    message(paste("Test", ii, "passed."))
  }
  
})


