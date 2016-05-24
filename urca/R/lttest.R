##
## lttest
##
lttest <- function(z, r){
  if(!(class(z)=="ca.jo")){
    stop("\nObject 'x' must be of class 'ca.jo'\n")
  }
  r <- as.integer(r)
  if(r >= z@P || r < 1){
    stop("\nProvide number of cointegration vectors in valid range.\n")
  }
  idx <- r + 1
  df <- length(idx:z@P)
  N <- nrow(z@Z0)
  test1 <- ca.jo(z@x, ecdet = "const", K=z@lag, season=z@season, dumvar=z@dumvar)
  lambda1 <- test1@lambda
  test2 <- ca.jo(z@x, ecdet = "none", K=z@lag, season=z@season, dumvar=z@dumvar)
  lambda2 <- test2@lambda
  teststat <- -N*sum(log((1-lambda1[idx:z@P])/(1-lambda2[idx:z@P])))
  pval <- 1 - pchisq(teststat, df)
  lttest <- as.matrix(t(c(teststat, pval)))
  colnames(lttest) <- c("test statistic", "p-value")
  rownames(lttest) <- "LR test"
  cat("LR-test for no linear trend\n")
  cat("\n")
  cat(paste("H0: H*2(r<=", r, ")\n", sep=""))
  cat(paste("H1: H2(r<=", r, ")\n", sep=""))
  cat("\n")
  cat("Test statistic is distributed as chi-square\n")
  cat(paste("with", df, "degress of freedom\n", sep=" "))
  print(round(lttest, 2))
}
