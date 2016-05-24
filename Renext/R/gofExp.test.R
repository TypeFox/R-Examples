gofExp.test <- function(x) {

  nx <- length(x)
  df <- nx - 1
  Tn <- sum(x)
  Mn <- mean(log(x))
  Bn <- 2*nx*(log(Tn/nx) - Mn)/(1 + (nx + 1)/(6*nx) )
  FBn <- pchisq(Bn, df = df)
  
  p.value <- 2*min(c(FBn, 1-FBn))

  #"cat("goodness-of-fit Bartlett'test for exponential\n")
  ##cat("pvalue", p.value, "\n")
  res <- list(statistic = Bn,
              df = df,
              p.value = p.value,
              method = paste("Bartlett gof for exponential"))
  
  invisible(res)
} 
