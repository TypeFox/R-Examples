zelen.test<-function (z, example=F, r=3) 
{
  # Zelen's test.
  # Based on chapter 10 of:
  #
  #   Nonparametric Statistical Methods, 3e
  #   Hollander, Wolfe & Chicken 
  #
  # The data z is an array of k 2x2 matrices.
  # Small data sets only!
  # Uses package "partitions".
  #
  # Inefficiently programmed by Eric Chicken, November 2012.
  
 
  if(example) z <- array(c(2, 1, 2, 5, 1, 5, 4, 1), dim=c(2, 2, 2))
  
  s <- sum(z[1, 1, ])
  k <- dim(z)[3]
  
  # blockparts is from package "partitions".  This is where large data
  # sets will be an issue.
  # Make sure that each part of the sum is no more than the column or
  # row margin total.
  bp <- numeric(0)
  for(i in 1:k) bp <- c(bp, min(sum(z[1,,i]),sum(z[,1,i])))
  a <- blockparts(bp, s)
  
  
  y <- numeric(0)
  for(i in 1:dim(a)[2])
  {
    is.tau.0 <- T
    x <- numeric(0)
    for(j in 1:k)
    {
      O.11 <- a[j, i]
      O.12 <- sum(z[1, , j]) - O.11
      O.21 <- sum(z[, 1, j]) - O.11
      O.22 <- sum(z[2, , j]) - O.21
      tau <- matrix(c(O.11, O.12, O.21, O.22), nrow=2, byrow=T)
      if(sum(tau == z[, , j]) < 4) is.tau.0 <- F
      n1 <- O.11 + O.12
      n2 <- O.21 + O.22
      n.1 <- O.11 + O.21
      n <- n1 + n2
      x.j <- choose(n1, O.11) * choose(n2, O.21) / choose(n, n.1)
      x <- c(x, x.j)
    }
    if(is.tau.0) tau.0 <- i
    y <- c(y, prod(x))
  }
  y <- y / sum(y)
  p <- sum(y[y<=y[tau.0]])
  
  cat("\n")
  cat("Zelen's test:") 
  cat("\n")
  cat(paste("P = ", round(p, r), sep=""))
  cat("\n")
}
