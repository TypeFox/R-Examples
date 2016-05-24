sen.adichie<-function (z, example=F, r=3) 
{
  # This code tests for parallel lines.
  # Based on chapter 9 of:
  #
  #   Nonparametric Statistical Methods, 3e
  #   Hollander, Wolfe & Chicken 
  #
  # z is a list of paired vectors.  Each item in the list is a set
  # of two paired vectors in the form of a matrix.  The first column
  # of each matrix is the x vector, the second in the y vector.
  #
  # Inefficiently programmed by Eric Chicken, October 2012.
  
  if(example)
  {
    # Example 9.5 data
    x1 <- x2 <- x3 <- x4 <- c(0, 1.5, 3, 4.5, 6)
    y1 <- c(0, 33.019, 111.314, 196.205, 230.658)
    y2 <- c(0, 131.831, 181.603, 230.07, 258.119)
    y3 <- c(0, 33.351, 97.463, 196.615, 217.308)
    y4 <- c(0, 8.959, 105.384, 211.392, 255.105)
    z <- list(cbind(x1, y1), cbind(x2, y2), cbind(x3, y3), 
              cbind(x4, y4))
  }
  
  k <- length(z)
  x.bar <- beta.bar.num <- beta.bar.den <- numeric(0)
  for(i in 1:k)
  {
    # (9.41)
    x.bar <- c(x.bar, mean(z[[i]][,1]))
    # (9.40)
    beta.bar.num <- c(beta.bar.num, 
                      sum((z[[i]][, 1] - x.bar[i]) * z[[i]][, 2]))
    beta.bar.den <- c(beta.bar.den, 
                      sum((z[[i]][, 1] - x.bar[i])^2))
    
  }
  beta.bar <- sum(beta.bar.num) / sum(beta.bar.den)
  # (9.42) - (9.45)
  
  V <- numeric(0)
  for(i in 1:k)
  {
    z[[i]][, 2] <- z[[i]][, 2] - beta.bar * z[[i]][, 1]
    T.i <- (z[[i]][, 1] - x.bar[i]) * rank(z[[i]][, 2])
    T.i <- sum(T.i) / (length(z[[i]][, 1]) + 1)
    V <- c(V, T.i^2 / beta.bar.den[i])
  }
  V <- 12 * sum(V)
  p <- pchisq(V, df=(k - 1), lower.tail=F)
  
  cat("\n")
  cat("Null: all slopes are equal") 
  cat("\n")
  cat(paste("V = ", round(V, r), ", P = ", round(p, r), sep=""))
  cat("\n")
  
  
}
