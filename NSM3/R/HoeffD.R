HoeffD<-function (x, y, example=FALSE) 
{
  # This will calculate Hoeffding's statistic D.
  # Follows section 8.6 of
  #
  #   Nonparametric Statistical Methods, 3e
  #   Hollander, Wolfe & Chicken 
  #
  # Uses the correction for ties given at (8.92).
  #
  # It is intended for small sample sizes n only.  For large n,
  # use the asymptotic equivalence of D to the Blum-Kliefer-Rosenblatt
  # statistic in the R package "Hmisc", command "hoeffd".
  #
  # Very inefficiently programmed by Eric Chicken, October 2012.
  
  if(example)
  {
    x <- c(7.1, 7.1, 7.2, 8.3, 9.4, 10.5, 11.4)
    y <- c(2.8, 2.9, 2.8, 2.6, 3.5, 4.6, 5.0)
  }
  
  n <- length(x)
  
  # phi* (8.93)
  phi <- function(a, b)
  {
    if(a < b) ans <- 1
    if(a == b) ans <- 1 / 2
    if(a > b) ans <- 0
    ans
  }
  
  # c.i (8.92)
  c.i <- numeric(0)
  for(i in 1:n)
  {
    cc.i <- numeric(0)
    for(j in 1:n) if(j !=i) 
      cc.i <- c(cc.i, phi(x[j], x[i]) * phi(y[j], y[i]))
    c.i <- c(c.i, sum(cc.i))
  }
  
  R.i <- rank(x)
  S.i <- rank(y)
  Q <- sum((R.i - 1) * (R.i - 2) * (S.i - 1) * (S.i - 2))
  R <- sum((R.i - 2) * (S.i - 2) * c.i)
  S <- sum(c.i * (c.i - 1))
  D <- Q - 2 * (n - 2) * R + (n - 2) * (n - 3) * S
  D / (n * (n - 1) * (n - 2) * (n - 3) * (n - 4))
}
