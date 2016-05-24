pHoeff<-function (n=5, reps=10000, r=4) 
{
  # This will approximate the distribution of Hoeffding's statistic D
  # under the null.  This code follows section 8.6 of
  #
  #   Nonparametric Statistical Methods, 3e
  #   Hollander, Wolfe & Chicken 
  #
  # reps is the number of Monte Carlo runs to produce.
  #
  # This calls HoeffD, a small bit of code that produces the value
  # of D without any inference.
  #
  # It is intended for small sample sizes n only.  For large n,
  # use the asymptotic equivalence of D to the Blum-Kliefer-Rosenblatt
  # statistic in the R package "Hmisc", command "hoeffd".
  #
  # Very inefficiently programmed by Eric Chicken, October 2012.
  
  D <- numeric(0)
  for(m in 1:reps)
  {
    # These are independent samples, no ties (with probability 1):
    x <- runif(n)
    y <- runif(n)
    D <- c(D, HoeffD(x, y))
  }
  
  D.values <- sort(unique(D))
  D.prob <- D.values
  for(i in 1:length(D.values)) 
    D.prob[i] <- length(D[D == D.values[i]]) / reps
  D <- cbind(D.values, D.prob, cumsum(D.prob), rev(cumsum(rev(D.prob))))
  rownames(D) <- rep("", length(D.values))
  colnames(D) <- c("d", "P(D = d)", "P(D <= d)", "P(D >= d)")
  round(D, r)
}
