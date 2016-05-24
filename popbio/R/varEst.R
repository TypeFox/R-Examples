varEst <- function (rates,weighted=1)
{
  times<-length(unique(rates[,2]))
  classes<-unique(rates[,1])
  n<-length(classes)
  results <- matrix(rep(0,3*n), nrow=n,
               dimnames= list(classes, c("total.var", "demo.var", "env.var")))  
  for (rate in 1:n)
  {
    minrow   <- (rate - 1) * times + 1
    maxrow   <-  rate * times
    N <- rates[minrow:maxrow, 3]
    m <- rates[minrow:maxrow, 4]
    p <- m / N
    if (weighted == 0)
    {
      demvar <- sum((p*(1 - p))/ N)/ times # equ.2 of Akcakaya 2002: unweighted mean of demographic variance
      totalvar <- var(p, na.rm = TRUE)    # unweighted total variance of the vital rate
    }
    else
    {
      demvar <- sum(p*(1 - p)) / sum(N)   # equ.3 of Akcakaya 2002: weighted average demographic variance
      wmeanp <- sum(m) / sum(N)           # weighted mean of the vital rate
      totalvar <- sum(N * (p - wmeanp)^2) / sum(N) # weigthed total variance (Kendall 1998, eq. 1)
    }
    envstochvar <- totalvar - demvar
    results[rate,] <- c(totalvar,demvar,envstochvar)
  }
  results
} 
