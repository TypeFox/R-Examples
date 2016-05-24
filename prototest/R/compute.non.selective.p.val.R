compute.non.selective.p.val <-
function (ts, x, type, df1, df2, mu, sigma, M){
  if (type == "ELR"){
    return (pchisq (ts, df=1, lower.tail=FALSE))
  }
  if (type == "ALR"){
    if (is.null(mu)){ # do not have an exact distribution for unknown mu case
      return (pchisq (ts, df=1, lower.tail=FALSE))
    }else{ # can compute exact distribution
      p.up = pchisq (M + sigma*sqrt(2*M*ts), df=M, lower.tail=TRUE, log.p=TRUE)
      p.down = pchisq (M - sigma*sqrt(2*M*ts), df=M, lower.tail=TRUE, log.p=TRUE)
      
      return (1 - exp(p.up + log(1 - exp(p.down-p.up))))
    }
  }
  if (type == "F"){
    return (pf (ts, df1=df1, df2=df2, lower.tail=FALSE))
  }
  if (type == "MS"){
    return (2*pnorm(ts, 0, 1, lower.tail=FALSE))
  }
}
