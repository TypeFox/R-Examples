### compute for all combinations of correlations corr and
### margin probabilities mp the value of P(A,B)

simul.commonprob <- function(margprob, corr=0,
                             method="integrate", n1=10^5,n2=10)
{
  lm <- length(margprob)
  lr <- length(corr)
  
  z <- array(0, dim=c(lm,lm,lr))

  method <- pmatch(method, c("integrate", "monte carlo"))
  if(is.na(method))
    stop("invalid method")

  for(k in 1:lr){
    sigma <- matrix(c(1,corr[k],corr[k],1), ncol=2)
    for(m in 1:lm){
      q1 <- qnorm(margprob[m])
      for(n in m:lm){
        cat(corr[k], margprob[m], margprob[n], ": ")
        q2 <- qnorm(margprob[n])

        if(corr[k]==-1){
          z[m,n,k] <- max(margprob[m] + margprob[n] -1, 0)
          cat("done\n")
        }
        else if(corr[k]==0){
          z[m,n,k] <- margprob[m] * margprob[n]
          cat("done\n")
        }
        else if(corr[k]==1){
          z[m,n,k] <- min(margprob[m], margprob[n])
          cat("done\n")
        }
        else if(margprob[m] * margprob[n] == 0){
          z[m,n,k] <- 0
          cat("done\n")
        }
        else if(margprob[m] == 1){
          z[m,n,k] <- margprob[n]
          cat("done\n")
        }
        else if(margprob[n] == 1){
          z[m,n,k] <- margprob[m]
          cat("done\n")
        }
        else if(method==1) {
            a <- pmvnorm(c(0, 0), c(Inf, Inf),
                         mean = c(q1, q2), sigma = sigma)
            z[m,n,k] <- if(attr(a, "msg") != "Normal Completion")
                NA
            else
                c(a)
        }
        else{
          x2 <- rep(0,n2)
          for(l in 1:n2){
            x1 <- rmvnorm(n1, mean = c(q1,q2), sigma=sigma)
            x2[l] <- sum( (x1[,1] > 0) & (x1[,2] > 0))/n1
          }
          z[m,n,k] <- mean(x2)
          cat("done\n")
        }
        z[n,m,k] <- z[m,n,k]
      }
    }
  }
  dimnames(z)<-list(margprob, margprob, corr)
  z
}
