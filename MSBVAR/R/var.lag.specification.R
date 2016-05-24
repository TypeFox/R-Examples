"var.lag.specification" <-
function(y, lagmax=20)
{ results <- matrix(0,nrow=lagmax,ncol=4)
  ldets <- matrix(0,nrow=lagmax,ncol=4)

  for (p in 1:lagmax)
    { x <- y[(lagmax-p+1):nrow(y),]
      tmp <- ar.ols(x, order.max=p, aic=F, intercept=T)
      no.param <- length(tmp$ar) + length(tmp$x.intercept)
      obs <- nrow(y) - lagmax
      ldets[lagmax-p+1,2] <- log(det(tmp$var.pred))
      ldets[lagmax-p+1,1] <- p
      aic <- ldets[lagmax-p+1,2] + no.param*2/obs
      bic <- ldets[lagmax-p+1,2] + no.param*log(obs)/obs
      hq <- ldets[lagmax-p+1,2] + 2*log(log(obs))*(no.param/obs)
      results[p,1] <- p
      results[p,2] <- aic
      results[p,3] <- bic
      results[p,4] <- hq
    }

  colnames(results) <- c("Lags","AIC","BIC", "HQ")

  # now do LR tests
  for (i in 1:(lagmax-1))
    { mcorr <- ncol(y)*ldets[i,1] + 1
      chi <- (obs-mcorr)*(ldets[i+1,2] - ldets[i,2])
      df <- ncol(y)^2
      chi.p <- 1-pchisq(chi,df)
      ldets[i,3] <- chi
      ldets[i,4] <- chi.p
    }

  colnames(ldets) <- c("Lags","Log-Det","Chi^2","p-value")
  output <- list(ldets=ldets,results=results)
  class(output) <- c("var.lag.specification")
  return(output)
}

