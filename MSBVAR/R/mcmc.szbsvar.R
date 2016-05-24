"mcmc.szbsvar" <- function(varobj, A0.posterior)
  {
    m<-dim(varobj$ar.coefs)[1]  # Capture the number of variablesbcoefs <- varobj$Bhat
    p<-dim(varobj$ar.coefs)[3]    # Capture the number of lags

    ncoef <- dim(varobj$B.posterior)[1]
    n0 <- varobj$n0
    n0cum <- c(0,cumsum(n0))
    N2 <- A0.posterior$N2

  # Get the covar for the coefficients
    XXinv <- chol(solve(varobj$Hpinv.posterior[[1]]))

  # storage for the sampled coefficients.
    B.sample <- matrix(0, nrow=N2, ncol=ncoef*m)

  # Loop for the sampling
    for(i in 1:N2)
    {
      # Set up the A0 for this iteration
        A0 <- A0.get(A0.posterior$A0.posterior, i)
        A0inv <- solve(A0)

        bj <- a2b(A0, varobj$Ui)

        F.draw <- matrix(0, ncoef, m)

        for(j in 1:m)
        { btmp <- bj[(n0cum[j]+1):(n0cum[(j+1)])]
          F.draw[,j] <- varobj$P.posterior[[j]]%*%(btmp)
        }

        F.draw <- F.draw + XXinv%*%matrix(rnorm(m*ncoef), ncoef, m)
        B.draw <- F.draw%*%(A0inv)

        B.sample[i,] <- matrix(B.draw, 1, ncoef*m)
        if (i%%1000==0)
        { cat("Monte Carlo Iteration = ", i, "\n"); }
    }

    # output / return
    output <- list(B.sample=B.sample)
    class(output) <- c("mcmc.posterior.BSVAR")
    return(output)

  }

