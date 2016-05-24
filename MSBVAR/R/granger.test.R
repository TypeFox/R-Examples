"granger.test" <- function(y, p)
  { m <- ncol(y)

    # Check that we have enough variables to sensibly do this.
    if(m<2)
      { stop(paste("Error: granger.test needs at least 2 variables"))
      }

    # Make objects to hold the results and names
    results <- matrix(0, m*(m-1),2)
    namelist <- vector(mode="character",m*(m-1))
    varnames <- dimnames(y)[[2]]
    
    k <- 0
    for(i in 1:m)
      { for (j in 1:m)
          {
            if(i==j) { next }
            Y <- embed(cbind(y[,i], y[,j]), p+1)
            X1 <- Y[, -(1:2)]
            X2 <- X1[, ((1:p)*2) - (1 %% 2)]
            restricted <- lm(Y[,1] ~ X2)
            unrestricted <- lm(Y[,1] ~ X1)
            
            ssqR <- sum(restricted$resid^2)
            ssqU <- sum(unrestricted$resid^2)

            ftest <- ((ssqR - ssqU)/p)/(ssqU/(nrow(Y) - 2*p - 1))

            # Save the results
            k <- k+1
            endog.name <- varnames[i]
            exog.name <- varnames[j]
            name <- paste(exog.name, "->", endog.name)
            namelist[k] <- name
            results[k,] <- c(ftest, 1 - pf(ftest, p, nrow(Y) - 2 * p - 1))

          }
      }
    rownames(results) <- namelist
    colnames(results) <- c("F-statistic", "p-value")
    return(results)
  }

