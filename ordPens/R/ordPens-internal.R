coding <- function(x, constant=TRUE, splitcod=TRUE)
  {
  # Converts ordinal/categorical data into split-coding, or dummy-coding.

  # Input:
  # x: data matrix with ordinal/categorical data
  # constant: should a constant added in the regression

  # Output:
  # Split/dummy-coded data-matrix
  # Reference-Category = 1
  n <- nrow(x)
  kx <- apply(x,2,max) - 1
  xds <- matrix(0,n,sum(kx))

  # Loop to run through the columns of x
  for (j in 1:ncol(x))
		{
      j1col <- ifelse(j>1,sum(kx[1:(j-1)]),0)
  		# Loop to run through the rows of x
  		for (i in 1:n)
  			{
          if (x[i,j] > 1)
            {
              if (splitcod)
                xds[i,j1col + 1:(x[i,j]-1)] <- 1
              else
                xds[i,j1col + (x[i,j]-1)] <- 1
            }
  		  }
    }
  ## Output
  if (constant)
    return(cbind(1,xds))
  else
    return(xds)
 	}



genRidge <- function(x, y, offset, omega, lambda, model, delta=1e-6, maxit=25)
  {
    coefs <- matrix(0,ncol(x),length(lambda))
    fits <- matrix(NA,nrow(x),length(lambda))
    l <- 1
    for (lam in lambda)
      {
        if (model == "linear")
          {
            yw <- y - offset
            chollam <- chol(t(x)%*%x + lam*omega)
            coefs[,l] <- backsolve(chollam,
            backsolve(chollam, t(x)%*%yw, transpose=TRUE))
            fits[,l] <- x%*%coefs[,l] + offset
          }
        else
          {
            ## penalized logistic/poisson regression
            conv <- FALSE
    
            # start values
            if (l > 1)
              {
                b.start <- coefs[,l-1]
              }
            else
              {
                if (model == "logit")
                  {
                    gy <- rep(log(mean(y)/(1-mean(y))),length(y)) - offset
                  }
                else
                  {
                    gy <- rep(log(mean(y)),length(y)) - offset
                  }
                chollam <- chol(t(x)%*%x + lam*omega)
                b.start <- backsolve(chollam,
                backsolve(chollam, t(x)%*%gy, transpose=TRUE))
              }
    
            # fisher scoring
            b.old <- b.start
            i <- 1
            while(!conv)
              {
                eta <- x%*%b.old + offset
                if (model == "logit")
                  {
                    mu <- plogis(eta)
                    sigma <- as.numeric(mu*(1-mu))
                  }
                else
                  {
                    mu <- exp(eta)
                    sigma <- as.numeric(mu)
                  }
                score <- t(x)%*%(y-mu)
                fisher <- t(x)%*%(x*sigma)
                choli <- chol(fisher + lam*omega)
                b.new <- b.old + backsolve(choli,
                backsolve(choli, (score - lam*omega%*%b.old), transpose=TRUE))
  
                if(sqrt(sum((b.new-b.old)^2)/sum(b.old^2)) < delta | i>=maxit)
                  {
                    # check the stop criterion
                    conv <- TRUE
                  }
                b.old <- b.new
                i <- i+1
              }  # end while
            coefs[,l] <- b.old   
            fits[,l] <- mu       
          }
        l <- l+1
      }
    rownames(fits) <- NULL
    colnames(fits) <- lambda
    rownames(coefs) <- NULL
    colnames(coefs) <- lambda
    
    return(list(fitted = fits, coefficients = coefs))   
  }



cd <- function(x)
  {
    n <- length(x)

    X <- matrix(1,n,1)

    Z <- matrix(rep(x,max(x)-1),n,max(x)-1)
    Z <- Z - matrix(rep(1:(max(x)-1),n),dim(Z)[1],dim(Z)[2],byrow=T)
    Z[Z<0] <- 0
    Z[Z>1] <- 1

    res <- list(X = X, Z = Z)
    return(res)
  }



ordAOV1 <- function(x, y, type = "RLRT", nsim = 10000, null.sample = NULL ,...){

  x <- as.numeric(x)

  ## check x
  if(min(x)!=1 | length(unique(x)) != max(x))
    stop("x has to contain levels 1,...,max")
    
  if(length(x) != length(y))
    stop("x and y have to be of the same length")

  k <- length(unique(x))
  cdx <- cd(x)
  X <- cdx$X
  Z <- cdx$Z

  # RLRT
  if (type == "RLRT")
  {
  # model under the alternative
  m <- gam(y ~ Z, paraPen=list(Z=list(diag(k-1))), method="REML")

  # null model
  m0 <- gam(y ~ 1, method="REML")

  # log-likelihood
  rlogLik.m <- -summary(m)$sp.criterion
  rlogLik.m0 <- -summary(m0)$sp.criterion
  rlrt.obs <- max(0, 2*(rlogLik.m - rlogLik.m0))

  # null distribution
  if (rlrt.obs != 0) {
      if (length(null.sample)==0)
        RLRTsample <- RLRTSim(X, Z, qr(cdx$X), chol(diag(k-1)), nsim = nsim, ...)
      else
        RLRTsample <- null.sample
        
      p <- mean(rlrt.obs < RLRTsample)
    }
  else
    {
      if (length(null.sample)==0)
        RLRTsample <- NULL
      else
        RLRTsample <- null.sample

      p <- 1
    }

  # return
  RVAL <- list(statistic = c(RLRT = rlrt.obs), p.value = p,
        method = paste("simulated finite sample distribution of RLRT.\n (p-value based on",
        length(RLRTsample), "simulated values)"), sample = RLRTsample)
  }

  # LRT
  else
  {
  # model under the alternative
  m <- gam(y ~ Z, paraPen=list(Z=list(diag(k-1))), method="ML")

  # null model
  m0 <- gam(y ~ 1, method="ML")

  # log-likelihood
  logLik.m <- -summary(m)$sp.criterion
  logLik.m0 <- -summary(m0)$sp.criterion
  lrt.obs <- max(0, 2*(logLik.m - logLik.m0))

  # null distribution
  if (lrt.obs != 0) {
      if (length(null.sample)==0)
        LRTsample <- LRTSim(X, Z, q=0, chol(diag(k-1)), nsim = nsim, ...)
      else
        LRTsample <- null.sample
        
      p <- mean(lrt.obs < LRTsample)
    }
  else
    {
      if (length(null.sample)==0)
        LRTsample <- NULL
      else
        LRTsample <- null.sample

      p <- 1
    }

  # return
  RVAL <- list(statistic = c(LRT = lrt.obs), p.value = p,
        method = paste("simulated finite sample distribution of LRT.\n (p-value based on",
        length(LRTsample), "simulated values)"), sample = LRTsample)
  }

  class(RVAL) <- "htest"
  return(RVAL)
}



ordAOV2 <- function(x, y, type = "RLRT", nsim = 10000, null.sample = NULL, ...){

  n <- length(y)
  p <- ncol(x)
  k <- apply(x, 2, max)

  ## check x
  if(min(x[,1])!=1 | length(unique(x[,1])) != max(x[,1]))
    stop("x has to contain levels 1,...,max")

  if(nrow(x) != length(y))
    stop("nrow(x) and length(y) do not match")

  cdx <- cd(x[,1])
  X <- cdx$X
  ZZ <- vector("list", p)
  ZZ[[1]] <- cdx$Z
  Z <- ZZ[[1]]
  DD <- vector("list", p)
  DD[[1]] <- diag(c(rep(1,k[1]-1),rep(0,sum(k[2:p]-1))))
  for (j in 2:p)
    {

      ## check x
      if(min(x[,j])!=1 | length(unique(x[,j])) != max(x[,j]))
        stop("x has to contain levels 1,...,max")

      ZZ[[j]] <- cd(x[,j])$Z
      Z <- cbind(Z,ZZ[[j]])
      if (j < p)
        DD[[j]] <- diag(c(rep(0,sum(k[1:(j-1)]-1)),rep(1,k[j]-1),rep(0,sum(k[(j+1):p]-1))))
      else
        DD[[j]] <- diag(c(rep(0,sum(k[1:(j-1)]-1)),rep(1,k[j]-1)))
    }
  RRVAL <- vector("list", p)

  # RLRT
  if (type == "RLRT")
  {
  # model under the alternative
  mA <- gam(y ~ Z, paraPen=list(Z=DD), method="REML")

  for (j in 1:p)
  {
  # null model
  Z0 <- matrix(unlist(ZZ[-j]),n,sum(k[-j]-1))
  D0 <- DD[-j]
  if(!is.list(D0))
    D0 <- list(D0)

  if (j == 1)
    out <- 1:(k[1]-1)
  else
    out <- sum(k[1:(j-1)]-1)+(1:(k[j]-1))

  for(j0 in 1:length(D0))
    {
      D0[[j0]] <- D0[[j0]][-out,-out]
    }
  m0 <- gam(y ~ Z0, paraPen=list(Z0=D0), method="REML")

  # log-likelihood
  rlogLik.mA <- -summary(mA)$sp.criterion
  rlogLik.m0 <- -summary(m0)$sp.criterion
  rlrt.obs <- max(0, 2*(rlogLik.mA - rlogLik.m0))

  # model that only contains the variance
  # ...

  # null distribution
  Z1 <- ZZ[[j]]
  if (rlrt.obs != 0) {
      if (length(null.sample) == 0)
        RLRTsample <- RLRTSim(X, Z1, qr(X), chol(diag(k[j]-1)), nsim = nsim, ...)
      else
        {
          if (length(null.sample) == p)
            RLRTsample <- null.sample[[j]]
          else
            stop("wrong number of null.sample elements")
        }

      p <- mean(rlrt.obs < RLRTsample)
    }
  else
    {
      if (length(null.sample) == 0)
        RLRTsample <- NULL
      else
        {
          if (length(null.sample) == p)
            RLRTsample <- null.sample[[j]]
          else
            stop("wrong number of null.sample elements")
        }

      p <- 1
    }

  # return
  RVAL <- list(statistic = c(RLRT = rlrt.obs), p.value = p,
        method = paste("simulated finite sample distribution of RLRT.\n (p-value based on",
        length(RLRTsample), "simulated values)"), sample = RLRTsample)

  class(RVAL) <- "htest"
  RRVAL[[j]] <- RVAL
  }
  }

  # LRT
  else
  {
  # model under the alternative
  mA <- gam(y ~ Z, paraPen=list(Z=DD), method="ML")

  for (j in 1:p)
  {
  # null model
  Z0 <- matrix(unlist(ZZ[-j]),n,sum(k[-j]-1))
  D0 <- DD[-j]
  if(!is.list(D0))
    D0 <- list(D0)

  if (j == 1)
    out <- 1:(k[1]-1)
  else
    out <- sum(k[1:(j-1)]-1)+(1:(k[j]-1))

  for(j0 in 1:length(D0))
    {
      D0[[j0]] <- D0[[j0]][-out,-out]
    }
  m0 <- gam(y ~ Z0, paraPen=list(Z0=D0), method="ML")

  # log-likelihood
  logLik.mA <- -summary(mA)$sp.criterion
  logLik.m0 <- -summary(m0)$sp.criterion
  lrt.obs <- max(0, 2*(logLik.mA - logLik.m0))

  # model that only contains the variance
  # ...

  # null distribution
  Z1 <- ZZ[[j]]
  if (lrt.obs != 0) {
      if (length(null.sample) == 0)
      LRTsample <- LRTSim(X, Z1, q=0, chol(diag(k[j]-1)), nsim = nsim, ...)
      else
        {
          if (length(null.sample) == p)
            LRTsample <- null.sample[[j]]
          else
            stop("wrong number of null.sample elements")
        }

      p <- mean(lrt.obs < LRTsample)
    }
  else
    {
      if (length(null.sample) == 0)
        LRTsample <- NULL
      else
        {
          if (length(null.sample) == p)
            LRTsample <- null.sample[[j]]
          else
            stop("wrong number of null.sample elements")
        }

      p <- 1
    }

  # return
  RVAL <- list(statistic = c(LRT = lrt.obs), p.value = p,
        method = paste("simulated finite sample distribution of LRT.\n (p-value based on",
        length(LRTsample), "simulated values)"), sample = LRTsample)

  class(RVAL) <- "htest"
  RRVAL[[j]] <- RVAL
  }
  }

  names(RRVAL) <- colnames(x)
  return(RRVAL)
}
