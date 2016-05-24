mtdt <- function(x,n.sim=0)
{
#
# lower triangular matrix 5/6/2004
#
  tril <- function(t)
  {
    a <- t
    a[upper.tri(t)] <- 0
    a
  }

#
# upper triangular matrix 5/6/2004
#
  triu <- function(t)
  {
    a <- t
    a[lower.tri(t)] <- 0
    a
  }

#
# Spielman & Ewens' statistic (diagnoal elements kept)
# SAS 20/07/1999
#
  T <- x
  T <- T-diag(diag(T))
  m <- dim(T)[1]
#
# Find non-zero off-diagonal elements of table.
#
  b <- matrix(rep(0,m*m),nrow=m)
  i <- j <- matrix(rep(0,m*(m-1)/2),ncol=1)
  for (ii in 1:m)
  {
     for (jj in 1:(ii-1)) b[ii,jj]=T[ii,jj]+T[jj,ii]
  }
  l <- 0
  for (ii in 1:(m-1))
  {
     for (jj in (ii + 1):m)
        if (b[jj,ii]>0)
        {
           l <- l + 1
           i[l] <- jj
           j[l] <- ii
        }
  }
#
# Count number of different heterozygotes.
#
  NH <- l
#
# Count number of each heterozygote.
#
  t0 <- T
  tc <- apply(t0,2,sum)
  tr <- apply(t0,1,sum)
  c0 <- (m-1)/m
  se0 <- c0*sum((tc-tr)^2/(tc+tr))
  if (n.sim>0)
  {
    C <- rep(0,NH)
    for (k in 1:NH) C[k] <- T[i[k],j[k]] + T[j[k],i[k]]
    X <- rep(0,n.sim)
    for (k in 1:n.sim)
    {
       for (L in 1:NH)
       {
          T[i[L],j[L]] <- rbinom(1,C[L],0.5) # sum(runif(C(L))<0.5)
          T[j[L],i[L]] <- C[L] - T[i[L],j[L]]
       }
       tc <- apply(T,2,sum)
       tr <- apply(T,1,sum)
       X[k] <- c0*sum((tc-tr)^2/(tc+tr))
    }
    MCp <- 0
    for (k in 1:n.sim) if (X[k]>=se0) MCp <- MCp + 1
    pSE <- MCp/n.sim
    sSE <- sqrt(pSE*(1-pSE)/n.sim)
    cat('Spielman-Ewens Chi-square and empirical p (se): ', se0, pSE, sSE, "\n")
  } else cat('Spielman-Ewens Chi-square: ', se0, "\n")
#
# Simulate tables and compute TDT chi-square statistics.
#
# Should diag(T) is kept, the statistic will be similar to Spielman-Ewens'
# se.check

  T <- x

# This is according to Mike Miller's Matlab program
# Produce inverse (IV) of variance-covariance matrix (V).  This is
# constant across repeated samples in the Monte Carlo simulation.

  V <- diag(apply(T,1,sum)+apply(T,2,sum))-tril(T)-t(triu(T))-t(tril(T))-triu(T)
  IV <- solve(V[1:(m-1),1:(m-1)])

  T0 <- T
  d0=apply(T0[1:(m-1),1:(m-1)],1,sum)-apply(T0[1:(m-1),1:(m-1)],2,sum)
  x0 <- t(d0)%*%IV%*%d0
  se.check <- x0

  T <- T - diag(diag(T))
  V <- diag(apply(T,1,sum)+apply(T,2,sum))-tril(T)-t(triu(T))-t(tril(T))-triu(T)
  IV <- solve(V[1:(m-1),1:(m-1)])

  d0=apply(T0[1:(m-1),1:(m-1)],1,sum)-apply(T0[1:(m-1),1:(m-1)],2,sum)
  st0 <- t(d0)%*%IV%*%d0
  if (n.sim>0)
  {
    C <- rep(0,NH)
    for (k in 1:NH) C[k] <- T[i[k],j[k]] + T[j[k],i[k]]
    X <- rep(0,n.sim)
    for (k in 1:n.sim)
    {
       for (L in 1:NH)
       {
          T[i[L],j[L]] <- rbinom(1,C[L],0.5) # sum(runif(C(L))<0.5)
          T[j[L],i[L]] <- C[L]-T[i[L],j[L]]
       }
       d <- apply(T[1:(m-1),1:(m-1)],1,sum)-apply(T[1:(m-1),1:(m-1)],2,sum)
       X[k] <- t(d)%*%IV%*%d
    }
    MCp <- 0
    for (k in 1:n.sim) if (X[k]>st0) MCp <- MCp + 1
    pST <- MCp/n.sim
    sST <- sqrt(pST*(1-pST)/n.sim)
    cat('Stuart Chi-square and p (se): ', st0, pST, sST,"\n")
  } else {
    cat('Stuart Chi-square ',st0, "\n")
    cat('Value of Chi-square if diagonal elements are kept: ',se.check,"\n")
  }
  if (n.sim>0) list(SE=se0,pSE=pSE,sSE=sSE,ST=st0,pST=pST,sST=sST)
  else list(SE=se0,ST=st0)
}
