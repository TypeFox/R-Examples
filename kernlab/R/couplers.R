## wrapper function for couplers
## author : alexandros karatzoglou

couple <- function(probin, coupler = "minpair")
{
  if(is.vector(probin))
    probin <- matrix(probin,1)
    m <- dim(probin)[1]
  
  coupler <- match.arg(coupler, c("minpair", "pkpd", "vote", "ht"))

#  if(coupler == "ht")
#    multiprob <- sapply(1:m, function(x) do.call(coupler, list(probin[x ,], clscnt))) 
#  else
    multiprob <- sapply(1:m, function(x) do.call(coupler, list(probin[x ,])))    

  return(t(multiprob))
}


ht <- function(probin, clscnt, iter=1000)
  {
    nclass <- length(clscnt)
    probim <- matrix(0, nclass, nclass)
    for(i in 1:nclass)
      for(j in 1:nclass)
        if(j>i)
          {
            probim[i,j] <- probin[i]
            probim[j,i] <- 1 - probin[i]
          }
  
    p <- rep(1/nclass,nclass)
    u <- matrix((1/nclass)/((1/nclass)+(1/nclass)) ,nclass,nclass)
    iter <- 0
  
    while(TRUE)
      {
        iter <- iter + 1
        stoperror <- 0

        for(i in 1:nclass){
          num <- den <- 0
          for(j in 1:nclass)
            {
              if (j!=i)
                {
                  num <- num + (clscnt[i] + clscnt[j]) * probim[i,j] 
                  den <- den + (clscnt[i] + clscnt[j]) * u[i,j]  
                }
            }
          alpha <- num/(den + 1e-308)
          p[i] <- p[i]*alpha
          stoperror <- stoperror + (alpha -1)^2
          if(0)
            {
              sum <- 0
              sum <- sum(p) + sum
              p <- p/sum
              for(ui in 1:nclass)
                for(uj in 1:nclass)
                  u[ui, uj] <- p[ui]/(p[ui] + p[uj])
            }
          else
            {
              for(j in 1:nclass)
                if (i!=j)
                  {
                    u[i,j] <- p[i]/(p[i] + p[j])
                    u[j,i] <- 1 - u[i,j]
                  }
            }
        }
        if(stoperror < 1e-3)
          break
        if(iter > 400)
          {
            cat("Too many iterations: aborting", probin, iter, stoperror, p)
            break
          }
      }
    ## normalize prob.
    p <- p/sum(p)
    return(p)
  }


minpair <- function(probin)
  {  ## Count number of classes and construct prob. matrix
    nclass <- (1+sqrt(1 + 8*length(probin)))/2
    if(nclass%%1 != 0) stop("Vector has wrong length only one against one problems supported")
    probim <- matrix(0, nclass, nclass)
    probim[upper.tri(probim)] <- probin
    probim[lower.tri(probim)] <- 1 - probin
    
    sum <- colSums(probim^2)
    Q <- diag(sum)
    Q[upper.tri(Q)] <- - probin*(1 - probin)
    Q[lower.tri(Q)] <- - probin*(1 - probin)
    SQ <- matrix(0,nclass +1, nclass +1)
    SQ[1:(nclass+1) <= nclass, 1:(nclass+1) <= nclass] <- Q
    SQ[1:(nclass+1) > nclass, 1:(nclass+1) <= nclass] <- rep(1,nclass)
    SQ[1:(nclass+1) <= nclass, 1:(nclass+1) > nclass] <- rep(1,nclass)
    
    rhs <- rep(0,nclass+1)
    rhs[nclass + 1] <- 1

    p <- solve(SQ,rhs)

    p <- p[-(nclass+1)]/sum(p[-(nclass+1)])
    return(p)
  }


pkpd <- function(probin)
  {  ## Count number of classes and constuct prob. matrix
    nclass <- k <- (1+sqrt(1 + 8*length(probin)))/2
    if(nclass%%1 != 0) stop("Vector has wrong length only one against one problems supported")
    probim <- matrix(0, nclass, nclass)
    probim[upper.tri(probim)] <- probin
    probim[lower.tri(probim)] <- 1 - probin
    
    probim[probim==0] <- 1e-300
    R <- 1/probim
    diag(R)  <-  0
    p <- 1/(rowSums(R) - (k-2))

    p <- p/sum(p)
    return(p)
  }
    
    
vote<- function(probin)
{
  nclass <- (1+sqrt(1 + 8*length(probin)))/2
  if(nclass%%1 != 0) stop("Vector has wrong length only one against one problems supported")
   
  votev <- rep(0,nclass)
  p <- 0
  for(i in 1:(nclass-1))
    {
      jj <- i+1
      for(j in jj:nclass)
        {
          p <- p+1
          votev[i][probin[i] >= 0.5] <- votev[i][probin[i] >= 0.5] + 1
          votev[j][probin[j] < 0.5] <- votev[j][probin[j] < 0.5] + 1
        }
    }
  
  p <- votev/sum(votev)
  return(p)
}


