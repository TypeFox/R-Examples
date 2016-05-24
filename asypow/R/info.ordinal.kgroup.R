info.ordinal.kgroup <- function(p, group.size=1) {
###-----------------------------------------------------------------------
###      Returns the information matrix for the raw ordinal model
###
### p: For a single group with n categories,
###        p[i] = prob(event occured in category i or less)
###    Dimension of p in this case is (n-1) since p[n] would
###    automatically be 1.
###
###    For k groups a matrix dimensioned k X (n-1).
###
### group.size: The relative number of observations in each group
###
###
### Returns: The information matrix from a single trial spread over
###          the groups.
###
###-----------------------------------------------------------------------

  if (any(apply(cbind(0,p,1),1, diff) < 0))
    stop("p must be increasing and between 0 and 1")
  if (is.vector(p)) p <- t(as.matrix(p))
  ngroups <- dim(p)[1]

  lngrpsz <- length(group.size)
  if (lngrpsz == 1) group.size <- rep(1,ngroups) else 
        if (ngroups != lngrpsz)
      stop ("\nnumber or rows of p and length of group.size must match")

  hess <- vector("list", ngroups)
  for(j in 1:ngroups) {
    pj <- p[j,]
    k <- length(pj)
    if (k == 0) next
    if (k == 1) hess[[j]] <- -1/(pj*(1-pj)) else {
        hess[[j]] <- matrix(0,k,k)

        p.at <- diff(c(0,pj,1))
        t1 <- 1/p.at[-(k+1)]
        t2 <- 1/p.at[-1]
        diag(hess[[j]]) <- -(t1 + t2 )
        for(i in 1:(k-1)) hess[[j]][i,i+1] <- hess[[j]][i+1,i] <-
          1/p.at[i+1]
      }
    hess[[j]] <- group.size[j] * hess[[j]]
  }

  info <- - k.blocks.info(hess)

  return(info/sum(group.size))
}
