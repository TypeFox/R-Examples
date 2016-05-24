ClusterInit <-
function(yh, ng=5) {
 
  seq2 <- function(low,up,n) seq(low,up, length.out=n+2)[-c(1,n+2)]  # help function for filling gaps

  # Try 1 to ng clusters, decide on number of clusters, fill in gaps
  km <- list(); wss <- numeric(ng); BIC <- numeric(ng)
  for (i in 1:ng){
    km[[i]] <- kmeans(yh, centers=i, nstart=10*i)  # try different starts, increasing numbers for more clusters
    wss[i] <- km[[i]]$tot.withinss

    clus.mu <- km[[i]]$centers             #i*1 matrix of cluster centers (transformed scale)
    clus.sd <- sqrt(wss[i]/(length(yh)-i)) #one numeric value: overall within-cluster sd (transformed scale)
    clus.p <- km[[i]]$size/length(yh)      #vector of fraction of observations in each cluster

    if (i==1) loglik <- sum(log(sapply(yh, dnorm, clus.mu, clus.sd))) 
    else loglik <- sum(log(t(sapply(yh, dnorm, clus.mu, clus.sd)) %*% clus.p))
    BIC[i] <- -2*loglik + log(length(yh))*(7*i)   # use more heavy penalty than ordinary BIC to avoid spurious clusters
  }
  
  nrclust <- which.min(BIC)
  extraclust <- ng-nrclust

  clus.mu <- km[[nrclust]]$centers #ng*1 matrix
  o <- order(clus.mu)
  clus.mu <- clus.mu[o] #vector of length ng

  clus.sd <- rep(sqrt(wss[nrclust]/(length(yh)-nrclust)), ng)

  if (extraclust > 0) {

  # Fill the gaps, assuming that chosen mu's so far are consecutive; 
  # so, fill only on the left side and/or on the right side.

    maxy <- pi/2
    lo.mu <- min(clus.mu)
    hi.mu <- max(clus.mu)
    lo.empty <- lo.mu
    hi.empty <- pi/2-hi.mu
    frac.lo <- lo.empty/(lo.empty+hi.empty)
    frac.hi <- 1-frac.lo

    n.lo <- round(frac.lo * extraclust)
    n.hi <- round(frac.hi * extraclust)
    
    #if ((n.lo+n.hi) > extraclust) n.hi <- n.hi - 1
    #if ((n.lo+n.hi) < extraclust) n.lo <- n.lo + 1
    
    if ((n.lo+n.hi) > extraclust) {
      # happens only if n.lo and n.hi both rounded up, i.e. exact n in both cases x.5
      # decreasing the smallest one tends to shift distribution to one side,
      # decreasing the largest tends to center the distribution
      # we go for the shift:
      if (n.lo < n.hi) n.lo <- n.lo - 1 else n.hi <- n.hi - 1
    } else {
      if ((n.lo+n.hi) < extraclust) {
        # happens only if n.lo and n.hi both rounded down, i.e. exact n in both cases x.5
        # increasing the largest one tends to shift distribution to one side,
        # increasing the smallest tends to center the distribution
        # we go for the shift:
        if (n.hi > n.lo) n.hi <- n.hi + 1 else n.lo <- n.lo + 1
      }
    }
    clus.mu <- c(seq2(0, lo.mu, n.lo), clus.mu, seq2(hi.mu, maxy, n.hi))  
  }

  list(clus.mu=clus.mu, clus.sd=clus.sd) #both on transformed scale
}
