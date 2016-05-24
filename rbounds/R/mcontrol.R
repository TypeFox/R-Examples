mcontrol <- function(y, grp.id, treat.id, group.size=3,
                     Gamma=4, GammaInc=1) {
  # Test Statistic
  grp.means <- tapply(y, grp.id, "mean")
  means <- kronecker(grp.means, matrix(1,group.size,1))
  values <- y - means
  a.ranks <- rank(values, ties.method="average")
  T <- sum(a.ranks[treat.id==1])
  # Housekeeping
  gamma.seq <- seq(1, Gamma, by = GammaInc)
  n.i <- group.size
  b <- length(unique(grp.id))
  mu <- matrix(0, b, 2)
  a <- seq(1, (n.i - 1), by = 1) 
  g.n <- length(gamma.seq)
  pvals <- matrix(NA, g.n, 2)
  
  for(k in 1:g.n) {
    # Max and Min Mean
    for(j in 1:b) {
      q <- sort(a.ranks[grp.id==j])
      denom.1 <- a + gamma.seq[k]*(n.i-a)
      denom.2 <- a + (1/gamma.seq[k])*(n.i-a)
      e.t <- function(q) {
        val.1 <- matrix(0,(n.i-1),1)
        val.2 <- matrix(0,(n.i-1),1)
        i <- 1
        for(i in 1:(n.i-1)) {
          num.1 <- sum(q[1:i]) + gamma.seq[k]*(sum(q[(i+1):n.i]))
          num.2 <- sum(q[1:i]) + (1/gamma.seq[k])*(sum(q[(i+1):n.i]))
          val.1[i,] <- num.1/denom.1[i]
          val.2[i,] <- num.2/denom.2[i]
        }
        ret <- c(max(val.1), min(val.2))
      }	
      
      mu[j,] <- e.t(q)
    }
    
    # Max and Min Variance
    V <- matrix(0, b, 2)
    for(j in 1:b) {
      q <- sort(a.ranks[grp.id==j])
      denom.1 <- a + gamma.seq[k]*(n.i-a)
      denom.2 <- a + (1/gamma.seq[k])*(n.i-a)
      e.v <-function(q) {
        val.1 <- matrix(0,(n.i-1),1)
        val.2 <- matrix(0,(n.i-1),1)
        i <- 1
        for(i in 1:(n.i-1)) {
          num.1 <- sum(q[1:i]^2) + gamma.seq[k]*(sum(q[(i+1):n.i]^2))
          num.2 <- sum(q[1:i]^2) + (1/gamma.seq[k])*(sum(q[(i+1):n.i]^2))
          val.1[i,] <- (num.1/denom.1[i]) - mu[j,1]^2
          val.2[i,] <- (num.2/denom.2[i]) - mu[j,2]^2
        }
        ret.v <- c(max(val.1), min(val.2))
      }	
      
      V[j,] <- e.v(q)
    }
    # Calculate P-values and Store
    z.up <- (T - sum(mu[,1]))/sqrt(sum(V[,1]))
    p.val.up <-1 - pnorm(z.up)
    z.lo <- (T - sum(mu[,2]))/sqrt(sum((V[,2])))
    p.val.low <- 1 - pnorm(z.lo)
    pvals[k,1] <- round(p.val.low, digits=4)
    pvals[k,2] <- round(p.val.up, digits=4)
  }
  
  pval <- pvals[1,1]
  bounds <- data.frame(gamma.seq, pvals)
  colnames(bounds) <- c("Gamma", "Lower bound", "Upper bound")
  
  msg <- "Rosenbaum Sensitivity Test for Wilcoxon Strat. Rank P-Value \n"
  note <- "Note: Gamma is Odds of Differential Assignment To
 Treatment Due to Unobserved Factors \n"
  
  Obj <- list(Gamma = Gamma, GammaInc = GammaInc, pval = pval,
              msg = msg, bounds = bounds, note = note)
  class(Obj) <- c("rbounds", class(Obj))
  
  Obj
}
