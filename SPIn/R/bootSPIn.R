bootSPIn <-
function(x, n.boot=50, conf = 0.95, bw = 0, lb = -Inf, ub = Inf, l=NA, u=NA){
	n.sims <- length(x)
	x <- sort(x)
	gaps <- x[2:n.sims] - x[1:(n.sims-1)]
	gaps <- c (gaps[1], gaps, gaps[n.sims-1])
	gap.bandwidth <- 2
	mean.gap <- rep (NA, n.sims)
	for (j in 1:n.sims){
		mean.gap[j] <- mean (gaps[max (1, j-(gap.bandwidth-1)) : min (n.sims+1, j+gap.bandwidth)])
	}
	theta.ordered.0 <- x
	if (lb != -Inf) 
        n.sims <- n.sims+1
    if (ub != Inf) 
        n.sims <- n.sims+1
	w.l <- matrix(0,nrow=n.sims,ncol=n.boot)
	w.u <- w.l
	for (i in 1:n.boot){
#		print(i)
		x <- sample(theta.ordered.0,n.sims,T)
		x <- x + mean.gap*runif (n.sims,-1,1)/20
		r <- SPIn(x,conf = conf, bw = bw, lb = lb, ub = ub, l=l, u=u)
		w.l[r$l.l:r$l.u,i] <- r$w.l
		w.u[r$u.l:r$u.u,i] <- r$w.u
	}
	x <- theta.ordered.0
	if (lb != -Inf) 
        x <- c(x, lb)
    if (ub != Inf) 
        x <- c(x, ub)
	x <- sort(x)
	w.l <- rowMeans(w.l)
	x1 <- w.l%*%x
	w.u <- rowMeans(w.u)
	x2 <- w.u%*%x
	l.ind <- which(w.l!=0)
	l.l <- l.ind[1]
	l.u <- l.ind[length(l.ind)]
	u.ind <- which(w.u!=0)
	u.l <- u.ind[1]
	u.u <- u.ind[length(u.ind)]

	hpd <- list(spin = c(x1, x2), conf = conf, x = x, w.l=w.l, w.u=w.u, l.l=l.l, l.u=l.u, u.l=u.l, u.u=u.u)
    class(hpd) <- "SPIn"
    return(hpd)
}
