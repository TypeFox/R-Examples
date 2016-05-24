ProDenICA <-
function(x, k=p, W0=NULL, whiten=FALSE, maxit = 20, thresh = 1e-7, restarts = 0,
	trace = FALSE, Gfunc=GPois, eps.rank=1e-7, ...)
{
  this.call=match.call()
  p <- ncol(x)
  n <- nrow(x)
  x <- scale(x, T, F)## x should have mean zero
  if(whiten){## First sphere the data
    sx <- svd(x)
    ##Get the effective rank
    condnum=sx$d;condnum=condnum/condnum[1]
    good=condnum >eps.rank
    rank=sum(good)
    if(k>rank){
      warning(paste("Rank of x is ",rank,"; k reduced from",k," to ",rank,sep=""))
      k=rank
    }
    x <- sqrt(n) * sx$u[,good]# no need to rotate via  %*% t(sx$v[,good])
    whitener=sqrt(n)*scale(sx$v[,good],FALSE,sx$d[good])
  }
  else whitener=NULL
### Get a random start if needed
        if(is.null(W0))	W0 <- matrix(rnorm(p * k), p, k) else   k=ncol(W0)
	W0 <- ICAorthW(W0)
###Initialization
	GS <- matrix(0., n, k)
	gS <- GS
	gpS <- GS
	s <- x %*% W0
	flist <- as.list(1.:k)
	for(j in 1.:k)
		flist[[j]] <- Gfunc(s[, j], ...)
	flist0 <- flist
	crit0 <- mean(sapply(flist0, "[[", "Gs"))
### can try some better starts; only evaluated at first iteration
	while(restarts) {
		W1 <- matrix(rnorm(p * k), p, k)
		W1 <- ICAorthW(W1)
		s <- x %*% W1
		for(j in 1.:k)
			flist[[j]] <- Gfunc(s[, j], ...)
		crit <- mean(sapply(flist, "[[", "Gs"))
		if(trace)
			cat("old crit", crit0, "new crit", crit, "\n")
		if(crit > crit0) {
			crit0 <- crit
			W0 <- W1
			flist0 <- flist
		}
		restarts <- restarts - 1.
	}
###Here is the loop
	nit <- 0
	nw <- 10
	repeat {
		nit <- nit + 1
		gS <- sapply(flist0, "[[", "gs")
		gpS <- sapply(flist0, "[[", "gps")
		t1 <- t(x) %*% gS/n
		t2 <- apply(gpS, 2, mean)
		W1 <- t1 - scale(W0, F, 1/t2)
		W1 <- ICAorthW(W1)
		nw <- amari(W0, W1)
		if(trace)
			cat("Iter", nit, "G", crit0, "crit", nw, "\n")
		W0 <- W1
		if((nit > maxit) | (nw < thresh))
			break
		s <- x %*% W0
		for(j in 1:k)
			flist0[[j]] <- Gfunc(s[, j], ...)
		crit0 <- mean(sapply(flist0, "[[", "Gs"))
	}
  rl=list(W = W0, negentropy = crit0,s= x %*% W0,whitener=whitener,call=this.call)
  rl$density=lapply(flist0,"[[","density")
  class(rl)="ProDenICA"
  rl
}

