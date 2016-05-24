icfit <-
function(L, R, initp = NA, minerror = 1e-06, maxcount = 1000)
{
##############################################################
# for discussion of algorithm and theory 
#see Gentleman and Geyer (1994) Biometrika 81:618-623
#############################################################
	n <- length(L)
	if(n != length(R))
		stop("length of the two interval vectors must be the same")
	theta <- sort(unique(c(L, R, 0, Inf)))
	k <- length(theta)	
	## allow L[i]==R[i], but must adjust it so that L[i] equals the 
## next smaller value of theta
	if(any(L == R)) {
		exacts <- sort(unique(R[R == L]))
		if(exacts[1] == 0)
			stop("L[i]==R[i]=0 for some i")
		for(j in 1:length(exacts))
			L[R == L & L == exacts[j]] <- theta[(1:k)[theta == 
				exacts[j]] - 1]
	}
	A <- matrix(0, n, k)
	for(i in 1:n) {
		A[i, L[i] < theta & theta <= R[i]] <- 1
	}
##################################################################
### perform primary reductions on A  
### see Aragon and Eberly (1992) J of Computational and Graphical
###     Statistics 1:129-140 for discussion of primary reduction
##################################################################
	colsums <- apply(A, 2, sum)
	pairmult <- rep(0, k - 1)
	mark.to.keep <- rep(TRUE, k) # change T to TRUE
	for(i in 1:(k - 1)) {
		pairmult[i] <- sum(A[, i] * A[, i + 1])
		if(pairmult[i] == colsums[i]) {
			if(colsums[i] < colsums[i + 1]) {
				mark.to.keep[i] <- FALSE # change F to False
			}
		}
		if(pairmult[i] == colsums[i + 1]) {
			if(colsums[i] >= colsums[i + 1]) {
				mark.to.keep[i + 1] <- FALSE # change F to False
			}
		}
	}
	A <- A[, mark.to.keep]	### come up with the initial estimates
## Replace inefficient code for more efficient code
# p <- matrix(1, n, k)
	if(any(is.na(initp))) {
## Replace inefficient code for more efficient code
#  for(i in 1:n) {
#   p[i,  ] <- A[i,  ]/sum(A[i,  ])
#  }
#  pbar <- apply(p, 2, mean)
		pbar <- apply(A/apply(A, 1, sum), 2, mean)
	}
	else {
		if(length(initp) != k)
			stop("initp not of proper length")
		if(sum(initp[mark.to.keep]) != 1) {
			warning("after primary reduction, sum of initp !=1")
			initp <- initp/sum(initp[mark.to.keep])
		}
		pbar <- initp[mark.to.keep]
	}
	error <- 1
	count <- 1
	u <- -1
	while(error > minerror & count < maxcount) {
### algorithm to improve initial estimates
		pbar[pbar < minerror] <- 0	
	## Replace inefficient code for more efficient code
#  for(i in 1:n) {
#   p[i,  ] <- A[i,  ] * (pbar/sum(pbar * A[i,  ]))
# }
#  pbar <- apply(p, 2, mean)
		temp <- A/as.vector(A %*% pbar)
		pbar <- apply(t(temp) * pbar, 1, mean)
		count <- count + 1
		d <- apply(temp, 2, sum)
		u <-  - d + n
		u[pbar > 0] <- 0
		error <- max(d + u - n)
	}
#### test the Kuhn-Tucker conditions
	if(any(u < 0))
		warning("problem with convergence, decrease minerror")
	if(count == maxcount)
		warning("problem with convergence, increase maxcount")
	temppbar <- rep(0, k)
	temppbar[mark.to.keep] <- pbar
	surv <- rep(0, k)
	for(i in 1:(k - 1)) {
		surv[i] <- sum(temppbar[(i + 1):k])
	}
	names(temppbar) <- as.character(theta)
	names(surv) <- as.character(theta)	
	## Since theta[k]==Inf, take off those values
	theta <- theta[ - k]
	surv <- surv[ - k]	
	## but leave it on the density so it sums to 1, and it may
## be reentered in as initp if needed
	out <- list(u = u, error = error, count = count, p = temppbar, time = 
		theta, surv = surv)
	out
}

