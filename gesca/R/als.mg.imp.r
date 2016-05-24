als.mg.imp <- function (Z, rawz0, W0, A0, W, A, V, I, PHT, nvar, nlv, ng, missingvalue, itmax, ceps)
{
	#---------------------------------------------------
	# ALS algorithm for GSCA (multiple group analysis)
	# Heungsun Hwang, Sunmee Kim
	# No calcuation of Kronecker products  
	# Last revised Jan 29, 2016 
	# Elements of A and W can be fixed to constants (0 > & < 1)
	# LS imputation of missing observations (moption  == 3)
	# Input arguments: Z = data matrix for all groups (normalized)
	#                  nvar = num of indicators per group
	#                  nlv = num of latents per group
	#                  ng = num of groups
	#                  itmax = maximum number of iterations (default = 100)
	#                  ceps = convergence tolerance (default = 0.00001)
	#---------------------------------------------------

	sizez <- ncol(Z)
	sizea <- ncol(A)
	aindex0 <- which(A0 >= 1)
	Psi <- Z%*%V%*%I
	Gamma <- Z%*%W
	it <- 0			# iteration counter
	imp <- 100000	# initial improvement
	f0 <- 10^10 	# initial function value
	
	while (it <= itmax && imp > ceps) {
		it <- it + 1
		# step 1: update A
		for (t in 1:sizea) {
			if ( !all(A0[,t] == 0) ) {
				H1 <- diag(1,sizea)
				H1[t,t] <- 0
				aindex <- which(A0[,t] >= 1) 	# free parameters
				if (length(aindex) != 0) {
					a <- A[,t,drop=FALSE]
					a[aindex] <- 0					# if fixed values, not all zeros
					e <- matrix(0,1,sizea)
					e[t] <- 1
					Y <- Psi - Gamma%*%(A%*%H1 - a%*%e)
					X <- Gamma[,aindex,drop=FALSE]
					A[aindex,t] <- solve(t(X)%*%X, t(X)%*%Y%*%t(e))
				}
			}
		}
		vecA <- A[aindex0]
		A[aindex0] <- PHT%*%vecA
		
		# step 2: update W
		kk = 0
		for (g in 1:ng) {
			k <- kk + 1
			kk <- kk + nlv
			p <- g*nvar + (g-1)*nlv
			s <- 0
			for (j in k:kk) {
				s <- s + 1
				windex <- which(W0[,j] == 99)
				w <- W[,j,drop=FALSE]
				w[windex] <- 0 # if fixed values, not all zeros
				beta <- I[p+s,,drop=FALSE] - A[j,,drop=FALSE]
				H2 <- diag(1,ng*nlv)
				H2[j,j] <- 0
				H3 <- diag(1,sizea*ng)
				H3[p+s, p+s] <- 0
				Delta <- W%*%H2%*%A - V%*%H3%*%I - w%*%beta
				Zp <- Z[,windex]
				temptheta <- solve(t(Zp)%*%Zp,t(Zp))%*%(Z%*%Delta)%*%t(beta)
				theta <- t(solve(t(beta%*%t(beta)),t(temptheta)))
				zw <- Zp%*%theta
				theta <- theta%*%(1/sqrt(t(zw)%*%zw))
				W[windex,j] <- theta
				V[windex,p+s] <- theta
			}
		}
		
		Q <- V%*%I - W%*%A
		for (j in 1:sizez) {
			zindex <- as.matrix(which(rawz0[,j] == missingvalue))
			if ( nrow(zindex) != 0 ) {
				z <- Z[,j,drop=FALSE]
				q <- Q[j,,drop=FALSE]
				H4 <- diag(1,sizez)
				H4[j,j] <- 0
				Y <- -Z%*%H4%*%Q
				mz <- (Y%*%t(q))%*%(1/q%*%t(q))
				z[zindex] <- mz[zindex]
				Z[,j] <- z%*%(1/sqrt(t(z)%*%z))
			}
		}
		
		Gamma <- Z%*%W
		Psi <- Z%*%V%*%I
		dif <- Psi-Gamma%*%A
		f <- sum(diag(t(dif)%*%dif))
		imp <- f0-f
		info <- c(it,f,f0,imp)
		f0 <- f
	}

	output.als.mg.imp <- list(W = W, A = A, Z = Z, Psi = Psi, Gamma = Gamma, f = f, it = it, imp = imp)
	output.als.mg.imp
}