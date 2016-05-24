als.mg.ho2 <- function (Z, W01, W02, A0, W1, W2, A, V, I, PHT, nvar, nlv1, nlv, ng, itmax, ceps)
{
	#---------------------------------------------------
	# ALS algorithm for GSCA with second-order latents (multiple group analysis)
	# Heungsun Hwang, Sunmee Kim
	# No calcuation of Kronecker products  
	# Last revised Jan 29, 2016 
	# Elements of A and W can be fixed to constants (0 > & < 1)
	# Z = block-diagonal matrix of "normalized" data for all groups
	#---------------------------------------------------

	sizea <- ncol(A)
	aindex0 <- which(A0 >= 1)
	W <- W1%*%W2
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
		
		#step 2: update W1
		W2A <- W2%*%A
		kk <- 0
		ii <- 0
		ll <- 0
		for (g in 1:ng) {
			k <- kk + 1
			kk <- kk + nlv1
			i <- ii + 1
			ii <- ii + nvar
			l <- ll + 1
			ll <- ll + nlv
			s <- 0
			for (j in k:kk) {
				s <- s + 1
				t <- sizea*(g-1) + 1
				tt <- sizea*(g-1) + sizea
				windex1 <- which(W01[,j] == 99)
				w1 <- W1[,j,drop=FALSE]
				w1[windex1] <- 0   # if fixed values, not all zeros
				m <- cbind(matrix(0,1,nvar), W2[j,l:ll,drop=FALSE])
				H2 <- diag(1, nlv1)
				H2[s,s] <- 0
				Dj <- cbind(diag(1,nvar), W1[i:ii,k:kk,drop=FALSE]%*%H2%*%W2[k:kk,l:ll,drop=FALSE])
				D <- V
				D[i:ii,t:tt] <- Dj
				beta <- m - W2A[j,,drop=FALSE]
				H3 <- diag(1,ng*nlv1)
				H3[j,j] <- 0
				Delta <- W1%*%H3%*%W2A - D%*%I - w1%*%beta
				Y <- Z%*%Delta
				X <- Z[,windex1]
				temptheta <- solve(t(X)%*%X,t(X)) %*% Y %*% t(beta)
				theta <- t(solve(t(beta%*%t(beta)),t(temptheta)))
				zw <- X%*%theta
				theta <- theta%*%(1/sqrt(t(zw)%*%zw))
				W1[windex1,j] <- theta
				V[i:ii,t:tt] <- cbind(diag(1,nvar), W1[i:ii,k:kk,drop=FALSE]%*%W2[k:kk,l:ll,drop=FALSE])
			}
		}
		
		# step 3: update W2
		ZW1 <- Z%*%W1
		kk <- 0
		ii <- 0
		ll <- 0
		for (g in 1:ng) {
			k <- kk + 1
			kk <- kk + nlv1
			i <- ii + 1
			ii <- ii + nvar
			l <- ll + 1
			ll <- ll + nlv
			s <- 0
			for (j in l:ll) {
				s <- s + 1
				t <- sizea*(g-1) + 1
				tt <- sizea*(g-1) + sizea
				windex2 <- which(W02[,j] == 99)
				if ( nrow(as.matrix(windex2)) != 0 ) {
					w2 <- W2[,j,drop=FALSE]
					w2[windex2] <- 0	# if fixed values, not all zeros
					e <- matrix(0,1,nlv)
					e[s] <- 1
					m <- cbind(matrix(0,1,nvar),e)
					H4 <- diag(1,nlv)
					H4[s,s] <- 0
					Dj <- cbind(diag(1,nvar),W1[i:ii,k:kk,drop=FALSE]%*%W2[k:kk,l:ll,drop=FALSE]%*%H4)
					D <- V
					D[i:ii,t:tt] <- Dj
					beta <- m - A[j,,drop=FALSE]
					H5 <- diag(1,ng*nlv)
					H5[j,j] <- 0
					Delta <- W1%*%W2%*%H5%*%A - D%*%I      
					Y <- Z%*%Delta - ZW1%*%w2%*%beta
					X <- ZW1[,windex2]
					temptheta <- solve(t(X)%*%X,t(X)) %*% Y %*% t(beta)
					theta <- t(solve(t(beta%*%t(beta)),t(temptheta)))
					zw <- X%*%theta
					theta <- theta%*%(1/sqrt(t(zw)%*%zw))
					W2[windex2,j] <- theta
					V[i:ii,t:tt] <- cbind(diag(1,nvar),W1[i:ii,k:kk,drop=FALSE]%*%W2[k:kk,l:ll,drop=FALSE])
				}
			}
		}
		
		W <- W1%*%W2
		Gamma <- Z%*%W
		Psi <- Z%*%V%*%I
		dif <- Psi-Gamma%*%A
		f <- sum(diag(t(dif)%*%dif))
		imp <- f0-f
		info <- c(it,f,f0,imp)
		f0 <- f	
	}

	output.als.mg.ho2 <- list(W = W, W1 = W1, W2 = W2, A = A, Psi = Psi, Gamma = Gamma, f = f, it = it, imp = imp)
	output.als.mg.ho2
}