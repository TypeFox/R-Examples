
click.sim <- function(n, int = c(5, 100), alpha, beta = NULL, gamma){

	if (n < 1) stop("Wrong sample size n...\n")
	if (int[1] > int[2]) stop("Lower bound is larger than upper bound in int...\n")
	if (int[1] < 1) stop("Incorrect sequence length in int...\n")

	K <- length(alpha)
	p <- dim(gamma)[1]
	
	S <- list()

	Nk <- rmultinom(1, n, alpha)

	for (k in 1:K){

		for (i in 1:Nk[k]){

			seq.length <- ceiling(runif(1, min = int[1] - 1, max = int[2]))

			if (is.null(beta)){

				new.x <- ceiling(runif(1) * p)

			} else {

				new.x <- which(rmultinom(1, 1, beta[k,]) == 1)

			}

			curr.seq <- new.x

			for (j in 1:(seq.length - 1)){

				new.x <- which(rmultinom(1, 1, gamma[new.x,,k]) == 1)
				curr.seq <- c(curr.seq, new.x)

			}

			S[[length(S) + 1]] <- curr.seq
		}

	}

	return(list(S = S, id = rep(1:K, Nk)))

}




click.read <- function(S){

	p <- max(sapply(S, max))

	n <- length(S)

	X <- array(rep(0, p*p*n), c(p, p, n))
	y <- rep(0, n)

	for (i in 1:n){

		y[i] <- S[[i]][1]

		for (j in 1:(length(S[[i]]) - 1)){

			r <- S[[i]][j]
			c <- S[[i]][j+1]

			X[r,c,i] <- X[r,c,i] + 1

		}

	}

	return(list(X = X, y = y))

}




merge.states <- function(Y, states){

	d <- max(states)

	Z <- matrix(rep(NA, d*d), ncol = d)

	ind1 <- states == d
	Z[d,d] <- sum(Y[ind1,ind1])

	for (i in 1:(d-1)){
		ind1 <- states == i
		Z[i,i] <- sum(Y[ind1,ind1])

		for (j in (i+1):d){
			ind2 <- states == j
			Z[i,j] <- sum(Y[ind1,ind2])
			Z[j,i] <- sum(Y[ind2,ind1])
		}
	}

	return(Z)

}




adjustment.logl <- function(states, X){

	G <- max(states)

	adj <- 0

	for (j in 1:G){

		ind <- states == j
		k <- sum(ind)

		adj <- adj + ifelse(k == 1, 0, -log(k) * sum(X[,ind,]))

	}

#	for (j in 1:G){
#
#		ind <- states == j
#		k <- sum(ind)
#
#		if (k != 1){
#
#			s <- 0
#			for (r in which(ind)){
#				s <- s + sum(y == r)
#			}
#
#			adj <- adj - log(k) * (sum(X[,ind,]) + s)			
#
#		}
#
#	}

	return(adj)

}




state.bic <- function(X, K, states0, new.state, eps, r, iter, min.gamma, scale.const, BIC.flag){

	p <- dim(X)[1]
	n.seq <- dim(X)[3]

	d1 <- max(states0) + 1

	states1 <- states0
	states1[new.state] <- d1

	M1 <- K * d1 * (d1 - 1) + K - 1

	Z1 <- apply(X, c(3), merge.states, states1)
	Z1 <- array(as.vector(Z1), c(d1, d1, n.seq))

	Q1 <- click.EM(X = Z1, K = K, eps = eps, r = r, iter = iter, min.gamma = min.gamma, scale.const = scale.const)
	adjust1 <- adjustment.logl(states1, X)

	L1 <- Q1$logl + adjust1 - n.seq * log(p) + n.seq * log(d1)

	if (BIC.flag){
		return(list(logl1 = L1, BIC1 = -2 * L1 + M1 * log(n.seq)))
	} else {
		return(list(logl1 = L1, AIC1 = -2 * L1 + 2 * M1))
	}
}



rearrange.states <- function(X, K, best.ll, best.states, eps, r, iter, min.gamma, scale.const, BIC.flag, silent){

	p <- dim(X)[1]
	n.seq <- dim(X)[3]

	states0 <- best.states

	#########################################
	# trying to move states to other groups of states
	if (!silent) cat("\nSTAGE B: Rearranging states...\n\n")
	#########################################

	d.max <- max(states0)

	repeat{

		improv.flag <- 0

		for (k in 1:d.max){

			ind1 <- which(states0 == k)

			if (length(ind1) > 1){

				for (i in ind1){
	
					for (j in (1:d.max)[-k]){

						states <- states0
						states[i] <- j
	
						X1 <- apply(X, c(3), merge.states, states)
						X1 <- array(as.vector(X1), c(d.max, d.max, n.seq))

						Q <- click.EM(X = X1, K = K, eps = eps, r = r, iter = iter, min.gamma = min.gamma, scale.const = scale.const)
	
						adjust <- adjustment.logl(states, X)

						ll <- Q$logl + adjust - n.seq * log(p) + n.seq * log(d.max)

						if (!silent) cat("\tStates:", states, "\tlogL =", ll, "\n")

						if (!is.nan(ll)){	
							if (ll > best.ll){
								best.ll <- ll
								best.states <- states
								improv.flag <- 1
							}
						}
					}

				}

			}

		}

		d <- max(best.states)
		M <- K * d * (d - 1) + K - 1

		if (BIC.flag){
			best.BIC <- -2 * best.ll + M * log(n.seq)
			if (!silent) cat("\tCurrent set of states:", best.states, "logL =", best.ll, "BIC =", best.BIC, "\n\n")
		} else {
			best.AIC <- -2 * best.ll + M * 2
			if (!silent) cat("\tCurrent set of states:", best.states, "logL =", best.ll, "AIC =", best.AIC, "\n\n")
		}
		states0 <- best.states

		if (improv.flag == 0){
			if (!silent) cat("\tStopped rearrangement...\n\n")
			break
		}

	}

	if (BIC.flag){
		return(list(states = best.states, ll = best.ll, BIC = best.BIC))
	} else {
		return(list(states = best.states, ll = best.ll, AIC = best.AIC))
	}

}


separate.states.BIC <- function(X, K, best.ll, best.states, eps, r, iter, min.gamma, scale.const, BIC.flag, silent){

	p <- dim(X)[1]
	n.seq <- dim(X)[3]

	d <- max(best.states)

	M <- K * d * (d - 1) + K - 1
	
	if (BIC.flag){
		best.BIC <- -2 * best.ll + M * log(n.seq)
		curr.BIC <- best.BIC
	} else {
		best.AIC <- -2 * best.ll + M * 2
		curr.AIC <- best.AIC
	}
		
	states0 <- best.states

	###################################################
	# checking the split of the state
	if (!silent) cat("\nSTAGE A: Separating states...\n\n")
	###################################################

	d.max <- max(states0)


	for (k in 1:d.max){

		ind1 <- which(states0 == k)

		if (length(ind1) > 1){
	
			for (j in ind1){

				Q <- state.bic(X = X, K = K, states0 = states0, new.state = j, eps = eps, r = r, iter = iter, min.gamma = min.gamma, scale.const = scale.const, BIC.flag = BIC.flag)
				states <- states0
				states[j] <- d.max + 1
				if (!silent) cat("\tStates:", states, "\tlogL =", Q$logl1, "\n")	

				if (BIC.flag){ 			
					if (!is.nan(Q$BIC1)){
						if (Q$BIC1 < best.BIC){
							best.ll <- Q$logl1
							best.states <- states
							best.BIC <- Q$BIC1
						}
					}
				} else {
					if (!is.nan(Q$AIC1)){
						if (Q$AIC1 < best.AIC){
							best.ll <- Q$logl1
							best.states <- states
							best.AIC <- Q$AIC1
						}
					}

				}

	
			}

		}

	}

	if (BIC.flag){
		if (!silent){
			if (best.BIC == curr.BIC){
				cat("\tNo improvement reached...\n\n")
			} else {
				cat("\tCurrent set of states:", best.states, "logL =", best.ll, "BIC =", best.BIC, "\n\n")
				states0 <- best.states
			}
		}
		return(list(states = best.states, ll = best.ll, BIC = best.BIC))

	} else {
		if (!silent){
			if (best.AIC == curr.AIC){
				cat("\tNo improvement reached...\n\n")
			} else {
				cat("\tCurrent set of states:", best.states, "logL =", best.ll, "AIC =", best.AIC, "\n\n")
				states0 <- best.states
			}
		}
		return(list(states = best.states, ll = best.ll, AIC = best.AIC))

	}

}




combine.states.BIC <- function(X, K, best.ll, best.states, eps, r, iter, min.gamma, scale.const, BIC.flag, silent){


	p <- dim(X)[1]
	n.seq <- dim(X)[3]

	d <- max(best.states)
	M <- K * d * (d - 1) + K - 1
	if (BIC.flag){
		best.BIC <- -2 * best.ll + M * log(n.seq)
		curr.BIC <- best.BIC
	} else {
		best.AIC <- -2 * best.ll + M * 2
		curr.AIC <- best.AIC
	}
	d0 <- d - 1
	states <- best.states

	if (!silent) cat("\nSTAGE A: Merging equivalence blocks...\n\n")

	for (i in 1:(d-1)){

		for (j in (i+1):d){

			states0 <- states
			states0[states0 == j] <- i

			states0 <- renumerate.seq(states0)

			if (d0 != 1){
				X1 <- apply(X, c(3), merge.states, states0)
			} else {
				X1 <- apply(X, c(3), sum)
			}

			X1 <- array(as.vector(X1), c(d0, d0, n.seq))

			Q <- click.EM(X = X1, K = K, eps = eps, r = r, iter = iter, min.gamma = min.gamma, scale.const = scale.const)

			adjust <- adjustment.logl(states0, X)

			ll <- Q$logl + adjust - n.seq * log(p) + n.seq * log(d0)
		 	if (!silent) cat("\tStates:", states0, "\tlogL =", ll, "\n")

			M <- K * d0 * (d0 - 1) + K - 1

			if (BIC.flag){
				curr.BIC <- -2 * ll + M * log(n.seq)
			} else {
				curr.AIC <- -2 * ll + M * 2
			}


			if (BIC.flag){ 			
				if (!is.nan(curr.BIC)){
					if (curr.BIC < best.BIC){
						best.ll <- ll
						best.states <- states0
						best.BIC <- curr.BIC
					}
				}
			} else {
				if (!is.nan(curr.AIC)){
					if (curr.AIC < best.AIC){
						best.ll <- ll
						best.states <- states0
						best.AIC <- curr.AIC
					}
				}
			}


		}

	}


	if (BIC.flag){
		if (!silent){
			if (best.BIC == curr.BIC){
				cat("\tNo improvement reached...\n\n")
			} else {
				cat("\tCurrent set of states:", best.states, "logL =", best.ll, "BIC =", best.BIC, "\n\n")
			}
		}
	
		return(list(states = best.states, ll = best.ll, BIC = best.BIC))

	} else {
		if (!silent){
			if (best.AIC == curr.AIC){
				cat("\tNo improvement reached...\n\n")
			} else {
				cat("\tCurrent set of states:", best.states, "logL =", best.ll, "AIC =", best.AIC, "\n\n")
			}
		}

		return(list(states = best.states, ll = best.ll, AIC = best.AIC))

	}


}



renumerate.seq <- function(x){

	a <- table(x)
	a <- as.numeric(names(a))
	y <- x

	for (i in 1:length(a)){
		ind <- x == a[i]
		y[ind] <- i
	}

	return(y)

}




click.forward <- function(X, K, eps = 1e-10, r = 100, iter = 5, bic = TRUE, min.gamma = 1e-3, scale.const = 1.0, silent = FALSE){

	if (K < 1) stop("Wrong number of mixture components K...\n")
	if (eps <= 0) stop("Wrong value of eps...\n")
	if (r < 1) stop("Wrong number of random restarts r...\n")
	if (iter < 1) stop("Wrong number of iterations iter...\n")
	if (min.gamma < 0) stop("Wrong lower bound min.gamma...\n")
	if (scale.const <= 0) stop("Wrong value of scale.const...\n")

	p <- dim(X)[1]
	n.seq <- dim(X)[3]

	best.states <- rep(1, p)
	best.ll <- -Inf ### -(sum(X) + n.seq) * log(p)


	repeat{

		B <- separate.states.BIC(X = X, K = K, best.ll = best.ll, best.states = best.states, eps = eps, r = r, iter = iter, min.gamma = min.gamma, scale.const = scale.const, BIC.flag = bic, silent)

		if (sum(B$states != best.states) == 0){
			break
		}

		best.states <- B$states
		best.ll <- B$ll


		A <- rearrange.states(X = X, K = K, best.ll = best.ll, best.states = best.states, eps = eps, r = r, iter = iter, min.gamma = min.gamma, scale.const = scale.const, BIC.flag = bic, silent)

		best.states <- A$states
		best.ll <- A$ll
		if (bic){
			best.BIC <- A$BIC
		} else {
			best.AIC <- A$AIC
		}

	}

	# the following part can be completely avoided
	# if we save all necessary parameters right away

	d.max <- max(best.states)

	X1 <- apply(X, c(3), merge.states, best.states)
	X1 <- array(as.vector(X1), c(d.max, d.max, n.seq))

	Q <- click.EM(X = X1, K = K, eps = eps, r = r, iter = iter, min.gamma = min.gamma, scale.const = scale.const)

	if (bic){
		ret <- list(z = Q$z, id = Q$id, alpha = Q$alpha, gamma = Q$gamma, states = best.states, logl = best.ll, BIC = best.BIC)
	} else {
		ret <- list(z = Q$z, id = Q$id, alpha = Q$alpha, gamma = Q$gamma, states = best.states, logl = best.ll, AIC = best.AIC)
	}

	class(ret) <- "search"
	return(ret)

}






click.backward <- function(X, K, eps = 1e-10, r = 100, iter = 5, bic = TRUE, min.gamma = 1e-3, scale.const = 1.0, silent = FALSE){

	if (K < 1) stop("Wrong number of mixture components K...\n")
	if (eps <= 0) stop("Wrong value of eps...\n")
	if (r < 1) stop("Wrong number of random restarts r...\n")
	if (iter < 1) stop("Wrong number of iterations iter...\n")
	if (min.gamma < 0) stop("Wrong lower bound min.gamma...\n")
	if (scale.const <= 0) stop("Wrong value of scale.const...\n")

	p <- dim(X)[1]
	n.seq <- dim(X)[3]

	Q <- click.EM(X = X, K = K, eps = eps, r = r, iter = iter, min.gamma = min.gamma, scale.const = scale.const)

	best.states <- 1:p
	best.ll <- Q$logl
	if (bic){
		best.BIC <- Q$BIC
	} else {
		best.AIC <- -2 * best.ll + 2 * (K * p * (p - 1) + K - 1)
	}

	
	repeat{

		if (max(best.states) == 2) break

		B <- combine.states.BIC(X = X, K = K, best.ll = best.ll, best.states = best.states, eps = eps, r = r, iter = iter, min.gamma = min.gamma, scale.const = scale.const, BIC.flag = bic, silent = silent)

		if (sum(B$states != best.states) == 0){
			break
		}

		best.states <- B$states
		best.ll <- B$ll


		A <- rearrange.states(X = X, K = K, best.ll = best.ll, best.states = best.states, eps = eps, r = r, iter = iter, min.gamma = min.gamma, scale.const = scale.const, BIC.flag = bic, silent = silent)
		best.states <- A$states
		best.ll <- A$ll
		if (bic){
			best.BIC <- A$BIC
		} else {
			best.AIC <- A$AIC
		}

	}

	# the following part can be completely avoided
	# if we save all necessary parameters right away

	d.max <- max(best.states)

	X1 <- apply(X, c(3), merge.states, best.states)
	X1 <- array(as.vector(X1), c(d.max, d.max, n.seq))

	Q <- click.EM(X = X1, K = K, eps = eps, r = r, iter = iter, min.gamma = min.gamma, scale.const = scale.const)

	if (bic){
		ret <- list(z = Q$z, id = Q$id, alpha = Q$alpha, gamma = Q$gamma, states = best.states, logl = best.ll, BIC = best.BIC)
	} else {
		ret <- list(z = Q$z, id = Q$id, alpha = Q$alpha, gamma = Q$gamma, states = best.states, logl = best.ll, AIC = best.AIC)
	}

	class(ret) <- "search"
	return(ret)

}



click.var <- function(X, y = NULL, alpha, beta = NULL, gamma, z){

	n <- dim(X)[3]
	p <- dim(X)[1]
	K <- length(alpha)

	if (is.null(beta)){
		M <- K - 1 + K * p * (p - 1)
	} else {
		M <- K - 1 + K * (p - 1) + K * p * (p - 1)
	}

	Inf.mat <- NULL

	for (i in 1:n){

		v <- z[i,1:(K-1)] / alpha[1:(K-1)] - z[i,K] / alpha[K]

		if (!is.null(beta) && !is.null(y)){

			for (k in 1:K){

				v1 <- rep(0, p-1)

				if (y[i] != p){
					v1[y[i]] <- z[i,k] / beta[k,y[i]]
				} else {
					v1[1:(p-1)] <- -z[i,k] / beta[k,y[i]]
				}

				v <- c(v, v1)

			}

		}

		
		for (k in 1:K){

			for (j in 1:p){

				v1 <- z[i,k] * (X[j,1:(p-1),i] / gamma[j,1:(p-1),k] - X[j,p,i] / gamma[j,p,k])
				v <- c(v, v1)

			}

		}

		
		Inf.mat <- cbind(Inf.mat, v)

	}

	Inf.mat <- Inf.mat %*% t(Inf.mat)

	V <- solve(Inf.mat)

	return(V)

}




click.predict <- function(M = 1, gamma, pr = NULL){

	p <- dim(gamma)[1]
	K <- dim(gamma)[3]

	if (is.null(pr)) pr <- rep(1 / K, K)

	if (M < 1) stop("Wrong number of transitions M...\n")
	if (length(pr) != K) stop("Incorrect length of probability vector pr...\n")

	state.prob <- matrix(0, p, p)
	Pr <- gamma

	for (k in 1:K) Pr[,,k] <- diag(p)
	for (k in 1:K){
		for (m in 1:M){
			Pr[,,k] <- Pr[,,k] %*% gamma[,,k]
		}
		state.prob <- state.prob + pr[k] * Pr[,,k]
	}

	return(state.prob)

}



















