mlm.influence <-
function (model, do.coef = TRUE, m=1, ...) 
{
#    wt.res <- weighted.residuals(model)
#    e <- na.omit(wt.res)
#    if (model$rank == 0) {
#        n <- length(wt.res)
#        sigma <- sqrt(deviance(model)/df.residual(model))
#        res <- list(hat = rep(0, n), coefficients = matrix(0, 
#            n, 0), sigma = rep(sigma, n), wt.res = e)
#    return(res)
#    }

	# helper functions
  vec <- function(M) {
  	R <- matrix(M, ncol=1)
  	if (is.vector(M)) return(R)
  	nn<-expand.grid(dimnames(M))[,2:1]
  	rownames(R) <- apply(as.matrix(nn), 1, paste, collapse=":")
  	R
  	}

	X <- model.matrix(model)
	data <- model.frame(model)
	Y <- as.matrix(model.response(data))
	r <- ncol(Y)
	n <- nrow(X)
	p <- ncol(X)
	labels <- rownames(X)
	call <- model$call
	
	B <- coef(model)
	E <- residuals(model)
	XPXI <- solve(crossprod(X))
	EPEI <- solve(crossprod(E))
	vB <- vec(t(B))
	S <- crossprod(E)/(n-p);
	V <- solve(p*S) %x% crossprod(X)
	
	subsets <- t(combn(n, m))
	nsub <- nrow(subsets) 
	
	# at this point, there are several choices for fitting:
	# (a) update(model, subset=!(1:n) %in% subsets[i,])
	# (b) lm.fit(X[rows,], Y[rows,])
	# (c) direct matrix calculation
	#
	# For each subset:  keep Beta=coefficients, Hat=hatvalues E=residuals, ...
	# I don't know how to use the qr() components of (a) or (b) to
	# calculate hatvalues, so I'm using direct calculation

	Beta <- as.list(rep(0, nsub))
	R <- L <- H <- Q <- as.list(rep(0, nsub))
	CookD <- as.vector(rep(0, nsub))

	for (i in seq(nsub)) {
		I <- c(subsets[i,])
		rows <- which(!(1:n) %in% I)
		XI <- X[rows,]
		YI <- Y[rows,]
		BI <- solve(crossprod(XI)) %*% t(XI) %*% YI
		EI <- (Y - X %*% BI)[I, , drop=FALSE]
		CookD[i] <- t(vec(B - BI)) %*% V %*% vec(B - BI)
		H[[i]] <- X[I, , drop=FALSE] %*% XPXI %*% t(X[I, , drop=FALSE])
		Q[[i]] <- EI %*% EPEI %*% t(EI)
		if (do.coef) Beta[[i]] <- BI
		L[[i]] <- if(m==1) H[[i]] / (1-H[[i]]) else H[[i]] %*% solve(diag(m) - H[[i]])
		R[[i]] <- if(m==1) Q[[i]] / (1-H[[i]]) 
					else {
					IH <- mpower(diag(m)-H[[i]], -1/2)
					IH %*% Q[[i]] %*% IH
					}
	}
	if(m==1) {
		H <- unlist(H)
		Q <- unlist(Q)
		L <- unlist(L)
		R <- unlist(R)
		subsets <- c(subsets)
	}
	result <- list(m=m, H=H, Q=Q, CookD=CookD, L=L, R=R, subsets=subsets, labels=labels, call=call)
	if (do.coef) result <- c(result, list(Beta=Beta))
	class(result) <-"inflmlm"
	result
}
