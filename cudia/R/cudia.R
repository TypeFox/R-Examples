cudia <- function(formula,data,K, ...) UseMethod("cudia")

cudiaDet<- function(formula,data,K){
	
	# model frame (input processing)
	cudia.mf<-model.frame(formula=formula,data=data)
	s <- model.response(cudia.mf,type="numeric")
	s.mean <- mean(s)
	s.centered <- (s - s.mean)
	x.label <- attr(attr(cudia.mf,"terms"),"term.labels")
	x.dim <- length(x.label)
	x <- (as.matrix(cudia.mf[,x.label]))
	x.sd <- apply(x,2,sd)
	x.mean <- colMeans(x)
	x <- t((t(x)-x.mean)/x.sd) # standardization
	div <- unique(s.centered)
	div.dim <- length(div)

	# model parameter initialization
	alpha <- matrix(1,K)
	pipk <- matrix(rdirichlet(div.dim,alpha),div.dim,K,dimnames=list(1:div.dim,1:K))
	theta <- matrix(rnorm(x.dim*K,sd=0.001),K,x.dim,dimnames=list(1:K,x.label))
	eta <- matrix(rnorm(K,sd=0.001),K,1,dimnames=list(1:K,"Aggr"))
	N <- nrow(x)
	z <- matrix(1/K,N,K)
	eps <- .Machine$double.eps
	cudia_imp.pre <- s.centered
	theta.pre <- theta
	eta.pre <- eta
	Np <- matrix(1,div.dim,1)
	for(div.idx in 1:div.dim){ Np[div.idx] <- length(which(s.centered==div[div.idx]))}

	# EM steps
    max_itr <- 250
	for (em_cnt in 1:max_itr){
		# E-step
		for (div.idx in 1:div.dim){
			aggr <- div[div.idx]
			sub.idx <- which(s.centered==aggr)
			x.sub <- x[sub.idx,]
			z.sub <- z[sub.idx,]
			# deterministic distance calculation
			for(k in 1:K){
				dist <- (t(x.sub) - theta[k,])**2
				if(x.dim>1){ dist <- colSums(dist) }
				adj <- 2*(aggr - pipk[div.idx,] %*% eta)*eta[k]/K/Np[div.idx]
				z.sub[,k] <- -dist+rep(adj,Np[div.idx])
			}
			# membership-assignment (z)
			for(row in 1:nrow(z.sub)){
				membership <- which.max(z.sub[row,])
				z.sub[row,] <- matrix(0,1,K)
				z.sub[row,membership] <- 1
			}
			z[sub.idx,]<- z.sub
			pipk[div.idx,] <- (colSums(z.sub) + 1.0)/(Np[div.idx]+K)
		}
	
		# M-step
		Nk <- (colSums(z) + 1.0)
		# theta estimation
		for (k in 1:K){	theta[k,] <- colSums(matrix(rep((z[,k]),x.dim),N,x.dim) * x / (Nk[k])) }
		# eta estimation: Ridge-regression for numerical stability
		eta <-matrix(lm.ridge(div~pipk-1,lambda=0.01)$coef,K,1,dimnames=list(1:K,"Aggr"))
	
		# convergence test using the estimated parameters: theta (1), eta (2)
		conv.t1 <- norm(theta.pre-theta)/K/x.dim
		conv.t2 <- norm(eta.pre-eta)/K
		if (conv.t1 < 1e-6 && conv.t2 < 1e-6){ cat("\n",em_cnt," iterations\n\n"); break;}	
		if (em_cnt==max_itr){ cat("\n",em_cnt," iterations\n\n") }
		theta.pre <- theta
		eta.pre <- eta
	}

	# pipk smoothing
	w.pri <- 0.8
	w.obs <- 0.2
	prior <- matrix(1/K,1,K)
	for (div.idx in 1:div.dim){
		posterior <- w.obs*pipk[div.idx,] + w.pri*prior
		pipk[div.idx,] <- posterior
	}

	# probabilistic z-estimation using Breg-div.
	for (div.idx in 1:div.dim){
		aggr <- div[div.idx]
		sub.idx <- which(s.centered==aggr)
		x.sub <- x[sub.idx,]
		z.sub <- z[sub.idx,]
		for(k in 1:K){
			breg.div <- -(t(x.sub) - theta[k,])**2
			if(x.dim>1){ breg.div <- colSums(breg.div) }
			z.sub[,k] <- pipk[div.idx,k] * exp(breg.div)
		}
		z[sub.idx,]<- z.sub/rowSums(z.sub)
	}

	# cudia imputation formula
	cudia_imp <- z %*% eta + s.mean
	
	# output formation
	Nk <- colMeans(z)
	x <- t(t(x)*x.sd + x.mean) # de-standardizing
	theta <- t(t(theta) * x.sd + x.mean) # de-standardizing
	
	list(indiv=x,fitted.values=cudia_imp,theta=theta,eta=eta,Nk=Nk,xlab=x.label)
}

cudia.default <- function(formula,data,K,...)
{
	
	est <- cudiaDet(formula,data,K)
	est$call <- match.call()
	class(est) <- "cudia"
	est
}

print.cudia <- function(x, ...)
{
	cat("Call:\n")
	print(x$call)
	cat("\n\nIndiv-level Cluster Coefficients:\n")
	print(x$theta)
	cat("\n\nAggr-level Cluster Coefficients:\n")
	print(x$eta)
	cat("\n\nCluster Sizes:\n")
	print(x$Nk)	
}

summary.cudia <- function(object, ...)
{
	cat("Call:\n")
	print(object$call)
	cat("\n\nIndiv-level Cluster Coefficients:\n")
	print(object$theta)
	cat("\n\nAggr-level Cluster Coefficients:\n")
	print(object$eta)
	cat("\n\nCluster Sizes:\n")
	print(object$Nk)	
}

plot.cudia <- function(x, ...)
{
	indiv.col <- ncol(x$indiv)
	par(mfrow=c(1,indiv.col))
	for (i in 1:indiv.col){
		plot(x$indiv[,i],x$fitted.values,xlab=x$xlab[i],ylab="CUDIA-imputed")	
	}
}


