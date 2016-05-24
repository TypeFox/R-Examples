#functions to include in FLAM R package
#updated 6/16/14

#########################################################################################################
### FUNCTIONS FOR DATA MANAGEMENT
#########################################################################################################

#########################################################################################################
#functions that take in matrix with columns being covariates (x) and returns a list of matrices,
#with each matrix P_j in the list corresponding to the permutation matrix for (i.e. P_j %*% x_j "sorts" x_j from least to greatest)

makePerm = function(x) {
	
	n = length(x)
	perm.mat = t(sapply(order(x), function(entry, n) {vec = rep(0,n); vec[entry]=1; return(vec)}, n=n))
	return(list(perm.mat))
}

makePermList = function(x) {
	
	perm.list = sapply(1:ncol(x),function(i,x) makePerm(x[,i]),x=x)
	return(perm.list)	
}

#########################################################################################################
#functions to convert f matrix to theta matrix and vice versa 

f.to.theta = function(f.mat, x, rank.x=NULL) {
	
	theta.mat = f.mat
	if (is.null(rank.x)) {
    for (j in 1:ncol(f.mat)) theta.mat[,j] = f.mat[rank(x[,j]),j]
    } else {
		theta.mat  = sapply(1:ncol(f.mat),function(index,rank.x,f.mat) f.mat[rank.x[,index],index], rank.x=rank.x, f.mat=f.mat)
	}
	
	return(theta.mat)
	
}

theta.to.f = function(theta.mat, x, order.x=NULL) {
	
	f.mat = theta.mat
	if (is.null(order.x)) for (j in 1:ncol(theta.mat)) f.mat[,j] = theta.mat[order(x[,j]),j] else {
		f.mat  = sapply(1:ncol(f.mat),function(index,order.x,theta.mat) theta.mat[order.x[,index],index], order.x=order.x, theta.mat=theta.mat)
	}
	
	return(f.mat)
	
}

#########################################################################################################
### END OF FUNCTIONS FOR DATA MANAGEMENT
#########################################################################################################


#########################################################################################################
### FUNCTIONS FOR FLAM
#########################################################################################################

#########################################################################################################
#main FLAM function: takes in the response (y), matrix with columns being covariates (x), 
#sequence of decreasing lambdas (desired penalty), and sequence of mixing parameter (alpha)
#family ("gaussian" for squared error loss and "binomial" for logistic loss)
#and method for "BCD"=block coordinate descent, "GGD"=generalized gradient descent, "GGD.backtrack"=GGD with backtracking (varying step size)
#method argument is currently ignored for family="binomial"
#tolerance is the convergence criterion for the objective

flam = function(x, y, lambda.min.ratio = 0.01, n.lambda = 50, lambda.seq = NULL, alpha.seq = 1, family = "gaussian", method = "BCD", tolerance = 10e-6) {	
	version = 3
	call = match.call()
	fit = list()
	
	if (is.null(nrow(x))) x = matrix(x,ncol=1)
	p = ncol(x); n = nrow(x)
	
	if (is.null(lambda.seq)) {
		max.lam = max(sapply(alpha.seq, maxLambda, y=y, x=x))
		lambda.seq = exp(seq(log(max.lam),log(max.lam*lambda.min.ratio),len=n.lambda))
	}
	
	#checks
	if (length(y)!=n) stop("The length of 'y' must equal the number of rows of 'x'")
	if (length(lambda.seq)==1) stop("Provide a sequence of decreasing values for 'lambda.seq'")
	if (!(family %in% c("gaussian","binomial"))) stop("'family' not recognized - must be 'gaussian' or 'binomial'")
	if (min(lambda.seq)<=0) stop("Values in 'lambda.seq' must be positive")
	if (min(alpha.seq)<0 | max(alpha.seq)>1) stop("Values in 'alpha.seq' must be in [0,1]")
	if (!(method %in% c("BCD","GGD","GGD.backtrack"))) stop("'method' not recognized - must be 'BCD', 'GGD', or 'GGD.backtrack'")

	#make sure lambda.seq is decreasing
	lambda.seq = sort(lambda.seq, decreasing=TRUE)
	#make alpha.seq increasing
	alpha.seq = sort(alpha.seq, decreasing=FALSE)

	rank.x = apply(x,2,rank)
	order.x = apply(x,2,order)
	
	#matrices and lists to store results
	fit$all.alpha = rep(alpha.seq, each=length(lambda.seq))
	fit$all.lambda = rep(lambda.seq, length(alpha.seq))
	n.tuning = length(fit$all.lambda)
	fit$f.hat.list <- fit$theta.hat.list <- vector("list", n.tuning)
	fit$beta0.hat.vec <- n.iter.vec <- L.vec <- rep(NA, n.tuning)
	fit$y.hat.mat = matrix(NA, nrow=n.tuning, ncol=n)
	obj.list = vector("list",n.tuning)
	
	for (l in 1:n.tuning) {
		
		if (l==1) initial.f.mat=NULL else if ((l %% length(lambda.seq))==1) initial.f.mat=fit$f.hat.list[[l-length(lambda.seq)]] else initial.f.mat=fit$f.hat.list[[l-1]]
		
		if (family=="gaussian" & method=="BCD") {		
			one.fit = flam.helper(y=y, lambda=fit$all.lambda[l], alpha=fit$all.alpha[l], x=x, initial.f.mat=initial.f.mat, n=n, p=p, order.x=order.x, rank.x=rank.x, tolerance=tolerance)			
		}
		
		if (family=="gaussian" & method=="GGD") {
			one.fit = flam.ggd.helper(y=y, lambda=fit$all.lambda[l], alpha=fit$all.alpha[l], x=x, rank.x=rank.x, order.x=order.x, n=n, p=p, initial.f.mat=initial.f.mat, tolerance=tolerance)
		}
		
		if (family=="gaussian" & method=="GGD.backtrack") {
			if (l==1) initial.L=1 else if ((l %% length(lambda.seq))==1) initial.L=L.vec[l-length(lambda.seq)] else initial.L=L.vec[l-1]
			one.fit = flam.ggd.helper_backtrack(y=y, lambda=fit$all.lambda[l], alpha=fit$all.alpha[l], x=x, n=n, p=p, initial.L=initial.L, rank.x=rank.x, order.x=order.x, initial.f.mat=initial.f.mat, tolerance=tolerance)
		}
		
		if (family=="binomial") one.fit = flam.logistic.helper(y=y, lambda=fit$all.lambda[l], alpha=fit$all.alpha[l], x=x, n=n, p=p, rank.x=rank.x, order.x=order.x, initial.f.mat=initial.f.mat, tolerance=tolerance)
		
		fit$f.hat.list[[l]] = one.fit$f.hat.mat
		fit$theta.hat.list[[l]] = one.fit$theta.hat.mat
		fit$beta0.hat.vec[l] = one.fit$beta0.hat
		fit$y.hat.mat[l,] = one.fit$y.hat
		if (method=="GGD.backtrack") L.vec[l] = one.fit$L
  
	}	

	fit$non.sparse.list = sapply(fit$theta.hat.list, function(theta.hat.mat) which(apply(theta.hat.mat, 2, function(col) sum(col^2))!=0))
	fit$num.non.sparse = sapply(fit$theta.hat.list, function(theta.hat.mat) sum(apply(theta.hat.mat, 2, function(col) sum(col^2)!=0)))

	fit$y = y; fit$x = x; fit$family = family; fit$method = method; fit$call = call; fit$tolerance = tolerance
	class(fit) = "flam"
	return(fit)
}

#########################################################################################################
# helper functions that are called by main FLAM function above

#########################################################################################################
# squared error loss, block coordinate descent
#########################################################################################################

flam.helper = function(y,lambda,alpha,x,initial.theta.mat=NULL,initial.f.mat=NULL,n,p,order.x,rank.x,tolerance) {
  
  #intialize matrix for storing theta hat estimates
  if (!is.null(initial.f.mat)) {
    initial.theta.mat = f.to.theta(initial.f.mat, x, rank.x)
  } else if (is.null(initial.theta.mat)) {
    initial.theta.mat = matrix(0, nrow=n, ncol=p)
  }
  beta0 = 0
  theta.mat = flamstep(initial.theta.mat,y,lambda,alpha,n,p,order.x,rank.x,beta0,tolerance)
  f.hat.mat = theta.to.f(theta.mat, x, order.x)

  return(list(f.hat.mat=f.hat.mat,theta.hat.mat=theta.mat,beta0.hat=beta0,y.hat=(beta0+rowSums(theta.mat)),x=x, family="gaussian", method="BCD"))
}

#########################################################################################################
# squared error loss, generalized gradient descent
#########################################################################################################

flam.ggd.helper = function(y,lambda,alpha,x,initial.theta.mat=NULL,initial.f.mat=NULL,rank.x, order.x,n,p,tolerance) {
		
	#intialize matrix for storing theta hat estimates
	if (!is.null(initial.theta.mat)) {
		theta.mat = initial.theta.mat
		} else if (!is.null(initial.f.mat)) {
			theta.mat = f.to.theta(initial.f.mat, x,rank.x)
		} else {
			theta.mat = matrix(0, nrow=nrow(x), ncol=p)
		}
		
	#matrix for storing theta hat estimates from previous iteration
	theta.mat_old = theta.mat
	beta0 = sum(y)/n
	
	L = p #Lipschitz constant

	#vector to store values of objective so convergence can be checked
	obj.vec = c()
	
	#number of iterations counter
	n.iter = 0

	converge = FALSE
	
	while (converge==FALSE & n.iter<1000) {
		
		n.iter = n.iter + 1
		
		resid.mat = theta.mat_old + matrix((y - beta0 - rowSums(theta.mat_old))/L,ncol=1) %*% matrix(1, nrow=1, ncol=p)	
		theta.mat = sapply(1:p, fit.est, resid.mat=resid.mat, tuning.parameter=alpha*lambda/L, rank.x=rank.x, order.x=order.x)
			
		theta.mat = apply(theta.mat, 2, soft.scale.est, tuning.parameter=(1-alpha)*lambda/L)
		
		#store value of objective
		theta.mat.reord = sapply(1:p, function(theta.mat, order.x, i) theta.mat[,i][order.x[,i]], theta.mat=theta.mat, order.x=order.x)
		sum1 = 0
		for (k in 2:n) sum1 = sum1 + sum(abs(theta.mat.reord[k,] - theta.mat.reord[k-1,]))

		sum2 = sum(sqrt(colSums(theta.mat^2)))
		
		if (p>1) sq.err = sum((y - (beta0 + rowSums(theta.mat)))^2) else sq.err = sum((y - (beta0 + theta.mat))^2)
		new.obj = 0.5 * sq.err + lambda * alpha * sum1 + lambda * (1-alpha) * sum2
		obj.vec = c(obj.vec, new.obj)

		#check for convergence
		if (length(obj.vec)>1) {
			dev = abs(obj.vec[length(obj.vec)] - obj.vec[length(obj.vec)-1]) / abs(obj.vec[length(obj.vec)-1])
			if (dev < tolerance) converge = TRUE
		}		
		theta.mat_old = theta.mat
	}
	
	#center theta_j's
	theta.mat = scale(theta.mat, center=T, scale=F)
	
	f.hat.mat = theta.to.f(theta.mat, x, order.x)
	
	return(list(f.hat.mat=f.hat.mat,theta.hat.mat=theta.mat,beta0.hat=beta0,y.hat=(beta0+rowSums(theta.mat)),objective=obj.vec,n.iter=n.iter,x=x, family="gaussian", method="GGD"))
}

#########################################################################################################
# squared error loss, generalized gradient descent with backtracking
#########################################################################################################

flam.ggd.helper_backtrack = function(y,lambda,alpha,x,initial.L,initial.theta.mat=NULL,initial.f.mat=NULL,n,p, rank.x, order.x, tolerance) {
		
	#intialize matrix for storing theta hat estimates
	if (!is.null(initial.theta.mat)) {
		theta.mat = initial.theta.mat
		} else if (!is.null(initial.f.mat)) {
			theta.mat = f.to.theta(initial.f.mat, x)
		} else {
			theta.mat = matrix(0, nrow=nrow(x), ncol=p)
		}
		
	#matrix for storing theta hat estimates from previous iteration
	theta.mat_old = theta.mat
	beta0 = sum(y)/n
	
	#vector to store values of objective so convergence can be checked
	obj.vec = c()
	
	#store values of L checked and whether a step was taken
	L.seq = initial.L
	maximizer = c()
	
	#number of iterations counter
	n.iter = 0

	converge = FALSE
	
	while (converge==FALSE & n.iter<1000) {
		
		n.iter = n.iter + 1
		
		#check if step should be taken by solving for theta.mat_star
		resid.mat = theta.mat_old + matrix((y - beta0 - rowSums(theta.mat_old))/L.seq[n.iter],ncol=1) %*% matrix(1, nrow=1, ncol=p)		
		theta.mat_star = sapply(1:p, fit.est, resid.mat=resid.mat, tuning.parameter=alpha*lambda/L.seq[n.iter], rank.x=rank.x, order.x=order.x)
		
		theta.mat_star = apply(theta.mat_star, 2, soft.scale.est, tuning.parameter=(1-alpha)*lambda/L.seq[n.iter])
		
		l_true = calc.sq.err(theta.mat_star, y-beta0)
		
		l_approx = calc.sq.err(theta.mat_old, y-beta0) + sum(diag(t(theta.mat_star-theta.mat_old) %*% (-beta0 -y + rowSums(theta.mat_old)) %*% matrix(1,nrow=1,ncol=p))) + sum((theta.mat_star-theta.mat_old)^2)*L.seq[n.iter]/2
		
		if (l_approx>=l_true) {
			
			maximizer = c(maximizer,1)
						
			#store value of objective
			theta.mat.reord = sapply(1:p, function(theta.mat_star, order.x, i) theta.mat_star[,i][order.x[,i]], theta.mat_star=theta.mat_star, order.x=order.x)
			sum1 = 0
			for (k in 2:n) sum1 = sum1 + sum(abs(theta.mat.reord[k,] - theta.mat.reord[k-1,]))
			sum2 = sum(sqrt(colSums(theta.mat_star^2)))
			new.obj = l_true + lambda * alpha * sum1 + lambda * (1-alpha) * sum2
			
			theta.mat = theta.mat_star
		
			obj.vec = c(obj.vec, new.obj)

			#check for convergence
			if (length(obj.vec)>1) {
				dev = abs(obj.vec[length(obj.vec)] - obj.vec[length(obj.vec)-1]) / abs(obj.vec[length(obj.vec)-1])
				if (dev < tolerance) converge = TRUE
			}
		
			theta.mat_old = theta.mat
			L.seq = c(L.seq, L.seq[n.iter])
						
		} else {
			
			maximizer = c(maximizer,0)
			L.seq = c(L.seq, L.seq[n.iter]*1.1) #1.1 is gamma value in notes
		}
	}
	
	#center theta_j's
	theta.mat = scale(theta.mat, center=T, scale=F)
	f.hat.mat = theta.to.f(theta.mat, x)
	firstL = (L.seq[which(maximizer==1)])[1]
	
	return(list(f.hat.mat=f.hat.mat,theta.hat.mat=theta.mat,beta0.hat=beta0,y.hat=(beta0+rowSums(theta.mat)),objective=obj.vec,n.iter=n.iter,n.steps=sum(maximizer),x=x,L=firstL, family="gaussian", method="GGD.backtrack"))
}

#########################################################################################################
# logistic loss, generalized gradient descent
#########################################################################################################

flam.logistic.helper = function(y,lambda,alpha,x,initial.theta.mat=NULL,initial.f.mat=NULL,n, p, rank.x, order.x, tolerance) {
		
	#intialize matrix for storing theta hat estimates
	if (!is.null(initial.f.mat)) {
	  initial.theta.mat = f.to.theta(initial.f.mat, x, rank.x)
	} else if (is.null(initial.theta.mat)) {
	  initial.theta.mat = matrix(0, nrow=n, ncol=p)
	}
  
  beta0 = 0
	theta.mat = flamsteplogistic(initial.theta.mat, y, lambda, alpha, n, p, order.x, rank.x, beta0, tolerance)
	
	#intercept
	beta0 = beta0 + sum(colMeans(theta.mat))
	theta.mat = scale(theta.mat, center=T, scale=F)
	
	f.hat.mat = theta.to.f(theta.mat, x)
	
	return(list(f.hat.mat=f.hat.mat,theta.hat.mat=theta.mat,beta0.hat=beta0,y.hat=expit(beta0+rowSums(theta.mat)),x=x, family="binomial", method="GGD"))
}

#########################################################################################################
# misc. functions that are called by helper FLAM functions
#########################################################################################################

calc.sq.err = function(theta.mat, y) {
	if (ncol(theta.mat)>1) sq.err = 0.5 * sum((y - rowSums(theta.mat))^2) else sq.err = 0.5 * sum((y - theta.mat)^2)
	return(sq.err)
}

soft.scale.est = function(theta.vec, tuning.parameter) {
	
	if (sum(theta.vec^2)!=0) theta.vec = max(0, 1 - (tuning.parameter/sqrt(sum(theta.vec^2)))) * theta.vec
	return(theta.vec)
}

fit.est = function(col, resid.mat, rank.x, order.x, tuning.parameter) {
	
	resid.vec = resid.mat[,col][order.x[,col]]
	
  betasol = tf_dp(n=length(resid.vec), y=resid.vec, lam=tuning.parameter)
  est = betasol[rank.x[,col]]
	return(est)
}

logit=function(p) {log(p/(1-p))}

expit=function(x) {exp(x)/(1+exp(x))}


#########################################################################################################
# predict function that provides predictions for new set of predictors (new.x) 
# based on fit of model returned by flam function (object)
#########################################################################################################

helper.predict.flam = function(x, new.x, fhat) {
	
	x.sort = sort(x)
	if (!is.na(match(new.x, x.sort))) {
		theta.hat = fhat[match(new.x, x.sort)]
		} else if (min(x.sort) > new.x) {
			theta.hat = fhat[1]
		} else if (max(x.sort) < new.x) {
			theta.hat = fhat[length(x.sort)]
		} else {
			#linearly interpolate for x not in original fit
			y1 = fhat[max(which(x.sort < new.x))]; y2 = fhat[min(which(x.sort > new.x))]
			x1 = x.sort[max(which(x.sort < new.x))]; x2 = x.sort[min(which(x.sort > new.x))]
			theta.hat = y1 + (new.x - x1) * (y2-y1)/(x2-x1)
		}
			
	return(theta.hat)
}

predict.flam = function(object, new.x, lambda, alpha, ...) {
	
	if (is.null(nrow(new.x))) new.x = matrix(new.x, ncol=1)
	new.n = nrow(new.x); new.p = ncol(new.x)
	
	#checks
	if (new.p!=ncol(object$x)) stop("'new.x' must have the same number of columns as 'object$x'")
	if (length(lambda)!=1 | length(alpha)!=1) stop("Provide a single value for both 'lambda' and 'alpha'")
	if (alpha<0 | alpha>1) stop("'alpha' must be in [0,1]")
	if (lambda<=0) stop("'lambda' must be positive")

	theta.hat.mat = matrix(NA, nrow=new.n, ncol=new.p)
	
	#tuning parameter index
	index.lambda = which(rep(lambda,length(object$all.lambda))==object$all.lambda)
	index.alpha = which(rep(alpha,length(object$all.alpha))==object$all.alpha)
	if (length(index.alpha)>0 & length(index.lambda)>0) index = intersect(index.lambda, index.alpha) else index = logical(0)
	
	#lambda and alpha not in both object$all.lambda and object$all.alpha (index is empty)
	
	if (length(index)==0) {
		
		n = nrow(object$x); p = ncol(object$x)
		rank.x = apply(object$x,2,rank)
		order.x = apply(object$x,2,order)
		
		#figure out initial f.mat
		initial.f.mat = NULL
		
		if (alpha %in% object$all.alpha) { #alpha in alpha.seq
			index.lambda = which(object$all.lambda>lambda)
			index.alpha = match(alpha, object$all.alpha)
			index = intersect(index.lambda, index.alpha)
			if (length(index)!=0) initial.f.mat=object$f.hat.list[[max(index)]]
		} else if (lambda %in% object$all.lambda) { #lambda in lambda.seq
			index.lambda = match(lambda, object$all.lambda) 
			index.alpha = which(object$all.alpha>alpha)
			index = intersect(index.lambda, index.alpha)
			if (length(index)!=0) initial.f.mat=object$f.hat.list[[max(index)]]
		} else { #neither in previous sequences
			index.lambda = which(object$all.lambda>lambda)
			index.alpha = which(object$all.alpha>alpha)
			index = intersect(index.lambda, index.alpha)
			if (length(index)!=0) initial.f.mat=object$f.hat.list[[max(index)]]			
		}
		
		#fit model
		if (object$family=="gaussian" & object$method=="BCD") {		
			fit = flam.helper(y=object$y, lambda=lambda, alpha=alpha, x=object$x, initial.f.mat=initial.f.mat, n=n, p=p, order.x=order.x, rank.x=rank.x, tolerance=object$tolerance)
		}
		
		if (object$family=="gaussian" & object$method=="GGD") {
			fit = flam.ggd.helper(y=object$y, lambda=lambda, alpha=alpha, x=object$x, rank.x=rank.x, order.x=order.x, n=n, p=p, initial.f.mat=initial.f.mat, tolerance=object$tolerance)
		}
		
		if (object$family=="gaussian" & object$method=="GGD.backtrack") {
			fit = flam.ggd.helper_backtrack(y=object$y, lambda=lambda, alpha=alpha, x=object$x, n=n, p=p, initial.L=1, rank.x=rank.x, order.x=order.x, initial.f.mat=initial.f.mat, tolerance=object$tolerance)
		}
		
		if (object$family=="binomial") fit = flam.logistic.helper(y=object$y, lambda=lambda, alpha=alpha, x=object$x, n=n, p=p, rank.x=rank.x, order.x=order.x, initial.f.mat=initial.f.mat, tolerance=object$tolerance)

		for (i in 1:new.n) {
			theta.hat.mat[i,] = sapply(1:new.p, function(j) helper.predict.flam(x=fit$x[,j], fhat=fit$f.hat.mat[,j], new.x=new.x[i,j]))
		}
		
		if (object$family=="gaussian") y.hat.new = fit$beta0.hat + rowSums(theta.hat.mat) else if (object$family=="binomial") y.hat.new = expit(fit$beta0.hat + rowSums(theta.hat.mat))

	} else {	
		
	#lambda in object$all.lambda and alpha in object$all.alpha
			
		for (i in 1:new.n) {
			theta.hat.mat[i,] = sapply(1:new.p, function(j) helper.predict.flam(x=object$x[,j], fhat=object$f.hat.list[[index]][,j], new.x=new.x[i,j]))
		}
		
		if (object$family=="gaussian") y.hat.new = object$beta0.hat.vec[index] + rowSums(theta.hat.mat) else if (object$family=="binomial") y.hat.new = expit(object$beta0.hat.vec[index] + rowSums(theta.hat.mat))
	}
	
	return(y.hat.new)
}

#########################################################################################################
# flamCV that does K-fold cross-validation to choose the value for lambda and alpha
# returns fit for value of cross-validated lambda and alpha
#########################################################################################################

flamCV.helper = function(lambda.seq,alpha.seq,y,x,train,method,family,tolerance) {

	fit = flam(y=y[train], lambda.seq=lambda.seq, alpha.seq=alpha.seq, x=x[train,], family=family, method=method, tolerance=tolerance)
	
	n.tuning = length(fit$all.lambda)
	error.vec = rep(NA, n.tuning)

	for (i in 1:n.tuning) {

		y.test.hat = predict.flam(object=fit, new.x=x[-train,], lambda=fit$all.lambda[i], alpha=fit$all.alpha[i])
		error.vec[i] = mean((y[-train] - y.test.hat)^2)
		
	}
	
	return(error.vec)
}

flamCV = function(x, y, lambda.min.ratio = 0.01, n.lambda = 50, lambda.seq=NULL, alpha=1, family="gaussian", method="BCD", fold=NULL, n.fold=NULL, seed=NULL, within1SE=T, tolerance=10e-6) {
	
	call = match.call()
	fit = list()
	
	if (!is.null(seed)) set.seed(seed)
	if (is.null(nrow(x))) x = matrix(x,ncol=1)
	n = nrow(x); p = ncol(x)
	
	if (is.null(lambda.seq)) {
		max.lam = maxLambda(y=y, x=x, alpha=alpha)
		lambda.seq = exp(seq(log(max.lam),log(max.lam*lambda.min.ratio),len=n.lambda))
	}
	
	#checks
	if (length(y)!=n) stop("The length of 'y' must equal the number of rows of 'x'")
	if (length(lambda.seq)==1) stop("Provide a sequence of decreasing values for 'lambda.seq'")
	if (!(family %in% c("gaussian","binomial"))) stop("'family' not recognized - must be 'gaussian' or 'binomial'")
	if (min(lambda.seq)<=0) stop("Values in 'lambda.seq' must be positive")
	if (alpha<0 | alpha>1) stop("Value for 'alpha' must be in [0,1]")
	if (!(method %in% c("BCD","GGD","GGD.backtrack"))) stop("'method' not recognized - must be 'BCD', 'GGD', or 'GGD.backtrack'")
	if (!is.null(n.fold)) if ((n.fold<1 | n.fold>n)) stop("'n.fold' must be between 1 and the length of 'y'")
	if (!is.null(fold) & length(fold)!=n) stop("The length of 'fold' must equal the length of 'y'")
	
	if (is.null(fold) & !is.null(n.fold)) fold = sample(rep(1:n.fold, ceiling(n/n.fold))[1:n], n) else if (is.null(fold) & is.null(n.fold)) fold = sample(rep(1:10, ceiling(n/10))[1:n], n)
	n.fold = length(unique(fold))

	#make sure lambda.seq is decreasing
	lambda.seq = sort(lambda.seq, decreasing=TRUE)
	whole.alpha.seq = rep(alpha, each=length(lambda.seq))
	whole.lambda.seq = lambda.seq
	n.tuning = length(whole.lambda.seq)

	cv.error.mat = matrix(NA, nrow=n.fold, ncol=n.tuning)

	for (k in 1:n.fold) {
		print(paste("fold: ",k,sep=""))
		cv.error.mat[k,] = flamCV.helper(y=y, x=x, train=which(fold!=k), lambda.seq=lambda.seq, alpha.seq=alpha, method=method, family=family, tolerance=tolerance)
	}
	
	fit$mean.cv.error = apply(cv.error.mat,2,function(col,n.vec,n) sum(col*n.vec)/n, n.vec=table(fold), n=length(fold))
	index.cv = min(which(fit$mean.cv.error==min(fit$mean.cv.error)))
	fit$se.cv.error = apply(cv.error.mat,2,sd)/sqrt(nrow(cv.error.mat))
	if (within1SE==T) {
		error.cutoff = fit$mean.cv.error[index.cv] + fit$se.cv.error[index.cv]
		index.cv = min(which(fit$mean.cv.error <= error.cutoff))
	}
	fit$lambda.cv = whole.lambda.seq[index.cv]
	fit$alpha = alpha
	fit$index.cv = index.cv
	fit$flam.out = flam(y=y, lambda.seq=lambda.seq, alpha.seq=alpha, x=x, family=family, method=method)
	
	fit$fold = fold; fit$n.fold = n.fold; fit$within1SE = within1SE; fit$call = call
	class(fit) = "flamCV"
	return(fit)
}


#########################################################################################################
### FUNCTION THAT PLOTS FUNCTION FITS  
#########################################################################################################

plot.flam = function(x, index, n.plot=10, predictor.indicators=NULL, predictor.labels=NULL, outcome.label="outcome", ticks=F, col="dodgerblue", n.panel.width=NULL, n.panel.height=NULL, ...) {
	
	if (!is.null(predictor.indicators) & !is.null(predictor.labels)) use.only = T else use.only = F
	
	#checks
	if (!(index %in% 1:length(x$all.lambda))) stop("Provide a valid 'index'")
	if (!is.null(predictor.indicators) & !is.null(predictor.labels) & length(predictor.indicators)!=length(predictor.labels)) stop("Lengths of 'predictor.indicators' and 'predictor.labels' must match when both are specified")
	if (is.null(predictor.indicators) & !is.null(predictor.labels) & ncol(x$x)!=length(predictor.labels)) stop("Length of 'predictor.labels' must match the number of columns of 'x$x'")
	
	f.hat.mat = x$f.hat.list[[index]]
	if (is.null(predictor.indicators)) {
		f.norms = apply(f.hat.mat, 2, function(col) sum(col^2))
		non.sparse = (f.norms != 0)
		ranks = length(f.norms) + 1 - rank(f.norms, ties.method="first") #samples ranks for decreasing order
		predictor.indicators = which(ranks<=n.plot & non.sparse==1)
	}
	n.plot = length(predictor.indicators)
	
	if (x$family=="gaussian") y.title = paste("Change in mean ",outcome.label,sep="") else y.title = paste("Change in log odds of ",outcome.label,sep="")

	min.f = min(f.hat.mat); max.f = max(f.hat.mat); range = max.f - min.f
	if (!is.null(n.panel.width) & !is.null(n.panel.height)) if ((n.panel.width*n.panel.height)<n.plot) n.panel.width <- n.panel.height <- NULL
	if (is.null(n.panel.width) | is.null(n.panel.height)) {
		if (n.plot<4) {n.panel.width = n.plot; n.panel.height = 1} else if (n.plot>4 & n.plot<9) {
			n.panel.width = ceiling(n.plot/2); n.panel.height = 2} else if (n.plot>9 & n.plot<16) {
			n.panel.width = ceiling(n.plot/3); n.panel.height = 3} else 
			n.panel.width <- n.panel.height <- ceiling(sqrt(n.plot))
	}
	
	col = rep(col, n.plot)
	
	par(mfrow=c(n.panel.height, n.panel.width))
	for (j in 1:n.plot) {
		
		if (use.only==T) x.title = predictor.labels[j] else if (is.null(predictor.labels[j])) x.title = paste("Predictor ",predictor.indicators[j],sep="") else x.title = predictor.labels[predictor.indicators[j]]
		var = predictor.indicators[j]
		plot(sort(x$x[,var]),f.hat.mat[,var],xlab=x.title,ylab=y.title,type="n",ylim=c(min.f-0.05*range,max.f+0.05*range))	
		if (ticks==T) points(x$x[,var],rep(min.f-0.04*range,nrow(x$x)),pch="-",col=rgb(col2rgb("gray")[1],col2rgb("gray")[2],col2rgb("gray")[3],100,maxColorValue=256))
		abline(h=0)
		points(sort(x$x[,var]),f.hat.mat[,var],type="l",col=col[j],lwd=2)	

	}
	cat("In most cases, only a subset of predictor fits are plotted - type '?plot.flam' for details.\n")
}

plot.flamSparsity = function(x, ...) {
	
	n.alpha = ncol(x$sparsity.mat)
	col.vec = c("seagreen1","orange","orchid","dodgerblue")
	col.vec = rep(col.vec, ceiling(n.alpha/4))
	line.type = rep(1:ceiling(n.alpha/4), each=4)
	
	par(mfrow=c(1,1))
	plot(x=1,type="n",xlim=c(min(log(x$lambda.seq)),max(log(x$lambda.seq))),ylim=c(0,100),ylab="Feature Sparsity (%)",xlab="log(Lambda)")
	for (i in 1:n.alpha) points(log(x$lambda.seq),x$sparsity.mat[,i],type="l",col=col.vec[i],lty=line.type[i],lwd=2)
	
	legend.vec = c()
	for (i in 1:n.alpha) legend.vec = c(legend.vec, paste("alpha=",x$alpha.seq[i],sep=""))
	legend("bottomright",bty="n",legend.vec,lty=line.type[1:n.alpha],col=col.vec[1:n.alpha])
}

plot.flamCV = function(x, showSE=T, ...) {
	
	min.y = min(x$mean.cv.error - x$se.cv.error); max.y = max(x$mean.cv.error + x$se.cv.error)
	
	par(mfrow=c(1,1))
	plot(x=1,type="n",xlim=c(min(log(x$flam.out$all.lambda)),max(log(x$flam.out$all.lambda))),ylim=c(min.y,max.y),ylab="Cross-validation Error",xlab="log(Lambda)")

	points(log(x$flam.out$all.lambda), x$mean.cv.error, cex=1.5, col="red",pch=16)
	
	if (showSE==T) arrows(x0=log(x$flam.out$all.lambda),x1=log(x$flam.out$all.lambda),y0=x$mean.cv.error-x$se.cv.error,y1=x$mean.cv.error+x$se.cv.error,angle=90,col="red",length=0.05,code=3,lwd=1.6)
	
	abline(v=log(x$flam.out$all.lambda[x$index.cv]),lty=2)
}

#########################################################################################################
### END OF FUNCTION THAT PLOTS FUNCTION FITS  
#########################################################################################################

#########################################################################################################
### FUNCTIONS THAT SUMMARIZE FUNCTION FITS  
#########################################################################################################

#summarizes overall call to flam (is.null(index)) or a specific fit (!is.null(index))
summary.flam = function(object, index=NULL, ...) {
	
	cat("Call: \n")	
	print(object$call)
	args = match.call()
	if (class(object)!="flam") stop("Provide 'object' of class 'flam'")
	if (is.null(index)) {
		cat("\nFLAM was fit using the tuning parameters:\n\n")
		cat("lambda: "); cat(round(unique(object$all.lambda),3))
		cat("\n\nalpha: "); cat(unique(object$all.alpha))	
		sparsity.vec = sapply(object$f.hat.list, function(f.hat.mat) 100*round(mean(apply(f.hat.mat, 2, function(col) sum(col^2)==0)),2))
		sparsity.mat = matrix(sparsity.vec, nrow=length(unique(object$all.lambda)))
		colnames(sparsity.mat) = unique(object$all.alpha); rownames(sparsity.mat) = round(unique(object$all.lambda),3)
		cat("\n\nThe percent sparsity for each combo of tuning parameters:\n")
		cat("(rows and columns correspond to values of lambda and alpha, respectively)\n\n")
		print(sparsity.mat, max=length(unique(object$all.alpha))*50)
		cat(paste("\nUse 'plot(summary(",args$object,"))' to plot this sparsity pattern.\n",sep=""))
		
		sparsity = list()
		sparsity$sparsity.mat = sparsity.mat
		sparsity$lambda.seq = unique(object$all.lambda)
		sparsity$alpha.seq = unique(object$all.alpha)
		class(sparsity) = "flamSparsity"
		invisible(sparsity)
		
	} else {
		
		cat("\nThe requested FLAM fit corresponds to:\n")
		cat(paste("lambda = ",object$all.lambda[index],sep=""))
		cat(paste("\nalpha = ",object$all.alpha[index],sep=""))	
		sparsity = 100*round(mean(apply(object$f.hat.list[[index]], 2, function(col) sum(col^2)==0)),2)
		n.nonsparse = sum(apply(object$f.hat.list[[index]], 2, function(col) sum(col^2)==0)!=1)
		cat(paste("\n\nThe percent sparsity was ",sparsity,"%.\n",sep=""))
		cat(paste("Thus ",object$num.non.sparse[index]," predictors were estimated to have a relationship with the outcome.\n",sep=""))
		cat("\nThe predictors with non-sparse fits:\n")
		cat(object$non.sparse.list[[index]])

		args.t = as.character(args)
		cat(paste("\n\nUse 'plot(",args.t[[2]],",",args.t[[3]],")' to plot the predictor fits.\n",sep=""))
	}
}

summary.flamCV = function(object, ...) {

	cat("Call: \n")	
	print(object$call)
	args = match.call()
	
	cat("\nFLAM was fit using the tuning parameters:\n\n")
	cat("lambda: "); cat(round(unique(object$flam.out$all.lambda),3))
	cat("\n\nalpha: "); cat(unique(object$flam.out$all.alpha))	
	
	cat(paste("\n\nCross-validation with K=",object$n.fold," folds was used to choose lambda.",sep=""))
	if (object$within1SE==T) {
		cat("\nLambda was chosen to be the largest value with CV error within one standard error of the minimum CV error.\n")
		cat(paste("\nThe chosen lambda was ",round(object$lambda.cv,3),".",sep=""))
		cat(paste("\nThis corresponds to ",object$flam.out$num.non.sparse[object$index.cv]," predictors having non-sparse fits.",sep=""))
		cat("\nThe predictors with non-sparse fits:\n")
		cat(object$flam.out$non.sparse.list[[object$index.cv]])
	} else {
		cat("\nLambda was chosen to be the value with the minimum CV error.\n")
		cat(paste("\nThe chosen lambda was ",round(object$lambda.cv,3),".\n",sep=""))
		cat(paste("\nThis corresponds to ",object$flam.out$num.non.sparse[object$index.cv]," predictors having non-sparse fits.",sep=""))
		cat("\nThe predictors with non-sparse fits:\n")
		cat(object$flam.out$non.sparse.list[[object$index.cv]])	
	}
	cat(paste("\n\nUse 'plot(",args$object,")' to plot CV error curve.\n",sep=""))
	cat(paste("Use 'plot(",args$object,"$flam.out,",args$object,"$index.cv)' to plot the predictor fits.\n",sep=""))

}

#########################################################################################################
### END OF FUNCTIONS THAT SUMMARIZE FUNCTION FITS  
#########################################################################################################

#########################################################################################################
### FUNCTION TO CALCULATE DEGREES OF FREEDOM  
#########################################################################################################

flamDOF = function(object, index) {
	
	#checks
	if (class(object)=="flam") {
		if (!(index %in% 1:length(object$all.lambda))) stop("Provide a valid 'index'")
		alpha = object$all.alpha[index]; lambda = object$all.lambda[index]
		f.hat.mat = object$f.hat.list[[index]]
	} else stop("'object' must be of class 'flam'")

	p = ncol(object$x); n = nrow(object$x)
	U = diag(n); for (i in 1:n) for (j in 1:n) if(j>i) U[i,j] = 1
	perm.list = makePermList(object$x)
	Zt = t(perm.list[[1]]) %*% U[,1:(n-1)]
	beta.vec = (solve(U) %*% f.hat.mat[,1])[1:(n-1),]
	for (k in 2:p) {
		Zt = cbind(Zt,t(perm.list[[k]]) %*% U[,1:(n-1)]) 
		beta.vec = c(beta.vec, (solve(U) %*% f.hat.mat[,k])[1:(n-1),])
	}
	Zt = scale(Zt,scale=F)

	if (alpha==1) {
		
	X_A = Zt[,which(beta.vec!=0)]
	df = qr(X_A)$rank + 1
		
	} else {
		
	dim = length(which(beta.vec!=0))
	penalty = matrix(0,nrow=dim,ncol=dim); counter = 1
	x = scale(U[,1:(n-1)],scale=F)
	for (k in 1:p) {
		if (length(which(beta.vec[(k-1)*(n-1) + 1:(n-1)]!=0))>0) {
			b = matrix(beta.vec[(k-1)*(n-1) + 1:(n-1)],ncol=1)
			x2 = x[,which(b!=0)]
			penalty[counter:(counter+length(which(b!=0))-1),counter:(counter+length(which(b!=0))-1)] = t(x2)%*%x2 / as.numeric((t(b)%*%t(x)%*%x%*%b)^0.5) - t(x2)%*%x%*%b%*%t(b)%*%t(x)%*%x2/as.numeric((t(b)%*%t(x)%*%x%*%b)^1.5)
			counter = counter + length(which(b!=0))
	}}
	X_A = Zt[,which(beta.vec!=0)]
	if (length(which(beta.vec!=0))==1) X_A = matrix(X_A,ncol=1)

	if (!is.null(ncol(X_A)) & ncol(X_A)>0) {
		if (any(!is.finite(t(X_A) %*% X_A + lambda * (1-alpha) * penalty))) {
			counter = 1
			for (k in 1:p) {
				if (length(which(beta.vec[(k-1)*(n-1) + 1:(n-1)]!=0))>0) {
					b = matrix(beta.vec[(k-1)*(n-1) + 1:(n-1)],ncol=1)
					x2 = x[,which(b!=0)]
					penalty[counter:(counter+length(which(b!=0))-1),counter:(counter+length(which(b!=0))-1)] = t(x2)%*%x2 / as.numeric((t(b)%*%t(x)%*%x%*%b)^0.5) - t(x2)%*%x%*%b%*%t(b)%*%t(x)%*%x2/as.numeric((t(b)%*%t(x)%*%x%*%b)^1.5)
					counter = counter + length(which(b!=0))
				}
			}
		}
		df = sum(diag(X_A %*% MASS::ginv(t(X_A) %*% X_A + lambda * (1-alpha) * penalty) %*% t(X_A))) + 1
	} else df = 1}
	
	return(df)	
}

#########################################################################################################
### END OF FUNCTION TO CALCULATE DEGREES OF FREEDOM  
#########################################################################################################

#########################################################################################################
### FUNCTIONS TO CALCULATE MAX LAMBDA TO CONSIDER  
#########################################################################################################

#########################################################################################################
#function to calculate minimum lambda_1 such that all the estimated f_j's are constant (alpha=1)
#########################################################################################################

maxLambda_a1 = function(y, x) {
  
  n = length(y)
  
  max.lambda.vec = sapply(1:ncol(x), function(i, y, x, n) maxLambda_a1_C_single(y - mean(y), x[,i], n), y=y, x=x, n=n)
  
  return(max(max.lambda.vec))
  
}

#########################################################################################################
#function to calculate minimum lambda_2 such that all the estimated f_j's are constant (alpha=0)
#########################################################################################################

maxLambda_a0 = function(y) {
	
	max.lambda.vec = sqrt(sum((y-mean(y))^2))
	
	return(max(max.lambda.vec))
	
}

#########################################################################################################
#function to calculate minimum lambda such that all the estimated f_j's are constant (any alpha)
#########################################################################################################

maxLambda = function(x, y, alpha) {
	
	if (is.null(nrow(x))) x = matrix(x,ncol=1)
	
	#checks
	if (length(y)!=nrow(x)) stop("The length of 'y' must equal the number of rows of 'x'")
	if (length(alpha)!=1) stop("Provide a single value for 'alpha'")
	if (alpha<0 | alpha>1) stop("'alpha' must be in [0,1]")

	max.lambda = min(maxLambda_a1(y,x)/alpha, maxLambda_a0(y)/(1-alpha))
	
  return(max.lambda)
}

#########################################################################################################
### END OF FUNCTIONS TO CALCULATE MAX LAMBDA TO CONSIDER  
#########################################################################################################

#########################################################################################################
### FUNCTIONS TO SIMULATE DATA  
#########################################################################################################

gen.theta = function(X.mat, fcns, index) {
	
	if (is.na(index)) {
		option = fcns; x = X.mat
	} else {
		option = fcns[index]
		x = X.mat[,index]
	}

	#f_j zero
	if (option==0) theta_j = rep(0,length(x))

	#scenario 1: all piecewise constant f_j
	if (option==11) {theta_j = 3 * (x<(-1)) + 2 + -7 * (x>0.5); theta_j = (theta_j - 0.1)/4.3232149}
	if (option==12) {theta_j = 12 * (x<(-0.2)) + -5 + 7 * (x>1.1); theta_j = (theta_j - 2.48)/4.8999837}
	if (option==13) {theta_j = 3 + -6 * (x<(-1.7) | x>0.8); theta_j = theta_j/3.0000150}
	if (option==14) {theta_j = -5 + 6 * (x>(-.7)) + 3 * (x>1.6); theta_j = (theta_j +.62)/3.4577044}

	#scenario 2: smooth f_j - SPAM paper functions
	if (option==21) {theta_j = -sin(1.5*x); theta_j = theta_j/0.6614151}
	if (option==22) {theta_j = x^3 + 1.5*(x-0.5)^2; theta_j = (theta_j - 3.500063)/4.8930086}
	if (option==23) {theta_j = -pnorm(x, mean=0.5, sd=0.8); theta_j = (theta_j + .4003183)/0.3874514}
	if (option==24) {theta_j = sin(exp(-0.5*x)); theta_j = (theta_j - .6195458)/0.2985293}
	
	#scenario 3: 2 f_j from scenario 1 + 2 f_j from scenario 2
	if (option==31) {theta_j = 3 * (x<(-1)) + 2 + -7 * (x>0.5); theta_j = (theta_j - .1)/4.3232149}
	if (option==32) {theta_j = 12 * (x<(-0.2)) + -5 + 7 * (x>1.1); theta_j = (theta_j - 2.48)/4.8999837}
	if (option==33) {theta_j = x^3 + 1.5*(x-0.5)^2; theta_j = (theta_j - 3.500063)/4.8930086}
	if (option==34) {theta_j = -pnorm(x, mean=0.5, sd=0.8); theta_j = (theta_j + .4003183)/0.3874514}
	
	#scenario 4: 'large n' with functions with large constant functions
	if (option==41) {theta_j = -5 * (x<0) + (10*(x)^2 - 5) * (x>=0 & x<1) + (sin((x-1)*20) + 5) * (x>=1); theta_j = (theta_j + 1.324793)/4.562398}
	if (option==42) {theta_j = 2 + 3*(cos((x+.75)*2*pi) - 1) * (x>=-.75 & x<(-.25)) + (-4 + (cos((x+.75)*2*pi) - 1)) * (x>=-.25 & x<=0.25) + -4*(x>.25); theta_j = (theta_j + 0.59994)/2.083371}
	if (option==43) {theta_j = (x<(-1)) * 3.125 + (x>1) * -3.125 + (x>=-1 & x<(-0.5))*(-50*(x+1)^5 + 3.125) + (x>=-0.5 & x<.5)*-50*(x)^5 + (x>=0.5 & x<=1)*(-50*(x-1)^5 + -3.125); theta_j = (theta_j)/2.752587}
	if (option==44) {theta_j = (x<(-1)) * ((cos((x+1+pi)*10-9*pi)+1)*5) + (x>1.5)* (cos((x-1.5)*20)-1); theta_j = (theta_j - 1.244387)/3.03563}
	
	return(theta_j)
}

sim.data = function(n, scenario, zerof, noise=1, family="gaussian") {

	data = list()
	
	#checks
	if (!(scenario %in% 1:4)) stop("'scenario' must be 1, 2, 3, or 4")
	if (!(family %in% c("gaussian","binomial"))) stop("'family' not recognized - must be 'gaussian' or 'binomial'")

	functions = scenario*10 + 1:4
	p = length(functions) + zerof
	functions = c(functions, rep(0, zerof))

	data$x = matrix(runif(n=n*p,min=-2.5,max=2.5), nrow=n, ncol=p)
	theta.mat = sapply(1:p,gen.theta,X.mat=data$x,fcns=functions)
	if (family=="gaussian") data$y = rowSums(theta.mat) + rnorm(n, mean=0, sd=noise)
	if (family=="binomial") data$y = rbinom(n=n,size=1,prob=expit(rowSums(theta.mat)))
	
	data$theta = theta.mat
	
	print("See example in '?sim.data' for code to plot generating functions.")
	return(data)
}

#########################################################################################################
### END OF FUNCTION TO SIMULATE DATA  
#########################################################################################################
