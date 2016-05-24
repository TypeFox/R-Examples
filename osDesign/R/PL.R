PL <-
function(nn0, nn1, x, N, case, group, cohort, alpha)
{
# S-function for "Pseudo Likelihood"  method  as developed by Breslow and Cain (1988)
	nntot0 <- sum(nn0)	# Phase I control Total
	nntot1 <- sum(nn1)	# Phase I case Total
	nstrata <- length(nn0)
  u <- sort(unique(group))
  strt.id <- outer(group,u,FUN="==")
  strt.id <- matrix(as.numeric(strt.id),nrow=length(group),ncol=length(u))
  n1 <- apply(strt.id * case, 2, sum)       # n1 = Phase II case sample sizes
	n0 <- apply(strt.id * (N - case), 2, sum) # n0 = phase II control sample sizes   
	if((any(n0 == 0)) || (any(n1 == 0)))
		stop("Zero cell frequency  at phase II")
	ofs <- log(n1/n0) - log(nn1/nn0)
	if(!cohort)
		ofs <- ofs + log(nntot1/nntot0) - log(alpha)
	ofs <- ofs[group]
	lp <-  - ofs
	pw <- rep(1, length(case))	
	#  Fitting model using standard GLM procedure
	m <- glm(cbind(case, N - case) ~ -1 + x + offset(ofs), family = binomial, weights = pw, x = TRUE, control = glm.control(epsilon = 9.9999999999999995e-07, maxit = 20))
	fv <- m$fitted.values
	nhat0 <- apply(strt.id * (1 - fv) * N, 2, sum)	
	#  Within group sum of fitted values for controls
	nhat1 <- apply(strt.id * fv * N, 2, sum)	
	# Within group sum of fited values
	# Adjusted covariance Matrix
	uj0 <-  - t(m$x) %*% (strt.id * (N - case) * m$fitted.values)	
	#   Within group sum of scores for controls
	uj1 <-  - t(m$x) %*% (strt.id * case * (1 - m$fitted.values))	
	#     Within group sum of scores for cases
	g0 <- t(m$x * (N - case) * fv^2) %*% m$x
	g1 <- t(m$x * case * (1 - fv)^2) %*% m$x
	a <- t(m$x) %*% (strt.id * (1 - fv) * fv * N)
	zz <- matrix(1, nrow = nstrata, ncol = nstrata)
	identity <- diag(rep(1, nstrata))
	b0 <- identity/as.vector(n0)
	b1 <- identity/as.vector(n1)
	bb0 <- (identity/as.vector(nn0))
	bb1 <- (identity/as.vector(nn1))
	if(!cohort)
	{
		bb0 <- bb0 - (zz/nntot0)
		bb1 <- bb1 - (zz/nntot1)
	}
	aba <- a %*% bb0 %*% t(a) + a %*% bb1 %*% t(a)
	info <- t(m$x) %*% diag(m$weight) %*% m$x
	info2 <- solve(summary(m)$cov.unscaled)
	ghat <- info - a %*% b0 %*% t(a) - a %*% b1 %*% t(a)
	ghate <- g1 + g0 - uj0 %*% b0 %*% t(uj0) - uj1 %*% b1 %*% t(uj1)
	cov <- summary(m)$cov.unscaled
	cove <- cov %*% (ghate + aba) %*% cov	
	# Empirical variance-covariance matrix
	cove <- cov %*% (ghate + aba) %*% cov	
	# Model based variance-covariance matrix
	covm <- cov %*% (ghat + aba) %*% cov
	return(list(coef = m$coef, covm = covm, cove = cove, fail = FALSE))
}
