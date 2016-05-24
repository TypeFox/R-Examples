WL <-
function(nn0, nn1, x, N, case, group, cohort, alpha)
{
	# S function for the weighted likelihood or Horvitz-Thompson type of estimator.        
	nntot0 <- sum(nn0)	# Phase I control Total
	nntot1 <- sum(nn1)	# Phase I case Total
  u <- sort(unique(group))
  strt.id <- outer(group,u,FUN="==")
  strt.id <- matrix(as.numeric(strt.id),nrow=length(group),ncol=length(u))
  n1 <- apply(strt.id * case, 2, sum)	      # n1= Phase II case sample sizes
	n0 <- apply(strt.id * (N - case), 2, sum) # n0 = phase II control sample sizes   
	if((any(n0 == 0)) || (any(n1 == 0)))
		stop("Zero cell frequency at phase II")
	strt.id <- rbind(strt.id, strt.id)
	N <- c(case, N - case)
	case <- c(rep(1, length(case)), rep(0, length(case)))
	group <- c(group, group)
	x <- rbind(x, x)
	nstrta <- length(nn0)
	ofs <- 0
	if(!cohort)
		ofs <- log(sum(nn1)/sum(nn0)) - log(alpha)
	ofs <- rep(ofs, length(case))
	pw0 <- nn0[group]/n0[group]
	pw1 <- nn1[group]/n1[group]
	pw <- case * pw1 + (1 - case) * pw0
	Nw <- N*pw
	lp <- rep(0, length(case))
	z <- rep(0, length(case))
	z[case == 1] <- Nw[case == 1]
	m <- glm(cbind(z, Nw - z) ~ -1 + x + offset(ofs), family = binomial, x=TRUE, 
	         control = glm.control(epsilon = 9.9999999999999995e-07, maxit = 20))
  #  weights=pw, start = lp, 
	uj0 <-  - t(m$x) %*% (strt.id * (1 - case) * m$fitted.values * N)	
	#   Within group sum of scores for controls
	uj1 <- t(m$x) %*% (strt.id * case * (1 - m$fitted.values) * N)	
	#     Within group sum of scores for cases
	g0 <- t(m$x * as.vector((1 - case) * m$fitted.values^2 * N * pw0^2)) %*% m$x           ## CHANGED 29TH OCT 2007
	g1 <- t(m$x * as.vector(case * (1 - m$fitted.values)^2 * N * pw1^2)) %*% m$x
	identity <- diag(rep(1, nstrta))
	if((any(n1 == 0)) || (any(n0 == 0)))
		stop("No cell should be empty in phase II")
  ##
  b0temp <- as.vector((nn0 * (nn0 - n0))/n0^3)                                           ## CHANGED 30TH OCT 2007
	b0 <- diag(b0temp, nrow=length(b0temp))
  ##
  b1temp <- as.vector((nn1 * (nn1 - n1))/n1^3)                                           ## CHANGED 30TH OCT 2007
	b1 <- diag(b1temp, nrow=length(b1temp))
  ##
  bb0temp <- as.vector((n0/(nntot0 * nn0)))                                              ## CHANGED 30TH OCT 2007
	bb0 <- diag(bb0temp, nrow=length(bb0temp))
  ##
  bb1temp <- as.vector((n1/(nntot1 * nn1)))                                              ## CHANGED 30TH OCT 2007
	bb1 <- diag(bb1temp, nrow=length(bb1temp))
	##
  cov <- summary(m)$cov.unscaled
	cove <- g1 + g0 - uj0 %*% b0 %*% t(uj0) - uj1 %*% b1 %*% t(uj1)
	if(!cohort)
		cove <- cove - uj0 %*% bb0 %*% t(uj0) - uj1 %*% bb1 %*% t(uj1)
	cove <- cov %*% cove %*% cov
	z <- list(coef = m$coef, cove = cove, fail = FALSE)
	z
}
