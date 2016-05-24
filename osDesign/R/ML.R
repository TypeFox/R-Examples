ML <-
function(nn0, nn1, x, N, case, group, cohort, alpha, maxiter = 100)
{
# S function to solve the concentrated lagrangian equations developed by 
# Breslow and Holubkov (1997). This function implements a modified
# Newton-Rhapson  algorithm to solve the equations. Gradients as computed by
# the authors were used to implement the NR algorithm.
	iter <- 0
	converge <- FALSE
	rerror <- 100
	nntot0 <- sum(nn0)	# Phase I control Total
	nntot1 <- sum(nn1)	# Phase I case Total
	nn <- nn0 + nn1
	nntot <- nntot0 + nntot1
	nstrata <- length(nn0)
	nobs <- length(case)
	ncovs <- ncol(x)
	stpmax <- 1 * (ncovs + nstrata)
 ca<-case
 co<-N-case
 u <- sort(unique(group))
 strt.id <- outer(group,u,FUN="==")
 strt.id <- matrix(as.numeric(strt.id),nrow=length(group),ncol=length(u))
	ee1 <- cbind(matrix(1, nstrata, 1), matrix(0, nstrata, nstrata))
	ee2 <- cbind(matrix(0, nobs, 1), strt.id)
	ee <- rbind(ee1, ee2)
 n1 <- apply(strt.id * case, 2, sum)
 # n1= Phase II case sample sizes
	n0 <- apply(strt.id * (N - case), 2, sum)
	#  n0 = phase II control sample sizes   
	if((any(n0 == 0)) || (any(n1 == 0)))
		stop("Zero cell frequency at phase II")
	n <- n0 + n1
	n1a <- c(nntot1, n1)
	n0a <- c(nntot0, n0)
	grpa <- c(rep(1, nstrata), (group + 1))	# Augmented group indicator
	na <- n0a + n1a
	yy <- c(nn1, case)
	identity <- diag(rep(1, nstrata))
	x0 <- cbind(identity, matrix(0, nrow = nstrata, ncol = ncovs))
	x <- cbind( - strt.id, x)	# Augmented covariate matrix
	xx <- rbind(x0, x)
	ofs <- log(n1a/n0a)
	ofs[1] <- ofs[1] - log(alpha)
	if(cohort)
		ofs[1] <- ofs[1] - log(nntot1/nntot0) + log(alpha)
	ofs <- ofs[grpa]
	repp <- c(nn1 + nn0, N)	# Augmented binomial denominator
	m <- glm(cbind(yy, repp - yy) ~ -1 + xx + offset(ofs), family = binomial, x = TRUE, control = glm.control(epsilon = 1e-10, maxit = 100))
	gamm.schill <- as.vector(m$coef)
	gamm0 <- gamm <- gamm.schill	# Initialize by Schill's estimates
	errcode <- 0
	while((iter < maxiter) && (rerror > 1e-10))
	{
		#   cat("No of iterations=", iter, "\n")
		l <- lagrange(gamm0, nstrata, nn1, nn0, n1a, n0a, grpa, repp, xx, ofs, yy)
		h <- lagrad(gamm0, nstrata, nobs, ncovs, nn1, nn0, n1a, n0a, ee, repp, grpa, n0, n1, xx, ofs)
 	p <-  - solve(h) %*% l	# Newton's direction
		sp2 <- sqrt(sum(p^2))

		if(sp2 > stpmax)
			p <- p * (stpmax/sp2)	
		# Scale if the attempted stepis too big
		gamm <- lnsrch(gamm0, 0.5 * sum(l^2), t(h) %*% l, p, nstrata, nn1, nn0, n1a, n0a, grpa, repp, xx, ofs, yy)	
		# Search for new value along the line of Newton's direction
		rerror <- max(abs(gamm - gamm0)/pmax(abs(gamm0), 0.10000000000000001))
		gamm0 <- gamm  

   iter <- iter + 1
	}
	if(iter < maxiter) converge <- TRUE	#cat("No of iterations=", iter, "\n")
	h <- lagrad(gamm0, nstrata, nobs, ncovs, nn1, nn0, n1a, n0a, ee, repp, grpa, n0, n1, xx, ofs)
	fvv <- xx %*% gamm + ofs
	fvv <- as.vector(exp(fvv)/(1 + exp(fvv)))
	t1 <- nn1 - nn * fvv[1:nstrata]
	mu <- 1 - t1/n1
	t1 <- c(0, t1)
	mu <- c(1, mu)
	qn <- repp * n0a[grpa]
	qd1 <- n0a[grpa]
	qd2 <- (1 - mu[grpa]) * (n1a[grpa] - na[grpa] * fvv)
	q <- qn/(qd1 + qd2)
	# 
	# Verify constraints: equations (8) and (9) of Breslow and Holubkov
	#qq <- q/na[grpa]
	#h0 <- round(as.vector(t(qq)%*%ee), 10)
	#h1 <- as.vector(round(t(n1a - na*t(ee)%*%(fvv*qq)),digits=4))
	#print("Check that q's sum to 1 within strata")
	#print(h0)
	#print("Check that constraints for i=1 hold")
	#print(h1)
	##
	#fail <- FALSE
	#if(sum(h0 != 1) > 0) fail <- TRUE
	#if(sum(h1 != 0) > 0) fail <- TRUE
	#if(fail == TRUE) cat("\n WARNING: ML estimates don't satisfy appropriate constraints\n")
	##
	gg <- t(xx) %*% ((1 - fvv) * fvv * q * ee)
	t00 <- (1 - fvv)^2 * q * ee
	d00 <- rep.int(1, nrow(t00)) %*% t00
	t01 <- (1 - fvv) * fvv * q * ee
	d01 <- rep.int(1, nrow(t01)) %*% t01
	t11 <- fvv * fvv * q * ee
	d11 <- rep.int(1, nrow(t11)) %*% t11
	d00 <- diag(as.vector(d00))
	d01 <- diag(as.vector(d01))
	d11 <- diag(as.vector(d11))
	tinfo <- (1 - fvv) * fvv * q * xx
	info <- t(xx) %*% tinfo
	cove <- solve(info)
	gig <- t(gg) %*% cove %*% gg
	zz1 <- cbind((gig + d00), ( - gig + d01))
	zz2 <- cbind(( - gig + d01), (gig + d11))
	zz <- rbind(zz1, zz2)
	zz <- solve(zz)
	gg <- cbind(gg,  - gg)
	gig <- cove %*% gg %*% zz %*% t(gg) %*% cove
	cov1 <- cove - gig	
	# Model covariance matrix using Atchinson and Silvey formula 
	# Computation of empirical variance-covariance matrix using within group scores
	ta <- n0a * (n1a - t1)
	tb <- na * t1
	g <- (ta[grpa] * fvv)/(ta[grpa] + tb[grpa] * (1 - fvv))
	uj0 <-  - t(m$x) %*% (ee * (repp - yy) * g)
	uj1 <- t(m$x) %*% (ee * yy * (1 - g))
	gg <- t(m$x * yy * (1 - g)^2) %*% m$x + t(m$x * (repp - yy) * g^2) %*% m$x
	#u1j0 <- -t(xx)%*%(ee*(repp-yy)*g)        
	#u1j1 <- t(xx)%*%(ee*yy*(1-g))
	#print("uj0:")
	#print(uj0)
	#print("uj1:")
	#print(uj1)
	# gg <- t(xx*((1-g)^2*yy + g^2*(repp-yy)))%*%xx
	gig <- gg - t(t(uj0)/n0a) %*% t(uj0) - t(t(uj1)/n1a) %*% t(uj1)	
	#gig <- gg - t( t(u1j0) / n0a)  %*% t(u1j0) -  t( t(u1j1) / n1a ) %*% t(u1j1)
	#print("gig")
	#print(gig)
	h <- solve(h)
	cov2 <- t(h) %*% gig %*% h	
	# Adjusting term for asymptotic variance-covariance matrix for prospective first stage sampling
	adjust <- (1/nntot0) + (1/nntot1)
	adjust <- matrix(adjust, nrow = (nstrata + 1), ncol = (nstrata + 1))
	if(cohort)
	{
		cov1[1:(nstrata + 1), 1:(nstrata + 1)] <- cov1[1:(nstrata + 1), 1:(nstrata + 1)] + adjust
		cov2[1:(nstrata + 1), 1:(nstrata + 1)] <- cov2[1:(nstrata + 1), 1:(nstrata + 1)] + adjust
	}
	cov1 <- cov1[(nstrata + 1):(nstrata + ncovs), (nstrata + 1):(nstrata + ncovs)]	
	# Exctract the covariance matrix for the regression coefficients
	cov2 <- cov2[(nstrata + 1):(nstrata + ncovs), (nstrata + 1):(nstrata + ncovs)]	
	# Exctract the covariance matrix for the regression coefficients
	# jc <- jc[(nstrata + 1):(nstrata + ncovs), (nstrata + 1):(nstrata + ncovs)]	
	delta <- gamm[1:nstrata]
	b <- gamm[(nstrata + 1):(nstrata + ncovs)]
	q <- q/na[grpa]
	
	#validation of constraints, calculate h_ij from eq. (9) in Breslow and Holubkov
	xx2<-xx[(nstrata+1):nrow(xx),]
	pred<-xx2%*%gamm
 help<-rep(0,nstrata)
	for (i in 1:nstrata)
	{
	   k<-xx2[,i]
    help[i]<--sum(k)
	} 
	help<-as.vector(help)
	#help<-nrow(xx2)/nstrata
 n1_new<-rep(n1,help)
 n0_new<-rep(n0,help)
 p1<-n1_new*exp(pred)/(n0_new+n1_new*exp(pred))
 xx3<-xx[1:nstrata,]
 pred2<-xx3%*%gamm
 nntot0_new<-rep(nntot0,nstrata)
 nntot1_new<-rep(nntot1,nstrata)
 if(cohort)
 {
 P<-exp(pred2)/(1+exp(pred2))
 }
 else 
 { 
 P<-nntot1_new*exp(pred2)/(nntot0_new+nntot1_new*exp(pred2))
 }
 T<-nn1-(nn1+nn0)*P
 T_new<-rep(T,help)
 q<-N/(n0_new+n1_new)*n1_new*n0_new/(n0_new*n1_new+T_new*(n1_new-(n1_new+n0_new)*p1))
 p0<-1-p1
 h1<-rep(0,nstrata)
 h0<-rep(0,nstrata)
 start<-1
 for (i in 1:nstrata)
 {
   sum_1<-0
   sum_0<-0
   end<-start+help[i]-1   
     sum_1<-t(p1[start:end])%*%q[start:end]
     sum_0<-t(p0[start:end])%*%q[start:end]
   h1[i]<-n1[i]-(n1[i]+n0[i])*sum_1
   h0[i]<-n0[i]-(n1[i]+n0[i])*sum_0
   start<-start+help[i]
 }
 
#	out <- list(coef = b, covm = cov1, cove = cov2, delta = delta, mu = mu, h0=h0, h1=h1, k=k)
	fail <- FALSE
	if(sum(round(h0, 10) != 0) > 0) fail <- TRUE
	if(sum(round(h1, 10) != 0) > 0) fail <- TRUE
	if(fail == TRUE) cat("\n WARNING: ML estimates don't satisfy appropriate constraints\n\n")
	out <- list(coef = b, covm = cov1, cove = cov2, delta = delta, mu = mu, fail = fail)
	out
}
