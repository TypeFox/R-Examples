####################################################################################################
# CREATE glmdmEST FUNCTION AS A WRAPPER (08/03/2011)
####################################################################################################
glmdmEST <- function(Y, X, family=gaussian, num.reps=1000, a1=3, b1=2, d=0.25, MM=15, VV=30,...) {

#########################################################################
# FUNCTIONS 
#########################################################################
# rdirichlet from LearnBayes (04/21/2011)
rdirichlet <- function (n, alpha) 
{
    k = length(alpha)
    z = array(0, dim = c(n, k))
    s = array(0, dim = c(n, 1))
    for (i in 1:k) {
        z[, i] = rgamma(n, shape = alpha[i])
        s = s + z[, i]
    }
    for (i in 1:k) {
        z[, i] = z[, i]/s
    }
    return(z)
}

# rdiscrete from e1071 (04/21/2011)
rdiscrete <- function (n, probs, values = 1:length(probs), ...) 
{
    sample(values, size = n, replace = TRUE, prob = probs)
}

# rinvgamma from MCMCpack (04/21/2011)
rinvgamma <- function (n, shape, scale = 1) 
{
    return(1/rgamma(n = n, shape = shape, rate = scale))
}

# rmvnorm from mvtnorm (04/21/2011)
rmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
    method = c("eigen", "svd", "chol")) 
{
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    sigma1 <- sigma
    dimnames(sigma1) <- NULL
    if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
        warning("sigma is numerically not symmetric")
    }
    method <- match.arg(method)
    if (method == "eigen") {
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% 
            t(ev$vectors)
    }
    else if (method == "svd") {
        sigsvd <- svd(sigma)
        if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }
    else if (method == "chol") {
        retval <- chol(sigma, pivot = TRUE)
        o <- order(attr(retval, "pivot"))
        retval <- retval[, o]
    }
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)
    retval
}

# rtnorm from msm (04/21/2011); FOR PROBIT MODEL
rtnorm <- function (n, mean = 0, sd = 1, lower = -Inf, upper = Inf) 
{
    if (length(n) > 1) 
        n <- length(n)
    mean <- rep(mean, length = n)
    sd <- rep(sd, length = n)
    lower <- rep(lower, length = n)
    upper <- rep(upper, length = n)
    lower <- (lower - mean)/sd
    upper <- (upper - mean)/sd
    ind <- seq(length = n)
    ret <- numeric(n)
    alg <- ifelse(lower > upper, -1, ifelse(((lower < 0 & upper == 
        Inf) | (lower == -Inf & upper > 0) | (is.finite(lower) & 
        is.finite(upper) & (lower < 0) & (upper > 0) & (upper - 
        lower > sqrt(2 * pi)))), 0, ifelse((lower >= 0 & (upper > 
        lower + 2 * sqrt(exp(1))/(lower + sqrt(lower^2 + 4)) * 
            exp((lower * 2 - lower * sqrt(lower^2 + 4))/4))), 
        1, ifelse(upper <= 0 & (-lower > -upper + 2 * sqrt(exp(1))/(-upper + 
            sqrt(upper^2 + 4)) * exp((upper * 2 - -upper * sqrt(upper^2 + 
            4))/4)), 2, 3))))
    ind.nan <- ind[alg == -1]
    ind.no <- ind[alg == 0]
    ind.expl <- ind[alg == 1]
    ind.expu <- ind[alg == 2]
    ind.u <- ind[alg == 3]
    ret[ind.nan] <- NaN
    while (length(ind.no) > 0) {
        y <- rnorm(length(ind.no))
        done <- which(y >= lower[ind.no] & y <= upper[ind.no])
        ret[ind.no[done]] <- y[done]
        ind.no <- setdiff(ind.no, ind.no[done])
    }
    stopifnot(length(ind.no) == 0)
    while (length(ind.expl) > 0) {
        a <- (lower[ind.expl] + sqrt(lower[ind.expl]^2 + 4))/2
        z <- rexp(length(ind.expl), a) + lower[ind.expl]
        u <- runif(length(ind.expl))
        done <- which((u <= exp(-(z - a)^2/2)) & (z <= upper[ind.expl]))
        ret[ind.expl[done]] <- z[done]
        ind.expl <- setdiff(ind.expl, ind.expl[done])
    }
    stopifnot(length(ind.expl) == 0)
    while (length(ind.expu) > 0) {
        a <- (-upper[ind.expu] + sqrt(upper[ind.expu]^2 + 4))/2
        z <- rexp(length(ind.expu), a) - upper[ind.expu]
        u <- runif(length(ind.expu))
        done <- which((u <= exp(-(z - a)^2/2)) & (z <= -lower[ind.expu]))
        ret[ind.expu[done]] <- -z[done]
        ind.expu <- setdiff(ind.expu, ind.expu[done])
    }
    stopifnot(length(ind.expu) == 0)
    while (length(ind.u) > 0) {
        z <- runif(length(ind.u), lower[ind.u], upper[ind.u])
        rho <- ifelse(lower[ind.u] > 0, exp((lower[ind.u]^2 - 
            z^2)/2), ifelse(upper[ind.u] < 0, exp((upper[ind.u]^2 - 
            z^2)/2), exp(-z^2/2)))
        u <- runif(length(ind.u))
        done <- which(u <= rho)
        ret[ind.u[done]] <- z[done]
        ind.u <- setdiff(ind.u, ind.u[done])
    }
    stopifnot(length(ind.u) == 0)
    ret * sd + mean
}


###############################
# Basic info from data
###############################
n <- nrow(X)
# For the prior of m, we use Gamma(a,b)
# Here ab=MM and ab^2=VV. 
Sh <- (MM^2)/VV
Sc <- VV/MM

#########################################################################
# RECOGNIZE FAMILY ARGUMENT from glm (08/03/2011)
#########################################################################
	call <- match.call()
	if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }

####################################################################################################
# RUN A GLM TO GET hat.beta and hat.var
####################################################################################################

suppressWarnings(glm.out <- glm(Y ~ X, family=family))  # GLM BASED ON THE SPECIFIED FAMILY. (08/03/2011)
hat.beta <- coefficients(summary(glm.out))[,1]
hat.var  <- (coefficients(summary(glm.out))[,2])^2

####################################################################################################
# FIX PRIOR DISTRIBUTIONS AND PRIOR PARAMETERS
####################################################################################################
X  <- cbind(1,X) #WITH DATA, INTERCEPT TERM NEEDS TO BE ADDED

# PRE-LOOP EFFICIENCIES
tau.sq        <- rinvgamma(1,shape=a1,scale=b1)
# If tau.sq is infinity; 1.7e+100 is reasonable big?!
if(is.finite(tau.sq)==FALSE) tau.sq <- 1.7e+100
tau           <- sqrt(tau.sq)
sigma.sq      <- 1
sigma         <- sqrt(sigma.sq)
sigma.beta.sq <- sigma.sq * d
p             <- length(hat.beta)     
aa            <- chol(solve(t(X) %*% X + diag(1/sigma.beta.sq, nrow=ncol(X))))
SIG           <- t(X)%*%X + diag(1/sigma.beta.sq, nrow=ncol(X))				# VAR/COVAR MATRIX
SIG.inv       <- solve(SIG)								# INVERT HERE FOR EFFICIENCY

Z             <- rep(NA,length=n)		      # Z TO FILL-IN DURING LOOPING (FOR PROBIT MODEL)

beta          <- rmvnorm(1, mean=hat.beta, sigma=diag(hat.var))				# STARTING VALUES FROM ABOVE
beta          <- rbind( beta,matrix(rep(NA,num.reps*ncol(beta)),ncol=ncol(beta)) )	# MATRIX TO FILL IN DURING LOOP
q             <- c(rdirichlet(1, rep(1, n)))					# PROBABILITY OF ASSIGNMENT
A.n 	      <- sample(1:n,n,replace=TRUE,q)				 		# CREATES THE ASSIGNMENTS.
A.K 	      <- table(A.n)						 		# CREATES THE LIST OF OCCUPIED
A.K.labels    <- as.numeric(names(A.K))  						# LOCATIONS OF OCCUPIED
K	          <- length(A.K)								# NUMBER OF OCCUPIED TABLES
B             <- rep(0,length=n)							# LENGTH n FOR ALL CASES
B[A.K.labels] <- 1									# 1 WHERE OCCUPIED >0
B             <- B*rnorm(n,0,tau)							# ASSIGN psi VALUES TO OCCUPIED
nk            <- rep(0,n); for (i in 1:n) nk[A.n[i]] <- nk[A.n[i]] + 1			# COUNTS OCCUPANTS AT TABLES
psi           <- B[A.n] 								# MAP TABLE psi'S TO CASES
like.K        <- 0									# LOG-LIKE STARTS AT ZERO
X.beta1       <- X%*%beta[1,]								# MULTIPLY X AND beta IN ADVANCE


# SETUP TO SAVE GIBBS SAMPLER VALUES (THIS IS USED IN glmdm.linear; THE ONE USED IN glmdm.probit IS A LITTLE DIFFERENT (08/16/2011))
tau.sq.post <- xi.post <- like.K.save <- m.save <- NULL 
psi.post <- matrix(rep(NA,num.reps*n),ncol=n) 
K.save <- rep(NA,length=num.reps)
####################################################################################################
# NEW POSTERIOR FOR m FUNCTION
####################################################################################################
initial.m <- m <-10 ; #initial value of m
like      <- function(m.v){gamma(m.v)/gamma(m.v+n) * dgamma(m.v, shape=Sh, scale=Sc) * m.v^K}
loglike   <- function(m.v){lgamma(m.v)-lgamma(m.v+n) + dgamma(m.v, shape=Sh, scale=Sc, log=TRUE) + K*log(m.v)}
loglike.s <- function(m.v){lgamma(m.v)-lgamma(m.v+n) + dgamma(m.v, shape=Sh, scale=Sc, log=TRUE) + K*log(m.v) + nu*log(m.v)}



#########################################################################
#########################################################################
# FROM THE lOG-LIKELIHOOD OF K TO THE END, THE CALCULATION DEPENDS ON DIFFERENT LINK FUNCTIONS (08/03/2011)
#########################################################################
#########################################################################

## BEGINING: IF GAUSSIAN WITH IDENTITY LINK FUNCTION IS TRUE
if (eval(family)$family=="gaussian" & eval(family)$link=="identity") {

	for (i in 1:n) { 
        like.K <- like.K + dnorm(Y[i], mean=X.beta1[i] + psi[i], log=TRUE) + dnorm(psi[i], mean=0, sd=tau, log=TRUE)
        like.K <- exp(like.K + sum(lgamma(A.K)))	
    }    

####################################################################################################
# CALL GIBBS FUNCTION
####################################################################################################
# MCMC LOOPING
for (M in 1:num.reps)  {
    #   if (M %% 100 == 0) print(paste("finished iteration",M))
        
    #print("M-H ON THE 'A' MATRIX")
    p.A.old  <- (m^K)
    f.y.old  <- like.K
    mult.old <- dmultinom(x=A.K, prob=q[A.K.labels])

    #print("CREATE NEW 'A' MATRIX, 'can' STANDS FOR CANDIDATE") 
    pq                <- nk +1   
    new.q             <- rdirichlet(1, pq)
    A.n.can           <- sample(1:n,n,replace=TRUE,new.q)                                           
    A.K.can           <- table(A.n.can)                                                           
    A.K.labels.can    <- as.numeric(names(A.K.can))                                             
    K.can             <- length(A.K.can)                                                       
    B                 <- rep(0,length=n)                                                  
    B[A.K.labels.can] <- 1                                                               
    B                 <- B*rnorm(n,0,tau)                                               
    nk.can            <- rep(0,n); for (i in 1:n) nk.can[A.n.can[i]] <- nk.can[A.n.can[i]] + 1
    psi.can           <- B[A.n.can]                                                        
    p.A.can           <- (m^K.can)  
    like.K.can        <- 0                                                                      
    X.betaM           <- X%*%beta[M,]                                                          
    for (i in 1:n)
	    like.K.can    <- like.K.can + dnorm(Y[i], mean=X.betaM[i] + psi.can[i], log=TRUE) + dnorm(psi.can[i], mean=0, sd=tau, log=TRUE)
    like.K.can        <- exp(like.K.can + sum(lgamma(A.K.can)))
    f.y.can           <- like.K.can
    mult.can          <- dmultinom(x=A.K.can, prob=new.q[A.K.labels.can])

    #print("UPDATE 'A' and 'K', CREATE LIKELIHOOD")
    p.ratio <- p.A.can/p.A.old; f.ratio <- f.y.can/f.y.old; mult.ratio <- mult.can/mult.old
    if (is.finite(p.ratio) == FALSE) p.ratio       <- 1
    if (is.finite(f.ratio) == FALSE) f.ratio       <- 1
    if (is.finite(mult.ratio) == FALSE) mult.ratio <- 1

    # NOW UPDATE GIVEN WE ARE DONE WITH THE METROPOLIS STEP (ACCEPTING 'can' OR NOT)

    #print("UPDATE rho")
    rho <- p.ratio * f.ratio * mult.ratio 
    if (rho>runif(1))   {
	    A.n        <- A.n.can 		# A.n <- 1:n for regular random effects model
   	    A.K        <- A.K.can
    	A.K.labels <- A.K.labels.can    
        K          <- K.can
        K.save[M]  <- K
        nk         <- nk.can
        psi        <- psi.can 
        like.K     <- like.K.can
    }

    #print("UPDATE 'beta', 'tau' AND 'eta'")
    mn            <- SIG.inv %*% (t(X) %*% (Y-psi)) 
    Mb            <- t(aa) %*% array(rnorm(p), c(p, 1)) + mn 
    beta[M+1, ]   <- t(Mb) 
    
    A.m1          <- matrix(0, nrow=n, ncol=n)
    for (i in 1:n){ A.m1[i,A.n[i]] <- 1}
    A.K.m         <- A.m1[,which(apply(A.m1,2,sum)>0)]
    S.eta.inv     <- 1/sigma.sq *t(A.K.m) %*% A.K.m + diag(1/tau.sq, nrow=K)
    S.eta         <- solve(S.eta.inv)
    bb            <- chol(S.eta) 
    meta          <- 1/sigma.sq * S.eta %*% t(A.K.m) %*% (Y-X%*%beta[M+1,])
    eta           <- t(bb) %*% array(rnorm(K), c(K, 1)) + meta
    psi           <- A.K.m %*% eta
    psi.post[M, ] <- t(psi)
    
    tau.sq        <- rinvgamma(1,shape=(K)/2+a1,scale=sum(eta^2)/2+b1)
    # If tau.sq is infinity; 1.7e+100 is reasonable big?!
    if(is.finite(tau.sq)==FALSE) tau.sq <- 1.7e+100
    tau           <- sqrt(tau.sq)
    tau.sq.post   <- c(tau.sq.post,tau.sq)

    #print("UPDATE 'm'") 
    m.hat.s <- m.hess.s <- L.m.s.hat <- NULL
    mle.m       <- optim(par=initial.m, fn=loglike, method ="L-BFGS-B", lower = 0.01, upper = Inf, 
		          hessian = TRUE, control=list(fnscale=-1))
	m.hat       <- mle.m$par
	m.hessian   <- mle.m$hessian 
    L.m.hat     <- loglike(m.hat)
    for (nu in 1:2){ 
 	mle.m.s     <- optim(par=initial.m, fn=loglike.s, method ="L-BFGS-B", lower = 0.01, upper = Inf, 
		          hessian = TRUE, control=list(fnscale=-1))
	m.hat.s     <- c(m.hat.s, mle.m.s$par)
	m.hess.s    <- c(m.hess.s, mle.m.s$hessian) 
	Lms.hat     <- loglike.s(mle.m.s$par)
	L.m.s.hat   <- c(L.m.s.hat, Lms.hat)
    }
        
    mean.m <- sqrt(m.hessian/m.hess.s[1])*exp(L.m.s.hat[1] - L.m.hat)
    var.m  <- sqrt(m.hessian/m.hess.s[2])*exp(L.m.s.hat[2] - L.m.hat) - mean.m^2
    Sha    <- (mean.m^2)/var.m
    Sca    <- var.m/mean.m  
    cand   <- rgamma(1, shape=Sha, scale=Sca) 
    m      <- cand*rho + m*(1-rho)
    m.save <- c(m.save,m) 
   if (M %% 100 == 0) print(paste("finished iteration",M))
}
   beta.post <- apply(beta[(round(num.reps/2)+1):(num.reps+1),],2,mean)
   beta.sd   <- apply(beta[(round(num.reps/2)+1):(num.reps+1),],2,sd)
   psi.hat   <- apply(psi.post[(round(num.reps/2)+1):(num.reps),],2,mean)
   Y.hat     <- c(X %*% beta.post) + c(psi.hat)
   resid     <- c(Y - Y.hat)
   
                       			
	} ## END: IF GAUSSIAN WITH IDENTITY LINK FUNCTION IS TRUE (BEGIN FROM LINE 234)


## BEGINING: IF GAUSSIAN WITH IDENTITY LINK FUNCTION IS NOT TRUE
else {
	## BEGINING: IF BINOMIAL WITH PROBIT LINK FUNCTION IS TRUE
	if (eval(family)$family=="binomial" & eval(family)$link=="probit") {
		
		for (i in 1:n) { 
    		like.K <- like.K + Y[i] * pnorm( X.beta1[i] + psi[i], log.p=TRUE) + (1-Y[i]) * (1 - pnorm( X.beta1[i] + psi[i] )) + dnorm(psi[i], mean=0, sd=tau, log=TRUE)
			like.K <- exp(like.K + sum(lgamma(A.K)))			
	    }
	    
####################################################################################################
# CALL GIBBS FUNCTION
####################################################################################################

# MCMC LOOPING
for (M in 1:num.reps)  {
          
    #print("M-H ON THE 'A' MATRIX")
    p.A.old  <- (m^K)
    f.y.old  <- like.K
    mult.old <- dmultinom(x=A.K, prob=q[A.K.labels])
 	  
    #print("CREATE NEW 'A' MATRIX, 'can' STANDS FOR CANDIDATE") 
    pq                <- nk +1   
    new.q             <- rdirichlet(1, pq)
    A.n.can           <- sample(1:n,n,replace=TRUE,new.q)                                           
    A.K.can           <- table(A.n.can)                                                           
    A.K.labels.can    <- as.numeric(names(A.K.can))                                             
    K.can             <- length(A.K.can)                                                       
    B                 <- rep(0,length=n)                                                  
    B[A.K.labels.can] <- 1                                                               
    B                 <- B*rnorm(n,0,tau)                                               
    nk.can            <- rep(0,n); for (i in 1:n) nk.can[A.n.can[i]] <- nk.can[A.n.can[i]] + 1
    psi.can           <- B[A.n.can]                                                        
    like.K.can        <- 0                                                                      
    X.betaM           <- X%*%beta[M,]                                                          
    for (i in 1:n)
        like.K.can    <- like.K.can + Y[i] * pnorm( X.betaM[i] + psi.can[i], log.p=TRUE) +
                         (1-Y[i]) * (1 - pnorm( X.betaM[i] + psi.can[i] )) +
                         dnorm(psi.can[i], mean=0, sd=tau, log=TRUE)
    like.K.can        <- exp(like.K.can + sum(lgamma(A.K.can)))

    p.A.can           <- (m^K.can)  
    f.y.can           <- like.K.can
    mult.can          <-dmultinom(x=A.K.can, prob=new.q[A.K.labels.can])

    #print("UPDATE 'A' and 'K', CREATE LIKELIHOOD")
    p.ratio <- p.A.can/p.A.old; f.ratio <- f.y.can/f.y.old; mult.ratio <- mult.can/mult.old
    if (is.finite(p.ratio) == FALSE) p.ratio       <- 1
    if (is.finite(f.ratio) == FALSE) f.ratio       <- 1
    if (is.finite(mult.ratio) == FALSE) mult.ratio <- 1

    # NOW UPDATE GIVEN WE ARE DONE WITH THE METROPOLIS STEP (ACCEPTING 'can' OR NOT)

    #print("UPDATE rho")
    rho <- p.ratio * f.ratio * mult.ratio 
    if (rho>runif(1))   {
    	A.n        <- A.n.can  # A.n <- 1:n for regular random effects model
   	    A.K        <- A.K.can
        A.K.labels <- A.K.labels.can    
        K          <- K.can
        nk         <- nk.can
        psi        <- psi.can 
        like.K     <- like.K.can
    }

    #print("UPDATE 'z': Truncated Normal")
    for (j in 1:n){
	mean = X[j,]%*%beta[M,] + psi[j]
        Z[j] <- rtnorm(1, mean=mean, sd=1, lower=-Inf, upper=0)
        if (Y[j]==1) Z[j] <- rtnorm(1, mean=mean, sd=1, lower=0, upper=Inf)
    }

    #print("UPDATE 'beta', 'tau' AND 'eta'")
    mn            <- SIG.inv %*% (t(X) %*% (Z-psi)) 
    Mb            <- t(aa) %*% array(rnorm(p), c(p, 1)) + mn 
    beta[M+1, ]   <- t(Mb) 
    
    A.m1          <- matrix(0, nrow=n, ncol=n)
    for (i in 1:n){ A.m1[i,A.n[i]] <- 1}
    A.K.m         <- A.m1[,which(apply(A.m1,2,sum)>0)]
    S.eta.inv     <- 1/sigma.sq *t(A.K.m) %*% A.K.m + diag(1/tau.sq, nrow=K)
    S.eta         <- solve(S.eta.inv)
    bb            <- chol(S.eta) 
    meta          <- 1/sigma.sq * S.eta %*% t(A.K.m) %*% (Z-X%*%beta[M+1,])
    eta           <- t(bb) %*% array(rnorm(K), c(K, 1)) + meta
    psi           <- A.K.m %*% eta
    psi.post[M, ] <- t(psi)

    tau.sq        <- rinvgamma(1,shape=(K)/2+a1,scale=sum(eta^2)/2+b1)
    # tau.sq can be "inf" and needs to be bound: large number: 1.7e+100
    if(is.finite(tau.sq)==FALSE) tau.sq <- 1.7e+100
    tau           <- sqrt(tau.sq)
    tau.sq.post   <- c(tau.sq.post,tau.sq)

    #print("UPDATE 'm'") 
    m.hat.s <- m.hess.s <- L.m.s.hat <- NULL
    mle.m       <- optim(par=initial.m, fn=loglike, method ="L-BFGS-B", lower = 0.01, upper = Inf, 
		          hessian = TRUE, control=list(fnscale=-1))
	m.hat       <- mle.m$par
	m.hessian   <- mle.m$hessian 
    L.m.hat     <- loglike(m.hat)
    for (nu in 1:2){ 
 	mle.m.s     <- optim(par=initial.m, fn=loglike.s, method ="L-BFGS-B", lower = 0.01, upper = Inf, 
		          hessian = TRUE, control=list(fnscale=-1))
	m.hat.s     <- c(m.hat.s, mle.m.s$par)
	m.hess.s    <- c(m.hess.s, mle.m.s$hessian) 
	Lms.hat     <- loglike.s(mle.m.s$par)
	L.m.s.hat   <- c(L.m.s.hat, Lms.hat)
   }
        
    mean.m <- sqrt(m.hessian/m.hess.s[1])*exp(L.m.s.hat[1] - L.m.hat)
    var.m  <- sqrt(m.hessian/m.hess.s[2])*exp(L.m.s.hat[2] - L.m.hat) - mean.m^2
    Sha    <- (mean.m^2)/var.m
    Sca    <- var.m/mean.m  
    cand   <- rgamma(1, shape=Sha, scale=Sca) 
    ##LL##		Using's Dominik work-around
    test   <- (loglike(cand) - loglike(m) + dgamma(m, shape=Sha, scale=Sca, log=TRUE) - dgamma(cand, shape=Sha, scale=Sca, log=TRUE)) 
    rho    <- (runif(1)<exp(test))    
    m      <- cand*rho + m*(1-rho)
    m.save <- c(m.save,m) 
     if (M %% 100 == 0) print(paste("finished iteration",M))
}

   beta.post <- apply(beta[(round(num.reps/2)+1):(num.reps+1),],2,mean)
   beta.sd   <- apply(beta[(round(num.reps/2)+1):(num.reps+1),],2,sd)
   psi.hat   <- apply(psi.post[(round(num.reps/2)+1):(num.reps),],2,mean)
   p.hat     <- c(X %*% beta.post) + c(psi.hat)
   suppressWarnings(Y.hat     <- c(pnorm(p.hat)))
   resid     <- c(Y - Y.hat)
   
	    
	} ## END: IF BINOMIAL WITH PROBIT LINK FUNCTION IS TRUE (BEGIN FROM LINE 353)
	
	## BEGINING: IF NOT GAUSSIAN WITH IDENTITY LINK AND IF NOT BINOMIAL WITH PROBIT LINK
	else {
		print(family)
        stop("'family' is not recognized; only 'gaussian(link=identity)' and 'binomial(link=probit)' are used now.")
		} 	## END: IF NOT GAUSSIAN WITH IDENTITY LINK AND IF NOT BINOMIAL WITH PROBIT LINK (BEGIN FROM LINE 487)

	} # END: IF GAUSSIAN WITH IDENTITY LINK FUNCTION IS NOT TRUE (BEGIN FROM LINE 352)

   if (sum(is.finite(m.save)) == num.reps)  {	# JG (6/1/11): IN CASE VALUES OF m.save GO TO INFINITY
      K.est     <- c()
      for (l in 1:length(m.save)) {
   	  K.est[l] <- m.save[l]*sum(1/seq(m.save[l],m.save[l]+n-1,by=1))	
      }
      K.post   <- mean(K.est[round(length(m.save)/2):length(m.save)])
      K.sd     <- sd(K.est[round(length(m.save)/2):length(m.save)])
   }
   else {
      K.post   <- mean(K.save[round(length(m.save)/2):length(m.save)])
      K.sd     <- sd(K.save[round(length(m.save)/2):length(m.save)])
	  print(K.post)
   }                  # (THIS IS USED IN glmdm.linear; THE ONE USED IN glmdm.probit DOES NOT HAVE "if" FUNCTION (08/16/2011))
   
   #out <- list(coefficients= beta.post, standard.error = beta.sd, fitted.values = Y.hat, residuals = resid, K.est = K.post, K.sd = K.sd)   
  #return(out)
   
   out <- list(coefficients= beta.post, standard.error = beta.sd, fitted.values = Y.hat, residuals = resid, K.est = K.post, K.sd = K.sd, call=call, formula=formula, family=family) # ADD CALL, FORMULA, AND FAMILY IN THE OUTPUT (08/15/2011)
   
   

  graph.summary <- function(in.object, alpha = 0.05, digits=3, ...) # DESIGNED BY NEAL, CODED BY JEFF
{
    lo <- in.object$coefficient - qnorm(1-alpha/2) * in.object$standard.error
    hi <- in.object$coefficient + qnorm(1-alpha/2) * in.object$standard.error
    out.mat <- round(cbind(in.object$coefficient, in.object$standard.error, lo, hi),digits)
    blanks <- "                                                                          "
    dashes <- "--------------------------------------------------------------------------"
    bar.plot <- NULL
    scale.min <- floor(min(out.mat[,3])); scale.max <- ceiling(max(out.mat[,4]))
    for (i in 1:nrow(out.mat))  {
        ci.half.length <- abs(out.mat[i,1]-out.mat[i,3])
        ci.start <- out.mat[i,1] - ci.half.length
        ci.stop <- out.mat[i,1] + ci.half.length
        bar <- paste("|",substr(dashes,1,ci.half.length), "o", substr(dashes,1,ci.half.length), "|", sep="", collapse="")
        start.buf <- substr(blanks,1,round(abs(scale.min - ci.start)))
        stop.buf <- substr(blanks,1,round(abs(scale.max - ci.stop)))
        bar.plot <- rbind( bar.plot, paste(start.buf,bar, stop.buf, sep="", collapse="") )
    }
    rr     <- in.object$residuals
    ff     <- in.object$fitted.values
    mss    <- sum((ff - mean(ff))^2)
    rss    <- sum((rr - mean(rr))^2)
    n      <- length(rr)
    p      <- length(lo)-1
    df.int <- 1
    rdf    <- n-p-df.int
    resvar <- rss/rdf
    se            <- sqrt(resvar)
    r.squared     <- mss/(mss + rss)
    adj.r.squared <- 1 - (1 - r.squared) * ((n - df.int)/rdf)
    fstatistic    <- c((mss/(p - df.int))/resvar)
    f.p.val       <- 1-pf(fstatistic,p,rdf)
    
    out.df <- data.frame( matrix(NA,nrow=nrow(out.mat),ncol=ncol(out.mat)), bar.plot[1:length(bar.plot)] )
    out.df[1:nrow(out.mat),1:ncol(out.mat)] <- out.mat
    CI.label <- paste( "CIs:", substr(blanks,1,abs(scale.min)-2-4),"ZE+RO",
        substr(blanks,1,abs(scale.max)-2), sep="", collapse="" )
    dimnames(out.df)[[1]] <- names(in.object$coefficient)
    dimnames(out.df)[[2]] <- c("Coef","Std.Err.", paste(1-alpha,"Lower"),paste(1-alpha,"Upper"),CI.label)
    
    out.df <- list(call=in.object$call, coefficients=out.df)   # (08/15/2011)
    
    print(out.df)

    if (in.object$family$family %in% c("gaussian")) {  # BEGINING: IF GAUSSIAN FAMILY IS TRUE
    cat("\n")
    cat( paste("Residual standard error: ",round(se, digits),"on",rdf,"degrees of freedom\n") )
    cat( paste("Multiple R-squared: ",round(r.squared, digits),"  Adjusted R-squared:",round(adj.r.squared, digits),"\n") )
    cat( paste("F-statistic: ",round(fstatistic,2),"on",p,"and",rdf,"DF", " p-value:", round(f.p.val, digits), "\n") )
    cat( paste("Estimated K: ",round(in.object$K.est,2),"  Std.Err. of K:", round(in.object$K.sd, 2), "\n") )
    cat("\n")
    cat( paste("Residual Sum of Squares: ",round(rss, digits),"\n") )    	
    }                # END: IF GAUSSIAN FAMILY IS TRUE
    
    else {                                # BEGINING: IF GAUSSIAN FAMILY IS NOT TRUE
    cat("\n")
    cat( paste("Estimated K: ",round(in.object$K.est,2),"Std.Err. of K", round(in.object$K.sd, 2), "\n") )
    cat("\n")    	
    }                # END: IF GAUSSIAN FAMILY IS NOT TRUE
    
}   # END TO graph.summary FUNCTION CALL
  graph.summary(out)
#  return(out)


} # END TO glmdmEST FUNCTION CALL



#########################################################################
# CREATE A GENERIC FUNCTION glmdm (08/03/2011)
#########################################################################
glmdm <- function(formula, family=gaussian, data, num.reps=1000, a1=3, b1=2, d=0.25, MM=15, VV=30, ...) {
	mf <- model.frame(formula=formula, data=data)
	x <- model.matrix(attr(mf, "terms"), data=mf)[,-1]  # REMOVE THE INTERCEPT TERM SINCE IT IS ADDED IN LINE 180.
	y <- model.response(mf)
	
	est <- glmdmEST(Y=y, X=x, family=family, data=data, num.reps=num.reps, a1=a1, b1=b1, d=d, MM=MM, VV=VV, ...)
	est$call <- match.call()
	est$formula <- formula
	class(est) <- "glmdm"
	est	
}


