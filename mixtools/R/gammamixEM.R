gammamixEM <-
function (x, lambda = NULL, alpha = NULL, beta = NULL, k = 2,  
          epsilon = 1e-08, maxit = 1000, maxrestarts=20, verb = FALSE) {
  x <- as.vector(x)
 tmp <- gammamix.init(x = x, lambda = lambda, alpha=alpha, beta=beta, k=k)
  lambda <- tmp$lambda 
  alpha <- tmp$alpha 
  beta <- tmp$beta
  theta <- c(alpha,beta)
  k <- tmp$k 
  iter <- 0
  mr <- 0
  diff <- epsilon+1
  n <- length(x)
  dens <- NULL
  dens <- function(lambda, theta, k){
	temp<-NULL
	alpha=theta[1:k]
	beta=theta[(k+1):(2*k)]
	for(j in 1:k){
	 temp=cbind(temp,dgamma(x,shape=alpha[j],scale=beta[j]))  }
	temp=t(lambda*t(temp))
	temp
	}
  old.obs.ll <- sum(log(apply(dens(lambda, theta, k),1,sum)))
  ll <- old.obs.ll
  gamma.ll <- function(theta, z,lambda, k) -sum(z*log(dens(lambda,theta,k)))
  while(diff > epsilon && iter < maxit){
	dens1=dens(lambda,theta,k)
	z=dens1/apply(dens1,1,sum)
	lambda.hat=apply(z,2,mean)
	out=try(suppressWarnings(nlm(gamma.ll,p=theta,lambda=lambda.hat,k=k,z=z)),
          silent=TRUE)
		if(class(out)=="try-error"){
			cat("Note: Choosing new starting values.", "\n")
			if(mr==maxrestarts) stop(paste("Try different number of components?","\n"))
			mr <- mr+1
  			tmp <- gammamix.init(x = x, k=k)
  			lambda <- tmp$lambda 
  			alpha <- tmp$alpha 
  			beta <- tmp$beta
  			theta <- c(alpha,beta)
  			k <- tmp$k 
  			iter <- 0
  			diff <- epsilon+1
  			old.obs.ll <- sum(log(apply(dens(lambda, theta, k),1,sum)))
  			ll <- old.obs.ll
  		} else{
	theta.hat=out$estimate
	alpha.hat=theta.hat[1:k]
	beta.hat=theta.hat[(k+1):(2*k)]
	new.obs.ll <- sum(log(apply(dens(lambda.hat, theta.hat, k),1,sum)))
	diff <- new.obs.ll-old.obs.ll
	old.obs.ll <- new.obs.ll
	ll <- c(ll,old.obs.ll)
	lambda=lambda.hat
	theta=theta.hat
	alpha=alpha.hat
	beta=beta.hat
	iter=iter+1
         if (verb) {
            cat("iteration =", iter, " log-lik diff =", diff, " log-lik =", 
                new.obs.ll, "\n")
          }
	}
   }
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
theta=rbind(alpha,beta)
rownames(theta)=c("alpha","beta")
colnames(theta)=c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=x, lambda = lambda, gamma.pars = theta, loglik = new.obs.ll, 
        posterior = z, all.loglik=ll, ft="gammamixEM")
    class(a) = "mixEM"
    a
  }	