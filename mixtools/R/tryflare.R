try.flare <- function (y, x, lambda = NULL, beta = NULL, sigma = NULL, alpha = NULL, nu=1,
		epsilon = 1e-04, maxit = 10000, verb = FALSE, restart=50) 
{

loglik <- function(res, sigma, lambda, alpha) {
  tmp <- lambda*dnorm(res,sd=sqrt(sigma)) + (1-lambda)*dexp(res,rate=alpha)
  sum(log(tmp))
}



Q=function(res, sigma, lambda, alpha, z) {  
  Q <- sum(z*log(lambda)) + sum((1-z)*log(1-lambda)) -
       log(2*pi*sigma)*sum(z)/2 - sum(z*res^2)/2/sigma +
       log(alpha)*sum(1-z) - alpha*sum((1-z)*res)
  Q
}


Z <- function(res, sigma, lambda, alpha) {
  z=rep(1, length(res))
  z[res>0] = lambda / (lambda+(1-lambda)* sqrt(2*pi*sigma) * alpha *
                       as.vector(exp(res[res>0]^2/2/sigma - alpha*res[res>0])))
  z
}


    x <- cbind(1, x)
    n <- length(y)
    p <- ncol(x)

    est <- flaremix.init(y=y, x=x, lambda=lambda, beta=beta, sigma=sigma, alpha=alpha)
    lambda <- est$lambda
    beta <- est$beta
    sigma <- est$sigma
    alpha <- est$alpha


    diff <- 1
    iter <- 0
    counts <- 0
	ll.counts<-0
    xbeta <- x %*% beta
    res <- y - xbeta


	dn <- dnorm(res,sd=sqrt(sigma))
	de <- dexp(res,rate=alpha)

	obsloglik <- loglik(res, sigma, lambda, alpha)
	ll<-obsloglik

  Q1 <- -Inf
  all.Q <- NULL

	z=Z(res,sigma,lambda,alpha)

    while (sum(abs(diff) > epsilon)>0 && iter < maxit) {


      iter=iter+1


temp=(solve(-1/sigma*t(x) %*% sweep(x, 1, z, "*") + nu*t(x) %*% sweep(x, 1, (1-z)/(y-x%*%beta), "*"))%*%(1/sigma*(apply(sweep(x,1,z*(y-x%*%beta),"*"),2,sum)) +alpha*apply(sweep(x,1,1-z,"*"),2,sum)))

m=1
while(m<restart){
j=1
while(j<restart){

beta.new <- try(beta - (1/2^j)*temp,silent=TRUE)

if(class(beta.new)=="try-error" || sum(is.na(beta.new))>0) beta.new=beta else beta.new=beta.new

	  xbeta.new <- x %*% beta.new
	  res.new <- y-xbeta.new


Q.beta <- Q(res.new,sigma,lambda,alpha,z)


if(Q.beta < Q1) j=j+1 else j=101

}

if(j==restart) stop(paste("Too many attempts at step-halving!","\n"))


	z.new=Z(res.new,sigma,lambda,alpha)


      lambda.new <- mean(z.new)


    sigma.new <- sum(z.new*(res.new^2))/sum(z.new)

  alpha.new <- sum(1-z.new[res.new>0])/sum((1-z.new[res.new>0])*res.new[res.new>0])




diff <- c(lambda.new,beta.new,sigma.new,alpha.new)-c(lambda,beta,sigma,alpha)

  	  z.new2=Z(res,sigma,lambda,alpha)
#		z.new2=z.new2/apply(z.new2,1,sum)
	  Q.new <- Q(res.new,sigma.new,lambda.new,alpha.new,z.new2)
q.diff=Q.new-Q1

if(q.diff<0) m=m+1 else m=101

}

if(m==restart) stop(paste("Too many attempts at step-halving!","\n"))

	  lambda <- lambda.new
	  beta <- beta.new
	  xbeta <- xbeta.new
	  res <- res.new
	  sigma <- sigma.new
	  alpha <- alpha.new
	  z<-z.new2

		newobsloglik <- loglik(res.new, sigma.new, lambda.new, alpha.new)
	ll<-c(ll, newobsloglik)


	  counts <- counts + (Q.new<Q1)
	  all.Q <- c(all.Q,Q.new)
		ll.counts <- ll.counts + (obsloglik>newobsloglik)




	  Q1 <- Q.new
	


            obsloglik <- newobsloglik



	if(verb==TRUE) cat("iteration=", iter, "diff=", diff, "log-likelihood", obsloglik, "\n")
	
	}

    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
#	par(mfrow=c(2,1))
#    plot(all.Q,type="l")
#    plot(ll,type="l")
    cat("number of iterations=", iter, "\n")
    a=list(x=x,y=y,posterior=cbind(z,1-z),lambda=c(lambda,1-lambda),beta=beta,sigma=sigma,alpha=alpha,loglik=obsloglik,all.loglik=ll,ft="flaremixEM")
	class(a)="mixEM"
	a

}


