boot.se <- function (em.fit, B = 100,
    arbmean = TRUE, arbvar = TRUE, N=NULL, ...)
{
    mix.type <- em.fit$ft
    if (mix.type == "regmixEM") {


#Start Here...

                k=length(em.fit$lambda)
                y=em.fit$y
                n=length(y)
                if(sum(em.fit$x[,1]==1)==n){
                    x=em.fit$x[,-1]
                }
                else x=em.fit$x
                    if (arbmean == FALSE) {
                      scale = em.fit$scale
                      beta = matrix(rep(em.fit$beta, k), ncol = k)
                    }
                    else {
                      scale = 1
                      beta = em.fit$beta
                    }
                  
                  xbeta.new=em.fit$x%*%beta
                  lambda = em.fit$lambda
                  if(arbvar==FALSE){
                    sigma = rep(em.fit$sigma,k)
                  } else{
                    sigma = scale*em.fit$sigma
                    }

                  j = 0
                lambda.bs=NULL
                beta.bs=NULL
                sigma.bs=NULL
                scale.bs=NULL
                while (j < B) {
                  j = j + 1
                  w=rmultinom(n,size=1,prob=lambda)
                  y.sim=sapply(1:n,function(i) rnorm(1,mean=xbeta.new[i,(w[,i]==1)],sd=sigma[w[,i]==1]) )
                  em.bs = try(regmixEM(y = y.sim, x = x, k = k,
                    arbmean = arbmean, arbvar = arbvar, lambda=em.fit$lambda, beta=em.fit$beta,
                    sigma=(scale*em.fit$sigma), ...), silent = TRUE)
                  if (class(em.bs) == "try-error" || em.bs$restarts!=0) {
                    j = j - 1
                  }

                  else {
                  if(arbmean==FALSE){
                    lambda.bs = cbind(lambda.bs,em.bs$lambda)
                    beta.bs = cbind(beta.bs,as.vector(em.bs$beta))
                    sigma.bs = cbind(sigma.bs,em.bs$sigma)
                    scale.bs = cbind(scale.bs,em.bs$scale)
                    }
                  else {
                    lambda.bs = cbind(lambda.bs,em.bs$lambda)
                    beta.bs = cbind(beta.bs,as.vector(em.bs$beta))
                    sigma.bs = cbind(sigma.bs,em.bs$sigma)
                  }
                  }

                }
    if(arbmean==FALSE){
    lambda.se=apply(lambda.bs,1,sd)
    beta.se=apply(beta.bs,1,sd)
    sigma.se=apply(sigma.bs,1,sd)
    scale.se=apply(scale.bs,1,sd)
    bs.list=list(lambda = lambda.bs, lambda.se = lambda.se, beta = beta.bs, beta.se = beta.se, 
    sigma=sigma.bs, sigma.se = sigma.se, scale=scale.bs, scale.se=scale.se)
    }
    else{
    lambda.se=apply(lambda.bs,1,sd)
    beta.se=matrix(apply(beta.bs,1,sd),ncol=k)
    sigma.se=apply(sigma.bs,1,sd)
    bs.list=list(lambda = lambda.bs, lambda.se = lambda.se, beta = beta.bs, beta.se = beta.se, 
    sigma=sigma.bs, sigma.se = sigma.se)
   }
}

#####
   if (mix.type == "regmixEM.mixed"){




#Start Here...

                k=length(em.fit$lambda)
                y=em.fit$y
		    x=em.fit$x
        if (length(y) != length(x))
            stop("Number of elements in lists for x and y must match!")
		    w=em.fit$w
		    p.z=em.fit$posterior.z
		    p.b=em.fit$posterior.beta
		    mu=em.fit$mu
		    R=em.fit$R
		    sigma=em.fit$sigma
		    lambda=em.fit$lambda
		    n=length(y)
		    p=nrow(mu)
		    n.i=sapply(y,length)
		    alpha=em.fit$alpha
                    
                    w.test <- NULL
			for(i in 1:length(w)){
				w.test <- c(w.test,as.vector(w[[i]]))
				}
		    if(sum(w.test)==0){
			 w1=NULL} else w1=w


		    if(length(R)==k) arb.R=TRUE else arb.R=FALSE
		    if(length(sigma)==k) arb.sigma=TRUE else arb.sigma=FALSE


                j = 0
                lambda.bs=NULL
                mu.bs=NULL
                sigma.bs=NULL
                R.bs=NULL
		    alpha.bs=NULL

                while (j < B) {
                  j = j + 1


		    k.bs=sapply(1:n, function(i) rmultinom(1,size=1,prob=lambda))
		    bs.i <- sample(1:n,n,replace=TRUE)
##
		    if(arb.R==FALSE){
			if(arb.sigma==FALSE){
		    		y.sim <- lapply(1:n, function(i) w[[i]] %*% 
                  alpha + as.vector(rmvnorm(1, mu = x[[i]] %*% mu[, 
                  (k.bs[, i] == 1)], sigma = (x[[i]]%*%R%*%t(x[[i]])  +diag(sigma,n.i[i])) )))
				} else y.sim <- lapply(1:n, function(i) w[[i]] %*% 
                  alpha + as.vector(rmvnorm(1, mu = x[[i]] %*% mu[, 
                  (k.bs[, i] == 1)], sigma = (x[[i]]%*%R%*%t(x[[i]])  +diag(sigma[(k.bs[,i] == 1)],n.i[i])) )))
		    } else{
				if(arb.sigma==FALSE){
					y.sim <- lapply(1:n, function(i) w[[i]] %*% 
                  alpha + as.vector(rmvnorm(1, mu = x[[i]] %*% mu[, 
                  (k.bs[, i] == 1)], sigma = (x[[i]]%*%R[(k.bs[, i] == 
                  1)][[1]]%*%t(x[[i]])  +diag(sigma,n.i[i])) )))   
					} else y.sim <- lapply(1:n, function(i) w[[i]] %*% 
                  alpha + as.vector(rmvnorm(1, mu = x[[i]] %*% mu[, 
                  (k.bs[, i] == 1)], sigma = (x[[i]]%*%R[(k.bs[, i] == 
                  1)][[1]]%*%t(x[[i]])  +diag(sigma[(k.bs[,i] == 1)],n.i[i])) )))   
			}

##

                  em.bs = try(regmixEM.mixed(y = y.sim, x = x, w=w1, sigma=sigma, mu=mu, alpha=alpha, R=R, 
					lambda=lambda, k = k, arb.R=arb.R, arb.sigma=arb.sigma, addintercept.fixed=FALSE,
					addintercept.random=FALSE, ...), silent = TRUE)
                  if (class(em.bs) == "try-error" || em.bs$restarts!=0) {
                    j = j - 1
                  }

                  else {

			lambda.bs = cbind(lambda.bs,em.bs$lambda)
			alpha.bs = cbind(alpha.bs,em.bs$alpha)
			sigma.bs = cbind(sigma.bs,em.bs$sigma)
			mu.bs = cbind(mu.bs,as.vector(em.bs$mu))
			R.bs = cbind(R.bs,as.vector(sapply(em.bs$R,c)))
			}


lambda.se=apply(lambda.bs,1,sd)
alpha.se=apply(alpha.bs,1,sd)
sigma.se=apply(sigma.bs,1,sd)
mu.se=matrix(apply(mu.bs,1,sd),ncol=k)
R.se1=apply(R.bs,1,sd)
	if(arb.R==TRUE){
		R.se2=matrix(R.se1,ncol=k)
		R.se=lapply(1:k,function(i) matrix(R.se2[,i],ncol=p))
		} else R.se=matrix(R.se1,ncol=p)

    bs.list=list(lambda = lambda.bs, lambda.se = lambda.se, alpha=alpha.bs, alpha.se=alpha.se,
			mu = mu.bs, mu.se = mu.se, sigma=sigma.bs, sigma.se = sigma.se, R=R.bs, R.se=R.se)

}



}


#####
# Commented out by DRH on 8-29-2008 due to absence of gammamixEM function
#    if (mix.type == "gammamixEM") {
#        x <- em.fit$x
#        n <- length(x)
#        k <- length(em.fit$lambda)
#        alpha <- em.fit$gamma.pars[1,]
#        beta <- em.fit$gamma.pars[2,]
#        lambda <- em.fit$lambda
#        j = 0
#        lambda.bs = NULL
#	  alpha.bs = NULL
#        beta.bs = NULL
#        while (j < B) {
#            j = j + 1
#		comp = sample(1:k,size=n,replace=T,prob=lambda)
#		x.sim = sapply(1:n,function(i) rgamma(1,shape=alpha[comp[i]],scale=beta[comp[i]]))
#            em.bs = try(gammamixEM(x = x.sim, k = k, lambda = lambda, 
#                alpha = alpha, beta = beta, ...), silent = TRUE)
#            if (class(em.bs) == "try-error") {
#                j = j - 1
#            }
#            else {
#                lambda.bs = cbind(lambda.bs, em.bs$lambda)
#                alpha.bs = cbind(alpha.bs, as.vector(em.bs$gamma.pars[1,]))
#                beta.bs = cbind(beta.bs, as.vector(em.bs$gamma.pars[2,]))
#            }
#        }
#        lambda.se = apply(lambda.bs, 1, sd)
#	  alpha.se = apply(alpha.bs, 1, sd)
#	  beta.se = apply(beta.bs, 1, sd)
#        bs.list = list(lambda = lambda.bs, lambda.se = lambda.se, 
#            alpha = alpha.bs, alpha.se = alpha.se, beta = beta.bs, beta.se = beta.se)
#    }

#####
    if (mix.type == "repnormmixEM") {


#Start Here...

                k=length(em.fit$lambda)
                y=em.fit$y
                m=nrow(y)
        n=ncol(y)

                    if (arbmean == FALSE) {
                      scale = em.fit$scale
                      mu = rep(em.fit$mu, k)
                    }
                    else {
                      scale = 1
              mu = em.fit$mu
                    }

                  lambda = em.fit$lambda
                  if(arbvar==FALSE){
                    sigma = rep(em.fit$sigma,k)
                  } else{
                    sigma = scale*em.fit$sigma
                    }

                  j = 0
                lambda.bs=NULL
                mu.bs=NULL
                sigma.bs=NULL
                scale.bs=NULL
                while (j < B) {
                  j = j + 1
                  w=rmultinom(n,size=1,prob=lambda)
                  y.sim=sapply(1:n,function(i) rnorm(m,mean=mu[w[,i]==1],sd=sigma[w[,i]==1]) )
                  em.bs = try(repnormmixEM(x = y.sim,k = k,
                    arbmean = arbmean, arbvar = arbvar, lambda=em.fit$lambda, mu=em.fit$mu,
                    sigma=(scale*em.fit$sigma), ...), silent = TRUE)
                  if (class(em.bs) == "try-error" || em.bs$restarts!=0) {
                    j = j - 1
                  }

                  else {
                  if(arbmean==FALSE){
                    lambda.bs = cbind(lambda.bs,em.bs$lambda)
                    mu.bs = cbind(mu.bs,as.vector(em.bs$mu))
                    sigma.bs = cbind(sigma.bs,em.bs$sigma)
                    scale.bs = cbind(scale.bs,em.bs$scale)
                    }
                  else {
                    lambda.bs = cbind(lambda.bs,em.bs$lambda)
                    mu.bs = cbind(mu.bs,as.vector(mu.bs$beta))
                    sigma.bs = cbind(sigma.bs,em.bs$sigma)
                  }
                  }

                }
    if(arbmean==FALSE){
    lambda.se=apply(lambda.bs,1,sd)
    mu.se=apply(mu.bs,1,sd)
    sigma.se=apply(sigma.bs,1,sd)
    scale.se=apply(scale.bs,1,sd)
    bs.list=list(lambda = lambda.bs, lambda.se = lambda.se, mu = mu.bs, mu.se = mu.se, 
    sigma=sigma.bs, sigma.se = sigma.se, scale=scale.bs, scale.se=scale.se)
    }
    else{
    lambda.se=apply(lambda.bs,1,sd)
    mu.se=apply(mu.bs,1,sd)
    sigma.se=apply(sigma.bs,1,sd)
    bs.list=list(lambda = lambda.bs, lambda.se = lambda.se, mu = mu.bs, mu.se = mu.se, 
    sigma=sigma.bs, sigma.se = sigma.se)
   }
}

####    
    if (mix.type == "mvnormalmixEM") {


#Start Here...

                k=length(em.fit$lambda)
                y=em.fit$x
                n=nrow(y)
        p=ncol(y)

                    if (arbmean == FALSE) {
                      mu = lapply(1:k,function(i) em.fit$mu)
                    }
                    else {
              mu = em.fit$mu
                    }

                  lambda = em.fit$lambda

                  if(arbvar==FALSE){
                    sigma = lapply(1:k, function(i) em.fit$sigma)
                  } else{
                    sigma = em.fit$sigma
                    }

                  j = 0
                lambda.bs=NULL
                mu.bs=NULL
                sigma.bs=NULL

                while (j < B) {
                  j = j + 1
                  w=rmultinom(n,size=1,prob=lambda)
                  y.sim=t(sapply(1:n,function(i) rmvnorm(1,mu=mu[w[,i]==1][[1]],sigma=sigma[w[,i]==1][[1]]) ))
                  em.bs = try(mvnormalmixEM(x = y.sim, k = k,
                    arbmean = arbmean, arbvar = arbvar, lambda=em.fit$lambda, mu=em.fit$mu,
                    sigma=em.fit$sigma, ...), silent = TRUE)
                  if (class(em.bs) == "try-error" || em.bs$restarts!=0) {
                    j = j - 1
                  }

                  else {

                    lambda.bs = cbind(lambda.bs,em.bs$lambda)
                    mu.bs1 = as.vector(sapply(em.bs$mu,c))
            mu.bs = cbind(mu.bs, mu.bs1)
                    if(arbvar==FALSE){
                    sigma.bs=cbind(sigma.bs,as.vector(em.bs$sigma))
                    } else{
                    sigma.bs1 = lapply(1:k, function(i) as.vector(em.bs$sigma[[i]]))
            sigma.bs2 = as.vector(sapply(sigma.bs1,c))
            sigma.bs = cbind(sigma.bs,sigma.bs2)
                  }
                  }

                }

    lambda.se=apply(lambda.bs,1,sd)
    mu.se1=apply(mu.bs,1,sd)
    if(arbmean==TRUE){
    mu.se = lapply(1:k,function(i) mu.se1[((i-1)*p+1):(i*p)])
    } else mu.se = mu.se1
    sigma.se1=apply(sigma.bs,1,sd)
    if(arbvar==TRUE){
    sigma.se=lapply(1:k, function(i) matrix(sigma.se1[((i-1)*(p^2)+1):(i*(p^2))], nrow=p,ncol=p))
    } else sigma.se=matrix(sigma.se1,nrow=p)
    bs.list=list(lambda = lambda.bs, lambda.se = lambda.se, mu = mu.bs, mu.se = mu.se, 
    sigma=sigma.bs, sigma.se = sigma.se)
 
}

####
    if (mix.type == "normalmixEM") {


#Start Here...

                k=length(em.fit$lambda)
                y=em.fit$x
                n=length(y)
                x=em.fit$x
                    if (arbmean == FALSE) {
                      scale = em.fit$scale
                      mu = rep(em.fit$mu, k)
                    }
                    else {
                      scale = 1
                      mu = em.fit$mu
                    }
                  
                  lambda = em.fit$lambda
                  if(arbvar==FALSE){
                    sigma = rep(em.fit$sigma,k)
                  } else{
                    sigma = scale*em.fit$sigma
                    }

                  j = 0
                lambda.bs=NULL
                mu.bs=NULL
                sigma.bs=NULL
                scale.bs=NULL
                while (j < B) {
                  j = j + 1
                  w=rmultinom(n,size=1,prob=lambda)
                  y.sim=sapply(1:n,function(i) rnorm(1,mean=mu[(w[,i]==1)],sd=sigma[w[,i]==1]) )
                  em.bs = try(normalmixEM(x = y.sim, k = k,
                    arbmean = arbmean, arbvar = arbvar, lambda=em.fit$lambda, mu=em.fit$mu,
                    sigma=(scale*em.fit$sigma), ...), silent = TRUE)
                  if (class(em.bs) == "try-error" || em.bs$restarts!=0) {
                    j = j - 1
                  }

                  else {
                  if(arbmean==FALSE){
                    lambda.bs = cbind(lambda.bs,em.bs$lambda)
                    mu.bs = cbind(mu.bs,as.vector(em.bs$mu))
                    sigma.bs = cbind(sigma.bs,em.bs$sigma)
                    scale.bs = cbind(scale.bs,em.bs$scale)
                    }
                  else {
                    lambda.bs = cbind(lambda.bs,em.bs$lambda)
                    mu.bs = cbind(mu.bs,as.vector(em.bs$mu))
                    sigma.bs = cbind(sigma.bs,em.bs$sigma)
                  }
                  }

                }
    if(arbmean==FALSE){
    lambda.se=apply(lambda.bs,1,sd)
    mu.se=apply(mu.bs,1,sd)
    sigma.se=apply(sigma.bs,1,sd)
    scale.se=apply(scale.bs,1,sd)
    bs.list=list(lambda = lambda.bs, lambda.se = lambda.se, mu = mu.bs, mu.se = mu.se, 
    sigma=sigma.bs, sigma.se = sigma.se, scale=scale.bs, scale.se=scale.se)
    }
    else{
    lambda.se=apply(lambda.bs,1,sd)
    mu.se=matrix(apply(mu.bs,1,sd),ncol=k)
    sigma.se=apply(sigma.bs,1,sd)
    bs.list=list(lambda = lambda.bs, lambda.se = lambda.se, mu = mu.bs, mu.se = mu.se, 
    sigma=sigma.bs, sigma.se = sigma.se)
   }
}

###
if (mix.type == "multmixEM") {


y<-em.fit$y
n<-nrow(y)
n.i<-apply(y,1,sum)
p<-ncol(y)
k<-length(em.fit$lambda)
theta<-em.fit$theta
lambda<-em.fit$lambda


                  j = 0
                lambda.bs=NULL
                theta.bs=NULL


                while (j < B) {
                  j = j + 1

                  w=rmultinom(n,size=1,prob=lambda)
                  y.sim=t(sapply(1:n,function(i) rmultinom(1,size=n.i[i],prob=theta[(w[,i]==1),]) ))

                  em.bs = try(multmixEM(y = y.sim, k = k, lambda=lambda, theta=theta, ...), silent = TRUE)
                  if (class(em.bs) == "try-error" || em.bs$restarts!=0) {
                    j = j - 1
                  }
        
            else{
            lambda.bs=cbind(lambda.bs,em.bs$lambda)
            theta.bs=cbind(theta.bs,as.vector(em.bs$theta))
            }

        }


        lambda.se=apply(lambda.bs,1,sd)
        theta.se=matrix(apply(theta.bs,1,sd),nrow=k)
        bs.list=list(lambda=lambda.bs, lambda.se=lambda.se, theta=theta.bs, theta.se=theta.se)

}

    
####



if (mix.type == "logisregmixEM") {

        y=em.fit$y
	n=length(y)
        if (is.null(N)) N=rep(1,n)
                k=length(em.fit$lambda)

                
                if(sum(em.fit$x[,1]==1)==n){
                    x=em.fit$x[,-1]
                } else x=em.fit$x
        logit <- function(x) 1/(1 + exp(-x))


        lambda<-em.fit$lambda
        beta<-em.fit$beta
        xbeta.new=em.fit$x%*%beta
        prob=logit(xbeta.new)

                  j = 0
                lambda.bs=NULL
                beta.bs=NULL


                while (j < B) {
                  j = j + 1

                  w=rmultinom(n,size=1,prob=lambda)

                  y.sim = sapply(1:n, function(i) rbinom(1, size=N[i], prob=prob[,(w[,i]==1)]))


                  em.bs = try(logisregmixEM(y = y.sim, x=x, N=N, k = k, lambda=lambda, beta=beta,...), silent = TRUE)
                  if (class(em.bs) == "try-error" || em.bs$restarts!=0) {
                    j = j - 1
                  }
        
            else{
            lambda.bs=cbind(lambda.bs,em.bs$lambda)
            beta.bs=cbind(beta.bs,as.vector(em.bs$beta))
            }

        }


        lambda.se=apply(lambda.bs,1,sd)
        beta.se=matrix(apply(beta.bs,1,sd),nrow=k)
        bs.list=list(lambda=lambda.bs, lambda.se=lambda.se, beta=beta.bs, beta.se=beta.se)

    


}
####




if (mix.type == "poisregmixEM") {


                k=length(em.fit$lambda)
                y=em.fit$y
                n=length(y)
                if(sum(em.fit$x[,1]==1)==n){
                    x=em.fit$x[,-1]
                } else x=em.fit$x


        lambda<-em.fit$lambda
        beta<-em.fit$beta
        xbeta.new=em.fit$x%*%beta
        prob=exp(xbeta.new)

                j = 0
                lambda.bs=NULL
                beta.bs=NULL


                while (j < B) {
                  j = j + 1

                  w=rmultinom(n,size=1,prob=lambda)

                  y.sim = sapply(1:n, function(i) rpois(1, lambda=prob[,(w[,i]==1)]))


                  em.bs = try(poisregmixEM(y = y.sim, x=x, k = k, lambda=lambda, beta=beta, ...), silent = TRUE)
                  if (class(em.bs) == "try-error" || em.bs$restarts!=0) {
                    j = j - 1
                  }
        
            else{
            lambda.bs=cbind(lambda.bs,em.bs$lambda)
            beta.bs=cbind(beta.bs,as.vector(em.bs$beta))
            }

        }


        lambda.se=apply(lambda.bs,1,sd)
        beta.se=matrix(apply(beta.bs,1,sd),nrow=k)
        bs.list=list(lambda=lambda.bs, lambda.se=lambda.se, beta=beta.bs, beta.se=beta.se)

}


    bs.list
    }
