regmix.mixed.init <- function(y,x,w=NULL,sigma=NULL,arb.sigma=TRUE,
	alpha=NULL,lambda=NULL,mu=NULL,R=NULL,arb.R=TRUE,k=2,mixed=FALSE,
	addintercept.fixed=FALSE,addintercept.random=TRUE){

    N <- length(y)
    n <- sapply(y, length)
    p <- ncol(x[[1]])

    if (addintercept.random) {
        x = lapply(1:N, function(i) as.matrix(x[[i]][,-1]))
    } else x=x

    if(is.null(w)==FALSE && sum(sapply(w,sum))!=0){
    q <- ncol(w[[1]])
    }

    if (mixed == TRUE && is.null(alpha)==TRUE) {
        if (addintercept.fixed) {
            w.1 = list()
		w.1 = lapply(1:N,function(i) w[[i]][,-1])
		lm.out = lapply(1:N,function(i) lm(y[[i]]~w.1[[i]]))
		alpha.hyp = apply(matrix(sapply(lm.out,coef),ncol=N),1,mean)
		sd.hyp = lapply(lm.out,anova)
		sd.hyp = mean(as.vector(sqrt(sapply(1:N,function(i) sd.hyp[[i]]$Mean[length(sd.hyp[[i]]$Mean)]))))
		alpha = rnorm(q, mean=alpha.hyp, sd=sd.hyp)
        } else {
            w.1 = w
		lm.out = lapply(1:N,function(i) lm(y[[i]]~w.1[[i]]-1))
		alpha.hyp = apply(matrix(sapply(lm.out,coef),ncol=N),1,mean)
		sd.hyp = lapply(lm.out,anova)
		sd.hyp = mean(as.vector(sqrt(sapply(1:N,function(i) sd.hyp[[i]]$Mean[length(sd.hyp[[i]]$Mean)]))))
		alpha = rnorm(q, mean=alpha.hyp, sd=sd.hyp)
	  }
    }
    if(mixed==FALSE) {
        alpha = 0
    }

y.x = lapply(1:N, function(i) cbind(y[[i]],x[[i]]))
a = order(sapply(1:N, function(i) mean(y[[i]])))
y.x = lapply(1:N,function(i) y.x[[a[i]]])

	y.x.bin.list = list()
	y.x.bin = list()
	for(j in 1:k){
		y.x.bin.list[[j]] <- y.x[max(1,floor((j-1)*N/k)):ceiling(j*N/k)]
		y.x.2 <- NULL
		for(i in 1:length(y.x.bin.list[[j]])){
		y.x.2 <- rbind(y.x.2,y.x.bin.list[[j]][[i]])
		}
		y.x.bin[[j]] <- y.x.2
	}


if(addintercept.random){
lm.out <- lapply(1:k, function(i) lm(y.x.bin[[i]][,1]~y.x.bin[[i]][,2:p]))
lm.out.beta <- lapply(1:k, function(j) lapply(1:length(y.x.bin.list[[j]]), function(i) lm(y.x.bin.list[[j]][[i]][,1]~y.x.bin.list[[j]][[i]][,2:p])))
beta <- lapply(1:k,function(j) matrix(sapply(lm.out.beta[[j]],coef),nrow=p))
} else {
lm.out <- lapply(1:k, function(i) lm(y.x.bin[[i]][,1]~y.x.bin[[i]][,2:(p+1)]-1))
lm.out.beta <- lapply(1:k, function(j) lapply(1:length(y.x.bin.list[[j]]), function(i) lm(y.x.bin.list[[j]][[i]][,1]~y.x.bin.list[[j]][[i]][,2:(p+1)]-1)))
beta <- lapply(1:k,function(j) matrix(sapply(lm.out.beta[[j]],coef),nrow=p))
}
	if(is.null(sigma)) {
		sigma.hyp = lapply(lm.out,anova)
		sigma.hyp = as.vector(sqrt(sapply(1:k,function(i) sigma.hyp[[i]]$Mean[length(sigma.hyp[[i]]$Mean)])))
		if(arb.sigma) {
			sigma=1/rexp(k,rate=sigma.hyp)
		} else {
			sigma.hyp=mean(sigma.hyp)
			sigma=1/rexp(1,rate=sigma.hyp)
		}
	}

	if(is.null(sigma)==FALSE && arb.sigma==TRUE) {
		k=length(sigma)
	}

	if(is.null(R)) {
		if(arb.R) {
			R.hyp = lapply(1:k,function(i) (apply(beta[[i]],1,var))^-1)
			R=lapply(1:k,function(i) diag(1/rexp(p,rate=R.hyp[[i]]),p))
		} else {
			R.hyp = apply(matrix(sapply(1:k,function(i) (apply(beta[[i]],1,var))^-1),ncol=k),2,mean)
			R = diag(1/rexp(p,rate=sigma.hyp),p)
	}
}
	if(is.null(R)==FALSE && arb.R==TRUE) {
		k=length(R)
	}


    if (is.null(mu)) {
        mu.hyp = lapply(1:k,function(i) apply(beta[[i]],1,mean))	  
	  mu = matrix(ncol=k,nrow=p)
		if(arb.R==TRUE){
		  for(j in 1:k){
		  mu[,j] = rmvnorm(1,mu=as.vector(mu.hyp[[j]]),sigma=R[[j]])
			} 
		} else {
		  for(j in 1:k){
		  mu[,j] = rmvnorm(1,mu=as.vector(mu.hyp[[j]]),sigma=R)
		} 
    }
} else k = ncol(mu)


    if (is.null(lambda)) {
        lambda = runif(k)
        lambda = lambda/sum(lambda)
    } else k = length(lambda)

list(sigma=sigma,alpha=alpha,lambda=lambda,mu=mu,R=R,k=k)

}