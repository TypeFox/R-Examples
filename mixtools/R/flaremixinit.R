flaremix.init <- function(y, x, lambda=NULL, beta=NULL, sigma=NULL, alpha=NULL){

n<-length(y)

if(is.null(lambda)){
	lambda=runif(2)
	lambda=lambda/sum(lambda)
}

lm.out=lm(y~x[,2])

if(is.null(beta)){
	beta=lm.out$coef
	beta[1]=beta[1]+mean(sort(lm.out$residuals)[(n-10):n])
	beta[2]=rnorm(1,mean=beta[2],sd=abs(beta[2]/10))
}

if(is.null(sigma)){
	sigma=rexp(1,rate=sqrt(1/anova(lm.out)$Mean[2]))
}

if(is.null(alpha)){
	a=1/sum(lm.out$res[lm.out$res>0])
	alpha=abs(rnorm(1,a))
}

list(lambda=lambda[1], beta=matrix(beta,ncol=1), sigma=sigma, alpha=alpha)

}
