bootcomp <- function(x,y,ncomp=2,ncincr=1,intercept=TRUE,nboot=1000,
		     ts1=NULL,ts2=NULL,sem.par=FALSE,verb=FALSE,
                     print.prog=TRUE,...)
{
#
# Function bootcomp to test for K versus K+INCR components in
# a mixture of linear regressions (with normal errors) via
# a bootstrap distribution for the likelihood ratio statistic.
# (NOTE:  This statistic DOES NOT have a chi-squared-nu
# distribution!!!) (Where --- irrelevantly --- nu is the number of
# parameters for any one component, equal to p+2 where
# p is the number of linear parameters in each component.)
#
# Screw-ups: 0 <--> bailed out, n components;
#            1 <--> didn't converge, n components;
#            2 <--> bailed out, n+ncincr components;
#            3 <--> didn't converge, n+ncincr components;
#            4 <--> lrs(n) > lrs(n+ncincr).

obj1 <- mixreg(x,y,ncomp=ncomp,intercept=intercept,
               theta.start=ts1,verb=verb,...)
obj2 <- mixreg(x,y,ncomp=ncomp+ncincr,intercept=intercept,
               theta.start=ts2,verb=verb,...)
theta1 <- obj1$theta
theta2 <- obj2$theta

if(sem.par) {
	m     <- matrix(unlist(theta1),ncol=length(theta1))
	nr    <- nrow(m)
	m     <- m[-c(nr-1,nr),]
	mu    <- if(intercept) cbind(1,x)%*%m else x%*%m
	resid <- y-mu
	prob  <- if(intercept)
			gfun(cbind(1,x),y,theta1)$gamma
		else
			gfun(x,y,theta1)$gamma
}

rslt <- list()
screw.ups <- list()
aic1 <- list()
aic2 <- list()
k <- 0
for(i in 1:nboot) {
	repeat {
		save.seed <- .Random.seed
		yboot <- if(sem.par) boot.sam(mu,resid,prob)
			        else simmix(theta1,intercept,x)$y
		tmp <- mixreg(x,yboot,ncomp=ncomp,theta.start=theta1,
                              intercept=intercept,verb=verb,...)
		l1 <- tmp$log.like
		a1 <- tmp$aic
		if(!tmp$converged) {
			k <- k+1
			if(is.na(l1)) screw.ups[[k]] <- list(seed=save.seed,
                                                             i=i,type=0)
			else screw.ups[[k]] <- list(seed=save.seed,i=i,type=1)
			next
                }
		tmp <- mixreg(x,yboot,ncomp=ncomp+ncincr,theta.start=theta2,
                              intercept=intercept,verb=verb,...)
		l2  <- tmp$log.like
		a2  <- tmp$aic
		if(!tmp$converged) {
			k <- k+1
			if(is.na(l2)) screw.ups[[k]] <- list(seed=save.seed,
                                                             i=i,type=2)
                        else screw.ups[[k]] <- list(seed=save.seed,i=i,type=3)
			next
                }
		if(l1 > l2) {
			k <- k+1
                        screw.ups[[k]] <- list(seed=save.seed,i=i,type=4)
                        next
		} else break
	}
	rslt[[i]] <- 2*(l2-l1)
	aic1[[i]] <- a1
	aic2[[i]] <- a2
	if(print.prog) cat(i,'\n')
}

rslt <- sort(unlist(rslt))
lrs  <- 2*(obj2$log.like-obj1$log.like)
pval <- sum(lrs<=rslt)/nboot
if(length(screw.ups) > 0) {
	seeds <- matrix(unlist(lapply(screw.ups,function(x){x[[1]]})),
                        byrow=TRUE,nrow=length(screw.ups))
	times <- unlist(lapply(screw.ups,function(x){x[[2]]}))
	types <- unlist(lapply(screw.ups,function(x){x[[3]]}))
	scrps <- list(seeds=seeds,times=times,types=types)
}
else scrps <- NULL
df <- ncincr*length(theta1[[1]])
xxx <- list(lrs=lrs,pval.boot=pval,lrs.boot=rslt,unlist(aic1),unlist(aic2),
            screw.ups=scrps,df=df)
names(xxx)[4] <- paste('aic',ncomp,sep='.')
names(xxx)[5] <- paste('aic',ncomp+ncincr,sep='.')
xxx
}
