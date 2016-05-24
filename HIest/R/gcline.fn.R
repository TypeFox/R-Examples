gcline.fn <- function(x,n,y,start=NULL,model="logit-logit",method="L-BFGS-B",iterations=99,SD=rep(0.01,length(start)),headstart=FALSE,Grid=TRUE){
	# requires nnet and stats
	require(stats)
	require(nnet)
	# this function is for fitting models, not identifying outliers
	if(!is.null(start)) start <- as.numeric(start)
	X.out <- is.na(y)
	if(sum(X.out)>0){
		X.in <- !is.na(y)
		x <- x[X.in]
		n <- n[X.in]
		y <- y[X.in]
		cat("\n",sum(X.out)," NA's were omitted\n")
		}
	if(model=="logistic"){
		method <- "glm.fit"
		Y <- cbind(y,n-y)
		fit.glm <- glm(Y~x,family="binomial")
		LL.fit <- logLik(fit.glm)
		k <- attr(LL.fit,"df")
		est <- list(par=fit.glm$coef,convergence=fit.glm$conv,value=LL.fit[1])
		null.LL <- NA
	}
	if(model=="multinom"){
		if(!is.matrix(y)){stop("error, multinomial regression requires a matrix of genotype counts")}
		method <- "nnet"
		fit.mlt <- multinom(y~x,trace=FALSE)
		LL.fit <- logLik(fit.mlt)
		k <- attr(LL.fit,"df")
		est <- list(par=summary(fit.mlt)$coeff,convergence=fit.mlt$conv,value=LL.fit[1])
		null.LL <- NA
	}
	if(model=="logit-logit"){
		if(is.null(start)){ start <- c(u=0,v=1)}
		k <- 2
		SD <- rep(0.1,2)
		LL.fn <- function(par,SD=SD,x,n,y){
			u <- par[1]
			v <- par[2]
			p <- x^v/(x^v + (1-x)^v*exp(u))
			p <- replace(p,x==1,1)
			p <- replace(p,x==0,0)
			p <- replace(p,p<=0,.Machine$double.xmin)
			p <- replace(p,p>=1,1-.Machine$double.neg.eps)
			sum(dbinom(y,n,p,log=TRUE))			
		}
		lower <- c(-Inf,0)
		upper <- c(Inf,Inf)
		GR <- function(par,SD=SD,x,n,y){
			newv <- rlnorm(1,par[1],SD[2])
			newu <- rnorm(1,par[2],SD[1])
			c(newu,newv)
			}
		null.LL <- LL.fn(par=c(u=0,v=1),SD,x,n,y) # straight line from (0,0) to (1,1)
	}
	if(model=="Beta"){
		if(is.null(start)) start <- c(mu=0.5,nu=2)
		SD <- rep(0.1,2)
		k <- 2
		LL.fn <- function(par,SD=SD,x,n,y){
			p <- pbeta(x,max(par[1]*par[2],.Machine$double.xmin),max((1-par[1])*par[2],.Machine$double.xmin))
			p <- replace(p,p<=0,.Machine$double.xmin)
			p <- replace(p,p>=1,1-.Machine$double.neg.eps)
			# dens <- dbinom(y,n,p,log=TRUE)
			# dens <- replace(dens,dens<(-.Machine$double.xmax),-.Machine$double.xmax)
			# sum(dens)
			sum(dbinom(y,n,p,log=TRUE))
			}
		lower <- c(.Machine$double.eps,.Machine$double.eps)
		upper <- c(1-.Machine$double.neg.eps,Inf)
		GR <- function(par,SD=SD,x,n,y){
			newa <- max(.Machine$double.xmin,rnorm(1,par[1],SD[1]))
			newa <- min(newa,1-.Machine$double.neg.eps)
			newb <- max(.Machine$double.xmin,rnorm(1,par[2],SD[2]))
			c(newa,newb)
			}
		null.LL <- LL.fn(par=c(mu=.5,nu=2),SD,x,n,y) # straight line from (0,0) to (1,1)
		}
	if(model=="Richards"){
		if(is.null(start)) start <- c(p1=1e+5,p2=1-1e+5,b=-2*log((1e+5-1)/1e+5),m=1/2)
		k <- 4
		SD <- rep(0.1,4)
		LL.fn <- function(par,SD=SD,x,n,y){
			p <- par[1]+(par[2]-par[1])/(1+exp(par[3]*(x-par[4])))
			p <- replace(p,p<=0,.Machine$double.xmin)
			p <- replace(p,p>=1,1-.Machine$double.neg.eps)
			sum(dbinom(y,n,p,log=TRUE))
			}
		lower <- c(-Inf,-Inf,-Inf,-Inf)
		upper <- c(Inf,Inf,Inf,Inf)
		GR <- function(par,SD=SD,x,n,y){
			p1 <- rnorm(1,par[1],SD[1])
			p2 <- rnorm(1,par[2],SD[2])
			b <- rnorm(1,par[3],SD[3])
			m <- rnorm(1,par[4],SD[4])
			c(p1,p2,b,m)
			}
			p1 <- 1e+5
		null.LL <- LL.fn(par=c(p1,1-p1,-2*log((p1-1)/p1),1/2),SD=SD,x,n,y) # a pretty straight line between (0,0) and (1,1)
		}
	if(model=="Barton"){
		if(is.null(start)) start <- c(a=0,b=0)
		SD <- rep(0.1,2)
		k <- 2
		LL.fn <- function(par,SD=SD,x,n,y){
			p <- x+2*x*(1-x)*(par[1]+par[2]*(2*x-1))
			p <- replace(p,p<=0,.Machine$double.xmin)
			p <- replace(p,p>=1,1-.Machine$double.neg.eps)
			sum(dbinom(y,n,p,log=TRUE))
			}
		lower <- c(-Inf,-Inf)
		upper <- c(Inf,Inf)
		GR <- function(par,SD=SD,x,n,y){
			newa <- rnorm(1,par[1],SD[1])
			newb <- rnorm(1,par[2],SD[2])
			c(newa,newb)
			}
		null.LL <- LL.fn(par=c(a=0,b=0),SD,x,n,y) # straight line from (0,0) to (1,1)
		}
	if(method=="L-BFGS-B"){
		est <- optim(par=start,fn=LL.fn,x=x,n=n,y=y,method=method,lower=lower,upper=upper,control=list(fnscale=-1))
		}
	if(method=="SANN"){
		if(headstart){
			start <- optim(par=start,fn=LL.fn,x=x,n=n,y=y,method="L-BFGS-B",lower=lower,upper=upper,control=list(fnscale=-1))$par
			}
		est <- optim(par=start,fn=LL.fn,gr=GR,SD=SD,x=x,n=n,y=y,method="SANN",control=list(fnscale=-1))
		}
	if(method=="mcmc"){
		if(headstart){
			start <- optim(par=start,fn=LL.fn,x=x,n=n,y=y,method="L-BFGS-B",lower=lower,upper=upper,control=list(fnscale=-1))$par
			}
		if(Grid & model=="Beta"){
			mu <- seq(from=0.02,to=0.90,length.out=10)
			nu <- 2^(0:9)/10
			LLG <- data.frame(mu=rep(mu,10),nu=rep(nu,each=10),LLik=NA)
			for(i in 1:100){
				LLG$LLik[i] <- LL.fn(as.numeric(LLG[i,1:2]),x=x,n=n,y=y)
			}
			start <- as.numeric(LLG[which.max(LLG$LLik),1:2])
		}
		chain <- matrix(nrow=iterations+1,ncol=length(start))
		chain[1,] <- start
		LLik <- LL.fn(par=start,x=x,n=n,y=y)
		for(i in 1:iterations){
			newpar <- GR(chain[i,],SD,x,n,y)
			newLL <- LL.fn(par=newpar,x=x,n=n,y=y)
			if(runif(1)<exp(newLL-LLik[i])){
				chain[i+1,] <- newpar
				LLik[i+1] <- newLL
				}else{
					chain[i+1,] <- chain[i,]
					LLik[i+1] <- LLik[i]
					}
			}
		colnames(chain) <- names(start)
		est <- list(par=chain[which.max(LLik),],value=max(LLik),convergence = data.frame(chain,LLik))
		}
	## estimation is done, now report
	estimates <- est$par
	convergence <- est$convergence
	names(estimates) <- names(start)
	N <- sum(n)
	lnL <- est$value
	AICc <- 2*k-2*lnL+2*k*(k+1)/(N-k-1)
	list(model=model,method=method,estimates=estimates,lnL=c(fitted=lnL,null=null.LL),k=k,AICc=AICc,convergence=convergence)
	}
