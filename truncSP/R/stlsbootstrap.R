funcSTLS <- function(mf,d,formula,beta,bet,point,direction,control)
{
	mf2 <- mf[d,]
	x <- model.matrix(formula,data=mf2)
	
	y <- matrix(mf2[,1])

	if(is.numeric(beta))
	{
		bet <- bet
	}
	else
	{
		if(beta=="ols")
		{
			bet <- lm(data=mf2)$coef
			bet <- matrix(bet)
		}

		if(beta=="ml")
		{
			mlcof<-mlcoef(mf2)
			bet <- matrix(mlcof[-length(mlcof)])
		}
	}

	
	opt <- optim(par=bet,fn=funcval.STLS,control=control,x=x,y=y)

	if(direction=="right")
	{
		opt$par <- -opt$par
	}

	if(point!=0)
	{
		opt$par[1,1] <- opt$par[1,1]+(point)
	}
	
	b <- opt$par

	return(b)
}


covar.bootSTLS <- function(mf,funcSTLS,R,formula,beta,bet,point,direction,control)
{
	bootrepl <- boot(data=mf, statistic=funcSTLS, R=R, formula=formula, beta=beta,bet=bet,point=point,direction=direction,control=control)$t	
	 
	bootmean <- apply(bootrepl,2,mean)

	bootmean <- matrix(bootmean)

	k <- ncol(bootrepl)

	M <- matrix(0,k,k)

	for(i in 1:R)
	{
		m <- bootrepl[i,]-bootmean
		M <- M + m%*%t(m)
	}


	covarboot <- (1/R)*M
	bootvals <- list(covarboot,bootrepl)
	return(bootvals)
}

bootconfintSTLS <- function(bootrepl,level)
{
	k <- ncol(bootrepl)
	bootconfint <- matrix(0,k,2)
	for(i in 1:k)
	{
		bootconfint[i,1] <- quantile(bootrepl[,i], (1-level)/2) 			
		bootconfint[i,2] <- quantile(bootrepl[,i], 1-(1-level)/2)			
	}
return(bootconfint)
}
	
