genWbMixIni <- function()
{
	p <- runif(1, 0.1, 0.9)
	sh <- runif(2, 5, 20)
	sc <- runif(2, 4, 10)
	return(c(p, 1-p, sh, sc))
}


simWbMix<- function(n, mixParm)
{
	# function to produce the mixture data
	flagVec <- rbinom(n, size=1, prob=mixParm[1,1])
	pswd <- rep(NA, n)
	n1 <- sum(flagVec)
	n2 <- n-n1
	
	if (n1>0)
		pswd[flagVec==1] <- rweibull(n1, mixParm[1,2], mixParm[1,3])
	if (n2>0)
		pswd[flagVec==0] <- rweibull(n2, mixParm[2,2], mixParm[2,3])
	return(pswd)
}

quanWbMix <- function(intProb, mixParm)
{
	#mixmx, matrix return from findmx, first column the prob, shape, scale
	if (nrow(mixParm)!=2||ncol(mixParm)!=3)
		stop("Wrong Input of the Parameters")
	
	nq <- length(intProb)
	quan.est <- rep(NA, nq)
	t.fun <- function(x, tprob, mixParm)
	{
		p1 <- mixParm[1,1]*pweibull(x, shape=mixParm[1,2], scale= mixParm[1,3])
		p2 <- mixParm[2,1]*pweibull(x, shape=mixParm[2,2], scale= mixParm[2,3])
		return(p1+p2-tprob)
	}
	for (i in 1:nq)
	{
		quan.est[i] <- uniroot(t.fun, interval=c(0, 15), tprob=intProb[i], mixParm=mixParm)$root
	}
	return(rbind(intProb, quan.est))
}

quanWbMix.int <- function(intProb, mixmx)
{
	#Internal version, with only one intProb.
	if (length(mixmx)!=6||length(intProb)!=1)
		stop("Wrong Input")
	mixmx <- matrix(mixmx, ncol=3)
	t.fun <- function(x, mixmx)
	{
		p1 <- mixmx[1,1]*pweibull(x, shape=mixmx[1,2], scale= mixmx[1,3])
		p2 <- mixmx[2,1]*pweibull(x, shape=mixmx[2,2], scale= mixmx[2,3])
		return(p1+p2-intProb)
	}
	
	quan.est <- uniroot(t.fun, interval=c(0, 15), mixmx=mixmx)$root
	return(quan.est)
}

# The Weibull mixture with censored Likelihood.
emCenWbMix.T1 <- function(dat, Cx=NULL, iniParam=NULL, useC=FALSE, conCr=1e-6, nIter=10000)
{
	if (is.null(Cx))
		Cx <- max(dat)+ 10
	n <- length(dat)
	dat <- dat[!is.na(dat)]
	wdt <- dat[dat<=Cx]
	if(any(wdt<=0))
		stop("Positive observations required!")
		
	m <- length(wdt)
	if (m<=8)
		stop("Too few observations")
	if (is.null(iniParam))
		iniParam <- genWbMixIni()
	else
	{
		if (length(iniParam)!=6||any(iniParam<=0))
			stop("Wrong input of the initial parameters")
		else
		{
			if (abs(sum(iniParam[1:2])-1)>1E-6)
				stop("Wrong input of the mixing proportions")
		}
	}
	tmix <- matrix(iniParam, ncol=3)
	if (useC)
	{
		tmle <-.Call("R2C_CWbMix", wdt, Cx, n, m, iniParam, conCr, nIter)
	}
	else
	{
		tmle <- emCenWeibullMixture(wdt, Cx, n, m, tmix, conCr, nIter)
	}
	#return(tmle)
	##The second is negative log-likelihood
	mx <- matrix(tmle[-c(1,2)], ncol=3)
	rownames(mx) <- c("Population 1", "Population 2")
	colnames(mx) <- c("Proportion", "Shape", "Scale")
	return(list(convergence=tmle[1], nllh=tmle[2], estimates=mx, iniParam=iniParam))
}

emCenWbMix.T2 <- function(dat, n, iniParam=NULL, useC=FALSE, conCr=1e-6, nIter=10000)
{
	if (any(is.na(dat)))
		stop("No NA allowed in the input!")
	if(any(dat<=0))
		stop("Positive observations required!")
	Cx<- max(dat)	
	m <- length(dat)
	if (m<=8)
		stop("Too few observations")
	if (is.null(iniParam))
		iniParam <- genWbMixIni()
	else
	{
		if (length(iniParam)!=6||any(iniParam<=0))
			stop("Wrong input of the initial parameters")
		else
		{
			if (abs(sum(iniParam[1:2])-1)>1E-6)
				stop("Wrong input of the mixing proportions")
		}
	}
	tmix <- matrix(iniParam, ncol=3)
	if (useC)
	{
		tmle <-.Call("R2C_CWbMix", dat, Cx, n, m, iniParam, conCr, nIter)
	}
	else
	{
		tmle <- emCenWeibullMixture(dat, Cx, n, m, tmix, conCr, nIter)
	}
	#return(tmle)
	mx <- matrix(tmle[-c(1,2)], ncol=3)
	rownames(mx) <- c("Population 1", "Population 2")
	colnames(mx) <- c("Proportion", "Shape", "Scale")
	return(list(convergence=tmle[1], nllh=tmle[2], estimates=mx, iniParam=iniParam))
}

emCenWeibullMixture <-function(wdt, Cx, n, m, tmix, conCr, nIter)
{
	em.ncou <- 0
	ll.flag <- 1000
	
	t.llh <- log.llh.p2.censor(wdt, Cx, n, m, tmix)
	# EM.algorithm here
	while(em.ncou <nIter & ll.flag > conCr)
	{
		# E step
		E.mx <- expmx.p2.censor(wdt, Cx, n, m, tmix)
		
		# M stp
		if (any(is.na(E.mx)))
		{
			return(c(5, -t.llh, rep(NA, 6)))
		}
		else
		{
			omix <- tmix
			
			tmix[1,1]<- mean(E.mx[,1])
			tmix[2,1] <- mean(E.mx[,2])

			tmvec1 <- p2max.weibull.censor(wdt, Cx, n, m, E.mx[,1], tmix[1,2], conCr*1E-4, nIter/10)
			tmvec2 <- p2max.weibull.censor(wdt, Cx, n, m, E.mx[,2], tmix[2,2], conCr*1E-4, nIter/10)
			if (tmvec1[1]!=0 | tmvec2[1]!=0)
			{
				return(c(max(tmvec1[1], tmvec2[1]), t.llh, rep(NA, 6)))
			}
			else
			{
				tmix[1,2:3] <- tmvec1[-1]
				tmix[2,2:3] <- tmvec2[-1]
				o.llh <- t.llh
				t.llh <- log.llh.p2.censor(wdt, Cx, n, m, tmix)
				if (t.llh==-Inf | is.na(t.llh))
				{
					return(c(5, rep(NA, 7)))
				}
				else
				{
					ll.flag <- abs(t.llh - o.llh)/abs(o.llh)
				}
				em.ncou <- em.ncou + 1
			}
		}
	}
	if (em.ncou<nIter)
	{
		return(c(0, -t.llh, tmix))
	}
	else
	{
		return(c(6, -t.llh, tmix))
	}
}

log.llh.p2.censor <- function(wdt, Cx, n, m, mx)
{
	###This is log-likelihood
	lsum1 <- sum(log(mx[1,1]*dweibull(wdt, mx[1,2], mx[1,3]) + mx[2,1] * dweibull(wdt, mx[2,2], mx[2,3])))
	lsum2 <-0
	if ((n-m)>0)
		lsum2 <- (n-m) * (log( mx[1,1]*pweibull(Cx, mx[1,2], mx[1,3], lower.tail=F) + mx[2,1] * pweibull(Cx, mx[2,2], mx[2,3], lower.tail=F)))
	return(lsum1 + lsum2)
}

expmx.p2.censor  <-function(wdt, Cx, n, m, mx)
{
	E.mx <- matrix(NA, n, 2)

	tvec1 <- mx[1, 1] * dweibull(wdt, mx[1,2], mx[1,3]) 
	tvec2 <- mx[2, 1] * dweibull(wdt, mx[2,2], mx[2,3])
	E.mx[1:m, 1] <-  tvec1
	E.mx[1:m, 2] <-  tvec2
	if(n>m)
	{
		c.vec <- c(NA, NA)
		for (j in 1:2)
		{
			c.vec[j] <- mx[j,1] * pweibull(Cx, mx[j,2], mx[j, 3], lower.tail=F) 
		}
		E.mx[-c(1:m), 1] <- c.vec[1]
		E.mx[-c(1:m), 2] <- c.vec[2]
	}
	Esum <- (E.mx[,1] + E.mx[,2])
	E.mx <- E.mx/Esum	
	return(E.mx)
}

p2max.weibull.censor<-function(wdt, Cx, n, m, E.x, t.al, conCr, nIter)
{
	flag=1000
	ncount=0
	
	ncE.x <- E.x[1:m]
	cwt <-1
	if(n>m)
		cwt <- E.x[(m+1)]
	
	while(flag>conCr & ncount<nIter)
	{
		ncount <- ncount+1
		np1 <- sum(ncE.x * (wdt^t.al) * log(wdt))
		dnp1<- sum(ncE.x * (wdt^t.al))
		
		np3 <- sum(ncE.x * log(wdt))		
		dnp3 <- sum(ncE.x)
		
		np2 <-  (n-m) * cwt * Cx^t.al * log(Cx)
		dnp2 <- (n-m) * cwt * Cx^t.al
		
		inv.al <- (np1 + np2)/(dnp1 + dnp2) - np3/dnp3
			
		if (!is.na(inv.al))
		{
			flag <- abs(1/inv.al - t.al)
			t.al <- 1/inv.al
			if(t.al <0)
			{
				return(c(3, NA, NA))
			}
		}
		else
		{
			return(c(1, NA, NA))
		}
	}
	if (ncount>=nIter)
	{
		return(c(2, NA, NA)) 
	}
	else
	{
		t.le <- sum(ncE.x) / (dnp1 + dnp2)
		t.le <- 1/(t.le^(1/t.al))
		return(c(0, t.al, t.le)) 
	}
}


