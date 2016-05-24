##For Davor R version
cenWbMLE.T1 <- function(dat, Cx=NULL, useC=FALSE, conCr=1e-9, nIter=1000)
{
	if (is.null(Cx))
		Cx <- max(dat)+3
	else
		if (Cx <0)
			stop("Postive Cx required")
	
	n <- length(dat)
	dat <- dat[!is.na(dat)]
	wdt <- dat[dat<=Cx]
	if(any(wdt<=0))
		stop("Positive observations required!")
		
	m <- length(wdt)
	if (m<=3)
		stop("Too few observations")
	
	if (useC)
	{
		tmle <-.Call("R2C_cenWeibullMLE", wdt, Cx, n, m, conCr, nIter)
	}
	else
	{
		tmle <- cenWeibullMLE(wdt, Cx, n, m, conCr, nIter)
	}
	#return(tmle)
	return(list(convergence=tmle[1], estimates=c(Shape=tmle[2], Scale=tmle[3])))
}


cenWbMLE.T2 <-function(dat, n, useC=FALSE, conCr=1e-9, nIter=1000)
{
	if (any(is.na(dat)))
		stop("No NA allowed in the input!")
	if(any(dat<=0))
		stop("Positive observations required!")
	
	m <- length(dat)
	if (m<=3)
		stop("Too few observations")
		
	Cx<-max(dat)
	
	if (useC)
	{
		tmle <-.Call("R2C_cenWeibullMLE", dat, Cx, n, m, conCr, nIter)
	}
	else
	{
		tmle <- cenWeibullMLE(dat, Cx, n, m, conCr, nIter)
	}
	#return(tmle)
	return(list(convergence=tmle[1], estimates=c(Shape=tmle[2], Scale=tmle[3])))
}

cenWeibullMLE <- function(wdt, max.dt, n, m, conCr, nIter)
{
	t.al <- mean(wdt)/sd(wdt)
	flag=100 # judge the difference of 1/alpha and last alpha.
	ncount =0 # count the loop
	
	while(flag>conCr & ncount<nIter)
	{
		ncount <- ncount+1
		l.al <- t.al # t.al, temp alpha
		np1 <- sum(( wdt ^ t.al ) * log(wdt))
		np2 <- (n-m) * (max.dt^t.al) * log(max.dt)
		d1 <- sum(wdt ^ t.al) + (n-m) * (max.dt ^ t.al)
		id1 <- sum(log(wdt))/m
		if (any(is.na(c(np1, np2, d1, id1))))
		{
			return(c(1, NA, NA))
		}
		else
		{
			inv.alpha <- (np1 + np2)/d1 - id1
			if (is.finite(inv.alpha))
			{
				t.al <- 1/inv.alpha
				if (t.al <0)
				{
					return(c(3, NA, NA))
				}
				else
				{
					flag = abs( (1/l.al) - inv.alpha)
				}
			}
			else
			{
				return(c(4, NA, NA))
			}
		}
	}
	if (ncount<nIter)
	{
		e.le <- m/( sum(wdt ^ t.al) + (n-m) * (max.dt ^ t.al) )
		n.le<-1/(e.le^(1/t.al))
		return(c(0, t.al, n.le))
	}
	else
	{
		return(c(2, NA, NA))
	}
}