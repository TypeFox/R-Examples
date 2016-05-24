# Asymptotic maximum likelihood estimator of the shape parameter
`.aml` <-
function(x,n,k)
{
x <- sort(x)
nk <- n-k
x1 <- x[(nk+1):n]
return(-(log(x1[1])-sum(log(x1))/k))
}

`.amle` <-
function(x,n,k)
{
x <- sort(x)
nk <- n-k
x1 <- x[(nk+1):n]
g <- -(log(x1[1])-sum(log(x1))/k)
sigma <- g*exp(log(x1[1]) + g*log(k/n))
return(c(g,sigma))
}


# Combined estimator of the shape parameter
`.combined` <-
function(x)
{
m <- mean(x)
maxi <- max(x)
g <- m/(m-maxi)
sigma <- -g*maxi
return(c(g,sigma))
}

# Test statistic for H_0^-
`.r1` <-
function(x,n)
{
x <- sort(x)
gamman <- .combined(x)[1]
Fn <- seq(1,n,1)/n
Finv <- (1-Fn)^(-gamman)-1
return(abs(cor(x[1:(n-1)],Finv[1:(n-1)])))
}


# Approximated null distribution of the test statistic for H_0^-
`.cc1` <-
function(n,gamman,J)
{
Fn <- seq(1,n,1)/n
r <- rep(NA,J)
for(i in 1:J)
{
r[i] <- .r1(rgp(n, shape = gamman, scale = 1),n)
}
return(r)
}



# Test statistic for H_0^+
`.r2` <-
function(x,n,Fn)
{
	x <- sort(x)
	gammap <- .aml(x,n,ceiling(.2*n))
	Finv <- (1-Fn)^(-gammap)-1
 	if(gammap <= .5)
	{
		return(abs(cor(x[1:(n-1)],Finv[1:(n-1)])))
	}
	if(gammap > .5)
	{
		return(abs(cor(log(x[1:(n-1)]),log(Finv[1:(n-1)]))))
	}
}
 

# Approximated null distribution of the test statistic for H_0^+
`.cc2` <-
function(n,gammap,Fn, J)
{
	r <- rep(NA,J)
	for( i in 1:J)
	{
		xb <- rgp(n, shape = gammap,scale = 1)    #bootstrap sample
		r[i] <- .r2(xb,n,Fn)
	}
	return(r)
}

