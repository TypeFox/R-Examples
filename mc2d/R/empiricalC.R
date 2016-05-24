#<<BEGIN>>
dempiricalC <- function(x, min, max, values, prob=NULL, log=FALSE)
#TITLE The Continuous Empirical Distribution
#NAME empiricalC
#KEYWORDS distribution
#DESCRIPTION
#Density, distribution function and random generation for a continuous empirical distribution.
#INPUTS
#{x, q}<<Vector of quantiles.>>
#{p}<<Vector of probabilities.>>
#{n}<<Number of random values. If \samp{length(n) > 1}, the length is taken to be the number required.>>
#{min}<<A finite minimal value.>>
#{max}<<A finite maximal value.>>
#{values}<<Vector of numerical values.>>
#[INPUTS]
#{prob}<<Optionnal vector of count or probabilities.>>
#{log, log.p}<<logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.>>
#{lower.tail}<<logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.>>
#DETAILS
#Given \eqn{p_{i}}{p_i}, the distribution value for \eqn{x_{i}}{x_i} with \samp{i} the rank \eqn{i = 0, 1, 2, \ldots, N+1},
# \eqn{x_{0}=min}{x_0 = min} and \eqn{x_{N+1}=max}{x_(N+1) = max} the density is:
# \deqn{f(x)=p_{i}+(\frac{x-x_{i}}{x_{i+1}-x_{i}})(p_{i+1}-p_{i})}{f(x) = p_i + (p_(i+1) - p_i)/(x_(i+1) - x_i) for x_i<=x<x_(i+1)}
# The \samp{p} values being normalized to give the distribution a unit area.
#
#\samp{min} and/or \samp{max} and/or \samp{values} and/or \samp{prob} may vary: in that case, 
#\samp{min} and/or \samp{max} should be vector(s). \samp{values} and/or \samp{prob} should be 
#matrixes, the first row being used for the first element of \samp{x}, \samp{q}, \samp{p} or the first random value, the
#second row for the second element of \samp{x}, \samp{q}, \samp{p} or random value, ...
#Recycling is permitted if the number of elements of \samp{min} or \samp{max} or the number of rows of \samp{prob} and \samp{values} are equal or
#equals one.</>
#SEE ALSO
#\code{\link{empiricalD}}
#VALUE
#\samp{dempiricalC} gives the density, \samp{pempiricalC} gives the distribution function,
#\samp{qempiricalC} gives the quantile function and \samp{rempiricalC} generates random deviates.
#EXAMPLE
#prob <- c(2,3,1,6,1)
#values <- 1:5
#par(mfrow=c(1,2))
#curve(dempiricalC(x, min=0, max=6, values, prob), from=-1, to=7, n=1001)
#curve(pempiricalC(x, min=0, max=6, values, prob), from=-1, to=7, n=1001)
#
### Varying values
#(values <- matrix(1:10,ncol=5))
### the first x apply to the first row 
### the second x to the second one
#dempiricalC(c(1,1),values,min=0, max=11)
#
#
###Use with mc2d 
#val <- c(100, 150, 170, 200)
#pr <- c(6, 12, 6, 6)
#out <- c("min", "mean", "max")
###First Bootstrap in the uncertainty dimension
###with rempirical D
#(x <- mcstoc(rempiricalD, type = "U", outm = out, nvariates = 30, values = val, prob = pr))
###Continuous Empirical distribution in the variability dimension
#mcstoc(rempiricalC, type = "VU", values = x, min=90, max=210)

#CREATED 08-02-20
#--------------------------------------------
{
  
	lx   <- length(x)
	if(lx == 0) return(numeric(0))
    
	if(is.vector(values)) values <- matrix(values,nrow=1)
	if(is.null(prob))     prob   <- rep(1,length(values[1,]))
	if(is.vector(prob))   prob   <- matrix(prob,  nrow=1)

	lmin <- length(min)
	lmax <- length(max)
    nrv  <-  nrow(values)
    nrp  <-  nrow(prob)
	lpar <- max(lmin,lmax,nrv,nrp) # nb of combination of parameters
	if(nrv != 1 && nrp != 1 && nrv != nrp)
		stop("values/prob should be vector(s), matrix(es) of 1 row or matrix(es) of the same number of rows")

	if(any(!is.finite(values)) || any(!is.finite(min)) || any(!is.finite(max)) || any(!is.finite(prob))) stop("values, prob, min and ax should be finite values") 
	if(any(min > max) || any(min > apply(values,1,min)) || max < apply(values,1,max)) stop("at least one min is not a minimum or max is not a maximum")
	if(any(prob < 0) || any(apply(prob,1,sum) == 0)) stop("Prob should be non negative and sum(prob) should be != 0")

	onerow <- function(x, min, max, values, prob){

		val2 <- sort(unique(values))
				
		probi <- tapply(prob,values,sum)
		probi <- probi / sum(probi)
		
		prob <- c(0,0,probi,0,0)
		val <- c(-Inf,min,val2,max,Inf)
  
		h <- c(val2,max) - c(min,val2)
		a <- c(0,probi)
		b <- c(probi,0)
		Integ <- sum(h*(a+b)/2)

		d <- rep(NA,length(x))
		d[x >= max | x <= min] <- 0
		lesquel <- which(!is.na(x) & x < max & x > min)
		quel <- findInterval(x[lesquel], val) + 1
		d[lesquel] <- prob[quel-1]+(x[lesquel]-val[quel-1])/(val[quel]-val[quel-1])*(prob[quel]-prob[quel-1])
		d <- d / Integ
		return(d)
	}
	
	if(lpar == 1) { #Shortcut if only one set of parameters
		d <- onerow(x, as.vector(min), as.vector(max), as.vector(values), as.vector(prob))} else 
		{ # launch lpar time the function	
		x      <- matrix(x, nrow = lpar)
		x      <- lapply(1:lpar, function(y) x[y,])
		prob   <- lapply(1:nrp,  function(x) prob[x,])
		values <- lapply(1:nrv,  function(x) values[x,])
		d <- as.vector(t(mapply(onerow, x, min, max, values, prob)))
		d <- d[1:max(lx,lpar)]
		}
	
	if(log) d <- log(d)
	if(any(is.na(d))) warning("NaN in dempiricalC")
  return(d)
  }

#<<BEGIN>>
pempiricalC <- function(q, min, max, values, prob=NULL, lower.tail = TRUE, log.p = FALSE)                 
#ISALIAS dempiricalC
#--------------------------------------------
{

   lq <- length(q)
   if(lq == 0) return(numeric(0))
   
	if(is.vector(values)) values <- matrix(values,nrow=1)
	if(is.null(prob))     prob   <- rep(1,length(values[1,]))
	if(is.vector(prob))   prob   <- matrix(prob,  nrow=1)

	lmin <- length(min)
	lmax <- length(max)
    nrv  <-  nrow(values)
    nrp  <-  nrow(prob)
	lpar <- max(lmin,lmax,nrv,nrp) # nb of combination of parameters
	if(nrv != 1 && nrp != 1 && nrv != nrp)
		stop("values/prob should be vector(s), matrix(es) of 1 row or matrix(es) of the same number of rows")

	if(any(!is.finite(values)) || any(!is.finite(min)) || any(!is.finite(max)) || any(!is.finite(prob))) stop("values, prob, min and ax should be finite values") 
	if(any(min > max) || any(min > apply(values,1,min)) || max < apply(values,1,max)) stop("at least one min is not a minimum or max is not a maximum")
	if(any(prob < 0) || any(apply(prob,1,sum) == 0)) stop("Prob should be non negative and sum(prob) should be != 0")

  onerow <- function(q,min,max,values,prob){

	val2 <- sort(unique(values))

	probi <- tapply(prob,values,sum)
	h <- c(val2,max) - c(min,val2)
	a <- c(0,probi)
	b <- c(probi,0)
	Integ <- sum(h*(a+b)/2)
	probi <- probi / Integ

	a <- c(0,probi)
	b <- c(probi,0)
	probcum <- cumsum(h*(a+b)/2)

	probi <- c(0,0,probi,0,0)
	probcum <- c(0,0,probcum,1)
	val <- c(-Inf,min,val2,max,Inf)

	p <- rep(NA,length(q))
	p[q >= max] <- 1
	p[q <= min] <- 0
	lesquel <- which(!is.na(q) & q < max & q > min)
	quel <- findInterval(q[lesquel], val) + 1

	p[lesquel] <- probcum[quel-1]+(q[lesquel]-val[quel-1])*
					(probi[quel-1]+((probi[quel]-probi[quel-1])*(q[lesquel]-val[quel-1])/(2*(val[quel]-val[quel-1]))))
	return(p)
	}
	
	if(lpar == 1) { #Shortcut if only one set of parameters
		p <- onerow(q, as.vector(min), as.vector(max), as.vector(values), as.vector(prob))} else 
		{ # launch lpar time the function	
		q      <- matrix(q, nrow = lpar)
		q      <- lapply(1:lpar, function(y) q[y,])
		prob   <- lapply(1:nrp,  function(x) prob[x,])
		values <- lapply(1:nrv,  function(x) values[x,])
		p <- as.vector(t(mapply(onerow, q, min, max, values, prob)))
		p <- p[1:max(lq,lpar)]
		}

	if(!lower.tail) p <- 1-p
	if(log.p) p <- log(p)
	if(any(is.na(p))) warning("NaN in pempiricalC")
	
  return(p)}


#<<BEGIN>>
qempiricalC <- function(p, min, max, values, prob=NULL, lower.tail = TRUE, log.p = FALSE)
#ISALIAS dempiricalC
#--------------------------------------------
{
	lp <- length(p)
	if(lp == 0) return(numeric(0))
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1-p

	if(is.vector(values)) values <- matrix(values,nrow=1)
	if(is.null(prob))     prob   <- rep(1,length(values[1,]))
	if(is.vector(prob))   prob   <- matrix(prob,  nrow=1)

	lmin <- length(min)
	lmax <- length(max)
	nrv  <-  nrow(values)
	nrp  <-  nrow(prob)
	lpar <- max(lmin,lmax,nrv,nrp) # nb of combinations of parameters
	if(nrv != 1 && nrp != 1 && nrv != nrp)
		stop("values/prob should be vector(s), matrix(es) of 1 row or matrix(es) of the same number of rows")

	if(any(!is.finite(values)) || any(!is.finite(min)) || any(!is.finite(max)) || any(!is.finite(prob))) stop("values, prob, min and ax should be finite values") 
	if(any(min > max) || any(min > apply(values,1,min)) || max < apply(values,1,max)) stop("at least one min is not a minimum or max is not a maximum")
	if(any(prob < 0) || any(apply(prob,1,sum) == 0)) stop("Prob should be non negative and sum(prob) should be != 0")

	onerow <- function(p, min, max, values, prob){ 
		val2 <- sort(unique(values))
		probi <- tapply(prob,values,sum)
		h <- c(val2,max) - c(min,val2)
		a <- c(0,probi)
		b <- c(probi,0)
		Integ <- sum(h*(a+b)/2)
		probi <- probi / Integ

		a <- c(0,probi)
		b <- c(probi,0)
		probcum <- cumsum(h*(a+b)/2)

		probi <- c(0,0,probi,0,0)
		probcum <- c(0,0,probcum,Inf)
		val <- c(-Inf,min,val2,max,Inf)

		q <- rep(NA,length(p))
		q[p > 1] <- NaN
		q[p == 1] <- max
		q[p < 0] <- NaN

		lesquel <- which(!is.na(p) & p < 1 & p >= 0)

		quel <- findInterval(p[lesquel],probcum)+1

		a <- (probi[quel]-probi[quel-1])/(val[quel]-val[quel-1])/2	# (c-a)/(2b)
		b <- probi[quel-1]											# a
		c <- probcum[quel-1]-p[lesquel]								# -S
		d <- b^2-4*a*c
		q[lesquel] <- ifelse(a == 0, -c/b + val[quel-1], (-b+sqrt(d))/2/a  + val[quel-1])
		return(q)
		}
	
	if(lpar == 1) { #Shortcut if only one set of parameters
		q <- onerow(p, as.vector(min), as.vector(max), as.vector(values), as.vector(prob))} else 
		{ # launch lpar time the function	
		p      <- matrix(p, nrow = lpar)
		p      <- lapply(1:lpar, function(y) p[y,])
		prob   <- lapply(1:nrp,  function(x) prob[x,])
		values <- lapply(1:nrv,  function(x) values[x,])
		q <- as.vector(t(mapply(onerow, p, min, max, values, prob)))
		q <- q[1:max(lp,lpar)]
		}
		
	if(any(is.na(q))) warning("NaN in qempiricalC")
	return(q)}

#<<BEGIN>>
rempiricalC <- function(n, min, max, values, prob=NULL)
#ISALIAS dempiricalC
#--------------------------------------------
{ 
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) cannot be negative in rempiricalc")
  if(length(min) > n) min <- min[1:n]
  if(length(max) > n) max <- max[1:n]
  if(is.vector(values) && length(values) > n) values <- values[1:n] 
  else if(is.matrix(values) && nrow(values) > n)   values <- values[1:n,]
  if(is.vector(prob) && length(prob) > n) prob <- prob[1:n]
  else if(is.matrix(prob) && nrow(prob) > n)   prob <- prob[1:n,]
  
  r <- qempiricalC(runif(n), min=min, max=max, values=values, prob=prob, lower.tail = TRUE, log.p = FALSE)
  return(as.vector(r))

  }
