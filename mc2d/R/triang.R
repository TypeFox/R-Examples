#<<BEGIN>>
dtriang <- function(x, min=-1, mode=0, max=1, log=FALSE)
#TITLE The Triangular Distribution
#NAME triangular
#KEYWORDS distribution
#DESCRIPTION
#Density, distribution function, quantile function and random generation
#for the triangular distribution with minimum equal to \samp{min}, mode equal \samp{mode}
#and maximum equal to \samp{max}.
#INPUTS
#{x,q}<<vector of quantiles.>>
#{p}<<vector of probabilities.>>
#{n}<<number of observations. If length(n) > 1, the length is taken to be the number required.>>
#[INPUTS]
#{min}<<vector of minima.>>
#{mode}<<vector of modes.>>
#{max}<<vector of maxima.>>
#{log, log.p}<<logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.>>
#{lower.tail}<<logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.>>
#DETAILS
#For the case of u := min == mode == max, there is no density in that case and dtriang will return NaN (the error condition) (Similarity with dunif).
#VALUE
#\samp{dtriang} gives the density, \samp{ptriang} gives the distribution function,
#\samp{qtriang} gives the quantile function, and \samp{rtriang} generates random deviates.

#EXAMPLE
#curve(dtriang(x, min=3, mode=5, max=10), from = 2, to = 11)
###no density when  min == mode == max
#dtriang(c(1,2,3),min=2,mode=2,max=2)
#CREATED 08-02-20
#--------------------------------------------
{
	if(length(x) == 0) return(numeric(0))
	min <- as.vector(min)
	mode <- as.vector(mode)
	max <- as.vector(max)
	
	# quel: x < mode or x = mode = max 
	xmaxmode <- (abs(x-max) < (.Machine$double.eps^0.5)) & (abs(max-mode) < (.Machine$double.eps^0.5)) 
	quel <- (x < mode) | xmaxmode  
	d <- ifelse(quel,
              2*(x-min)/((mode-min)*(max-min)),
	            2 *(max-x)/((max-mode)*(max-min)))

	d[x < min | x > max] <- 0
	d[mode < min | max < mode] <- NaN

	# For min = mode = max: provide an error like in dunif
	xminmodemax <- (abs(min-max)) < (.Machine$double.eps^0.5)
	d[xminmodemax] <- NaN
	
	if(log) d <- log(d)
	if(any(is.na(d))) warning("NaN in dtriang")
  return(d)}

#<<BEGIN>>
ptriang <- function(q,min=-1,mode=0,max=1,lower.tail = TRUE, log.p = FALSE)
#ISALIAS dtriang
#--------------------------------------------
{
	if(length(q) == 0) return(numeric(0))
	# quel: q < mode or q = mode = max 
	min <- as.vector(min)
	mode <- as.vector(mode)
	max <- as.vector(max)
	qmaxmode <- (abs(q-max) < (.Machine$double.eps^0.5)) & (abs(max-mode) < (.Machine$double.eps^0.5)) 
	quel <- (q < mode) | qmaxmode  
	p <- ifelse(quel,
              (q-min)^2 / ((mode-min)*(max-min)),
	             1 - ((max-q)^2/((max-mode)*(max-min))))
	#if q = max = mode = min
	qminmodemax <- qmaxmode & (abs(q - min) < .Machine$double.eps^0.5)
	p[qminmodemax] <- 1

	p[q < min] <- 0
	p[q > max] <- 1
	p[mode < min | max < mode] <- NaN
    if(!lower.tail) p <- 1-p
    if(log.p) p <- log(p)
	if(any(is.na(p))) warning("NaN in ptriang")
  return(p)}

#<<BEGIN>>
qtriang <- function(p, min=-1, mode=0, max=1, lower.tail=TRUE, log.p=FALSE)
#ISALIAS dtriang
#--------------------------------------------
{
	if(length(p) == 0) return(numeric(0))
	min <- as.vector(min)
	mode <- as.vector(mode)
	max <- as.vector(max)
    if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1-p
	quel <- p <= (mode-min)/(max-min)
	q <- ifelse(quel,
              min + sqrt(p*(mode-min)*(max-min)),
              max - sqrt((1-p)*(max-min)*(max-mode)))
	
	#if max = min (and then = mode)
	minmodemax <- (abs(min-max) < (.Machine$double.eps^0.5)) 
	q <- ifelse(rep(minmodemax, length.out=length(p)), min, q)
	
	q[p < 0 | p > 1] <- NaN
	q[mode < min | max < mode] <- NaN
	if(any(is.na(q))) warning("NaN in qtriang")
  return(q)}


#<<BEGIN>>
rtriang <- function(n, min=-1, mode=0, max=1)
#ISALIAS dtriang
#--------------------------------------------
{ if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) can not be negative in rtriang")
  U <- runif(n)
  ow <- options(warn=-1)
  q <- qtriang(U,	min=min,
					mode=mode,
					max=max, lower.tail = TRUE, log.p = FALSE)
  options(ow)
  if(any(is.na(q))) warning("NaN in rtriang")
   
  return(q)}
