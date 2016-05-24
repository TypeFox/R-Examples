dweibull3 <- function(x,shape,scale,threshold)                                  #### WEIBULL3
{
 if(missing(x))
  stop("x must be a vector")
 if(missing(threshold))
  threshold=0
 if(missing(shape))
  shape=1
 if(missing(scale))
  scale=1
 temp=function(x)
{
	if(x>=threshold)
		return((shape/scale)*(((x-threshold)/scale)^(shape-1))*exp(-((x-threshold)/scale)^shape))
	else
		return(0)
}
 return(unlist(lapply(x,temp)))
}

pweibull3 <- function(q,shape,scale,threshold)
{
if(missing(q))
  stop("q must be a vector")
 if(missing(threshold))
  threshold=0
 if(missing(shape))
  shape=1
 if(missing(scale))
  scale=1
 temp=function(q)
 {
	if(q>=threshold)
		return(1-exp(-((q-threshold)/scale)^shape))
	else
		return(0)
 }
 return(unlist(lapply(q,temp)))
}

qweibull3 <- function(p,shape,scale,threshold,...)
{
 if(missing(p))
  stop("p must be a vector")
 if(missing(threshold))
  threshold=0
 if(missing(shape))
  shape=1
 if(missing(scale))
  scale=1
 myfun = function(x,p) pweibull3(q = x,
         threshold = threshold, scale = scale, shape = shape) - p
 temp=function(p)
 {
  return(uniroot(f=myfun,lower=threshold,upper=threshold+10000000,p=p,...)$root)       #solve myfun=0
 }
 return(unlist(lapply(p,temp)))
}

dlnorm3 <- function(x,meanlog,sdlog,threshold)                                  #### LOG-NORM3
{
 if(missing(x))
  stop("x must be a vector")
 if(missing(threshold))
  threshold=0
 if(missing(meanlog))
  meanlog=0
 if(missing(sdlog))
  sdlog=1
 temp=function(x)
 {
	if(x>threshold)
		return((1/(sqrt(2*pi)*sdlog*(x-threshold)))*exp(-(((log((x-threshold))-meanlog)^2)/(2*(sdlog)^2))))
	else
		return(0)
 }
 return(unlist(lapply(x,temp)))
}

plnorm3 <- function(q,meanlog,sdlog,threshold)
{
if(missing(q))
  stop("q must be a vector")
 if(missing(threshold))
  threshold=0
 if(missing(meanlog))
  meanlog=0
 if(missing(sdlog))
  sdlog=1
 temp=function(q)
 {
	if(q>threshold)
		return(pnorm((log((q-threshold))-meanlog)/sdlog))
	else
		return(0)
 }
 return(unlist(lapply(q,temp)))
}

qlnorm3 <- function(p,meanlog,sdlog,threshold,...)
{
 if(missing(p))
  stop("p must be a vector")
 if(missing(meanlog))
  meanlog=0
 if(missing(threshold))
  threshold=0
 if(missing(sdlog))
  sdlog=1
 myfun = function(x, p) plnorm3(q = x,
         meanlog = meanlog, sdlog = sdlog, threshold = threshold) - p
 temp=function(p)
 {
	return(uniroot(f=myfun,lower=threshold,upper=threshold+10000000,p=p,...)$root)       #solve myfun=0
 }
 return(unlist(lapply(p,temp)))
}

dgamma3 <- function(x,shape,scale,threshold)                                    #### GAMMA3
{
 if(missing(x))
  stop("x must be a vector")
 if(missing(threshold))
  threshold=0
 if(missing(shape))
  shape=1
 if(missing(scale))
  scale=1
 temp=function(x)
{
	if(x>=threshold)
		return( dgamma(x-threshold,shape,scale) )
		return(0)
}
 return(unlist(lapply(x,temp)))
}

pgamma3 <- function(q,shape,scale,threshold)
{
if(missing(q))
  stop("q must be a vector")
 if(missing(threshold))
  threshold=-2
 if(missing(shape))
  shape=1
 if(missing(scale))
  scale=1
 temp=function(q)
 {
	if(q>=threshold)
		return(pgamma(q-threshold,shape,scale))
	else
		return(0)
 }
 return(unlist(lapply(q,temp)))
}

qgamma3 <- function(p,shape,scale,threshold,...)
{
 if(missing(p))
  stop("p must be a vector")
 if(missing(threshold))
  threshold=0
 if(missing(shape))
  shape=1
 if(missing(scale))
  scale=1
 myfun = function(x,p) pgamma3(q = x,
         threshold = threshold, scale = scale, shape = shape) - p
 temp=function(p)
 {
  return(uniroot(f=myfun,lower=threshold,upper=threshold+10000000,p=p,...)$root)       #solve myfun=0
 }
 return(unlist(lapply(p,temp)))
}