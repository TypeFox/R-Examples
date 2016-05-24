#Support functions
############################################
#Simplification functions
L <- function(t,sigma,k)
{
	y <- sigma^2/(2*k)*(1-exp(-2*k*t))
	return(y)
}
K <- function(t,k)
{
	y <- (1-exp(-k*t))/k
	return(y)
}
Tau <- function(t,vol,volActions,k,rho)
{
	y <- (vol^2/k^2+2*rho*vol*volActions/k+volActions^2)*t
	y <- y+2*(vol^2/k^2-rho*volActions*vol/k)*K(t,k)
	y <- y-L(t,vol,k)/k^2
	return(y)
}
############################################
#LMN model functions
A <- function(alpha,beta,sigma,t)
{
	phi <- sqrt(2*sigma^2+beta^2)
	kappa <- (beta+phi)/(beta-phi)
	y1 <- exp(alpha*(beta+phi)/sigma^2*t)
	y2 <- ((1-kappa)/(1-kappa*exp(phi*t)))^(2*alpha/sigma^2)
	y <- y1*y2
	return(y)
}
B <- function(beta,sigma,t)
{
	phi <- sqrt(2*sigma^2+beta^2)
	kappa <- (beta+phi)/(beta-phi)
	y1 <- (beta-phi)/sigma^2
	y2 <- 2*phi/(sigma^2*(1-kappa*exp(phi*t)))
	y <- y1+y2
	return(y)
}
C <- function(eta,t)
{
	y <- exp(eta^2*t^3/6)
}
G <- function(alpha,beta,sigma,t)
{
	phi <- sqrt(2*sigma^2+beta^2)
	kappa <- (beta+phi)/(beta-phi)
	y1 <- alpha/phi*(exp(phi*t)-1)
	y2 <- exp(alpha*(beta+phi)/sigma^2*t)
	y3 <- ((1-kappa)/(1-kappa*exp(phi*t)))^(2*alpha/sigma^2+1)
	y <- y1*y2*y3
	return(y)
}
H <- function(alpha,beta,sigma,t)
{
	phi <- sqrt(2*sigma^2+beta^2)
	kappa <- (beta+phi)/(beta-phi)
	y1 <- exp((alpha*(beta+phi)+phi*sigma^2)/sigma^2*t)
	y2 <- ((1-kappa)/(1-kappa*exp(phi*t)))^(2*alpha/sigma^2+2)
	y <- y1*y2
	return(y)
}
############################################
#Bond support function
Beta <- function(nc,t,T)
{
	y <- floor(nc*(T-t)+1)
	return(y)
}
Alpha <- function(nc,i,T)
{
	y <- T-i/nc
	return(y)
}
CouponTimeValue <- function(couponRate,nc,t,T)
{
	betat <- Beta(nc,t,T)
	alpha <- Alpha(nc,betat,T)
	return(couponRate*(t-alpha+3/360))
}
#Forward rates extraction
#ZC shall be a vector containing all ZC rates from file
############################################
ForwardExtraction <- function(ZC,horizon)
{
	t <- seq(from = 1/12, to=horizon+13/12, by=1/12)
	n  <- length(t)
	ForwardRates <- (ZC[2:n]*t[2:n] - ZC[1:(n-1)]*t[1:(n-1)])*12
	ForwardRates <- ForwardRates[!is.na(ForwardRates)]
	#reduction of the vector
	forwardFinal <- rep(0,(horizon+1))
	forwardFinal <- ForwardRates[(1:(horizon+1))*12]  
	return(forwardFinal)
}

############################################
ZCExtraction <- function(ZC,horizon)
{
	#ZC(0)=0,ZC(1)=R(0,1/12),...
	ZCRate <- ZC[1:(12*horizon+1)]
	return(ZCRate)
}