NPL.bands <- function(x,conf.level=c(0.95,0.99)){
# Function to compute nonparametric likelihood confidence bands
# as described by Owen, A. (1997) JASA 90:516-521.
# Requires a vector of numeric values 
# and a confidence level (either 0.95 or 0.99).
# The function stores the limits in a list of three variables
# lower and upper along with the distinct values  of x.
  if (!is.numeric(x)) stop("argument must be numeric")
  yn <- table(x) # table of frequencies of distinct value (sorted)
  yi <- as.numeric(names(yn)) # distinct values
  cn <- as.numeric(cumsum(yn)) # cumulated counts
  nn <- rep(sum(yn)+1,length(cn))
  p <- as.numeric(cn/nn)
  if (conf.level==0.95)
    {
      lambda <- ifelse(nn<=100,(3.0123+0.4835*log(nn)-0.00957*(log(nn))^2
                         -0.001488*(log(nn))^3),
                       (3.0806+0.4894*log(nn)-0.02086*(log(nn))^2))}
  else {if (conf.level==0.99)
          {lambda <- ifelse(nn<=100,(-4.626-0.541*log(nn)+0.0242*(log(nn))^2),
                            (-4.71-0.512*log(nn)+0.0219*(log(nn))^2))}
  else stop("Must be either 0.95 or 0.99")}
  lambda <- sqrt(2*lambda)
  phi <- pbeta(p,1/3,1/3)
  se <- 1/(5.3*sqrt(nn*(p*(1-p))^(1/3)))
  phiu <- phi + lambda*se
  phiu <- ifelse(phiu>1,1,phiu)
  phil <- phi - lambda*se
  phil <- ifelse(phil<0,0,phil)
  pu <- qbeta(phiu,1/3,1/3)
  pl <- qbeta(phil,1/3,1/3)
  list(x=yi,lower=pl,upper=pu)
}

  
