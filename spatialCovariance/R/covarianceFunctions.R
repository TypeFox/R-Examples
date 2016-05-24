defineK <- function(class,K)
  {
    ## current set of covariance functions
    ## analytic results exist for the power and ldt models
    ## the powerNI and ldtNI are for numerical integration

    ## K gets passed in only when class=misc
    
    if(class=="matern") K <- function(d,params) (params[1]*d)^params[2]*besselK(params[1]*d,params[2])
    if(class=="bess1") K <- function(d,params) (params[1]*d)*besselK(params[1]*d,1)
    if(class=="exp") K <- function(d,params) exp(-params[1]*d)
    ##if(class=="exp") K <- function(d,params) sqrt(params[1]*d)*besselK(params[1]*d,0.5)
    if(class=="bess0") K <- function(d,params) besselK(params[1]*d,0)
    
    if(class=="powerNI") K <- function(d,params) d^(-params[1])/params[1]
    ## params < 2 above
    if(class=="powerNI") K <- function(d,params) -d^(2*params[1])/params[1]
    ## range for params is params > -1 in power model

    if(class=="spline") K <- function(d,params) d^2*log(d)
    if(class=="ldtNI") K <- function(d,params) -log(d)
    if(class=="cauchy") K <- function(d,params) 1/(params[2])*(1+abs(params[1]*d)^params[3])^(-params[2]/params[3])
    ## params are (lambda, beta, alpha)

    ## derivatives of covariance functions
    if(class=="dmatern1") K <- function(d,params) d * (d*params[1])^(-1+params[2])*params[2]*besselK(params[1]*d,params[2])+d/2*(d*params[1])^(params[2])*(-besselK(params[1]*d,-1+params[2])-besselK(params[1]*d,1+params[2]))

    if(class=="dbess1") K <- function(d,params) d*besselK(params[1]*d,1)+d/2*(d*params[1])*(-besselK(params[1]*d,0)-besselK(params[1]*d,2))
    if(class=="dexp") K <- function(d,params) -d*exp(-params[1]*d)
    if(class=="dexp2") K <- function(d,params) d*d*exp(-params[1]*d)

    if(class=="dbess0") K <- function(d,params) -d*besselK(params[1]*d,1)
    if(class=="dpowerNI") K <- function(d,params) d^(2*params[1])/params[1]^2-2*d^(2*params[1])*log(d)/params[1]

    cov.f <- function(d,rw,cw,ax,bx,i,j,params,K)
      {
        K(d,params)*f(d,rw,cw,ax,bx,i,j)
      }
    
    list(K=K,cov.f=cov.f)
  }


