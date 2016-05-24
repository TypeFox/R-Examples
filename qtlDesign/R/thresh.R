##################################################################
# function for calculating tail probabilities of genome scans
##################################################################
# t = threshold on LOD (base 10) scale
# G = genome length in centiMorgans
# cross = cross type, "bc" or "f2"
# type = type of LOD score
#      "1": one-dim scan
#      "2o": two-dim scan, full model
#      "2i": two-dim scan, interaction only
#      "2a": two-dim scan, additive only
#      "3a": three-dim scan, additive only    
# d = marker spacing in centiMorgans
# cov.dim = dimension of interacting covariate

tailprob <- function(t,G,cross,type="1",d=0.01,cov.dim=0)
  {

    # convert to Morgans
    G <- G/100
    d <- d/100
    
    # convert LOD to chi-squared scale
    t <- t*log(10)*2

    # local drift rate
    # mu <- switch(cross,bc=sqrt(2*t*d),f2=sqrt(3*t*d))
    mu <- switch(cross,bc=sqrt(4*t*d),f2=sqrt(6*t*d))    

    # correction for terms not contributing to drift

    # backcross
    if((cross=="bc")&&(type=="2o"))
      {
        mu <- mu*sqrt(2/3)
      }

    if((cross=="bc")&&(type=="2x"))
      {
        mu <- mu
      }

    # intercross
    if((cross=="f2")&&(type=="2o"))
      {
        mu <- mu*sqrt(6/8)
      }

    if((cross=="f2")&&(type=="2x"))
      {
        mu <- mu
      }

    # c <- 2*mu^2*nu(2*mu)
    # c <- mu^2*nu(mu)/2    
    # c <- tau(mu)
    
    
    # calculate space volume and clump volume from search dimension

    # one-dimensional search
    if(regexpr("1",type)==1)
      {
        c <- mu^2*nu(mu)/2 
        S <- (G/d)
      }

    # two-dimensional search
    if(regexpr("2",type)==1)
      {
	if(type=="2x")
	{
	if(cross=="bc")
  	  c <- (mu^2*nu(mu)/2)*((mu/2)^2*nu(mu/2)/2)
	if(cross=="f2")
  	  c <- (mu^2*nu(mu)/2)*((mu*2/3)^2*nu(mu*2/3)/2)
	}
	else
	{
	c <- mu^2*nu(mu)/2 	
        c <- c^2
	}
        S <- (G/d)^2/2
      }

    if(regexpr("3",type)==1)
      {
        c <- c^3
        S <- (G/d)^3/6
      }

    


    # degrees of freedom of marginal test statistic
    df <- (1+cov.dim) * switch(paste(cross,type,sep="."),
                               bc.1=1, bc.2i=1, bc.2a=2, bc.2o=3, bc.2x=2,
                               f2.1=2, f2.2i=4, f2.2a=4, f2.2o=8, f2.2x=6,
                               bc.3a=3,f2.3a=6)

    # put everything together: search space volume, (inverse of)
    # expected clump size and marginal tail probability

    lam <- S*c*pchisq(t,df,low=FALSE)
    # lam <- 2*S*c*dchisq(t,df)
    
    # exponentiate for tail probability
    1-exp(-lam)
  }


##################################################################
# function for calculating thresholds
##################################################################
# G = genome length in centiMorgans
# cross = cross type, "bc" or "f2"
# type = type of LOD score
#      "1": one-dim scan
#      "2o": two-dim scan, full model
#      "2i": two-dim scan, interaction only
#      "2a": two-dim scan, additive only
#      "3a": three-dim scan, additive only
# p = vector of genome-wide significance levels
# d = marker spacing in centiMorgans
# cov.dim = dimension of interacting covariate
# interval = LOD interval to be searched for threshold
    

# calculate thresholds from tail probabilities
thresh <- function(G,cross,type="1",p=c(0.10,0.05,0.01),
                        d=0.01,cov.dim=0,interval=c(1,40))
  {

    thresh <- vector(mode="numeric",length=length(p))
    for( i in 1:length(p) )
      {
        thresh[i] <- uniroot(function(x){
          tailprob(x,G,cross,type,d,cov.dim) - p[i]},interval=interval)$root
      }
    thresh 

  }

# the nu function of Siegmund

nu <- function(mu,k.lim=10000)
  {
    if(abs(mu)<0.1)
      {
        res <- -0.583*mu
      }
    else
      {
        k <- 1:k.lim
        zzz <- pnorm( -0.5 * mu * sqrt( k ) ) / k
        res <- -2*sum(zzz) + log(2) - 2*log(mu)
      }
    exp(res)
  }

tau <- function(mu)
  {
    2*mu^2*nu(2*mu)
  }
