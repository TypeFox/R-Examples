################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.3                                        Peter Reichert 05.10.2014 #
#                                                                              #
################################################################################


# ==============================================================================
# utility aggregation functions
# ==============================================================================


utility.aggregate.add <- function(u,par)  # par[i]: weight of u[i]
{
  # check input:
  
  if ( length(u) != length(par) )
  {
    warning("Length of utilities/values and weights not equal: ",
            length(u)," ",length(par))
    return(NA)
  }
  ind <- which(!is.na(u))
  if ( length(ind) == 0 ) return(NA)
  if ( sum( par < 0 ) > 0 )
  {
    warning("Parameter of additive aggregation smaller than zero")
    return(NA)
  }
  
  # calculate aggregated value
  
  s <- sum(par[ind])
  if ( s <= 0 ) return(NA)
  u.agg <- sum(par[ind]*u[ind])/s
  
  return(as.numeric(u.agg))
}


utility.aggregate.min <- function(u,par=NA)
{
  # check input:
  
  ind <- which(!is.na(u))
  if ( length(ind) == 0 ) return(NA)
  
  # calculate aggregated value
  
  u.agg <- min(u[ind])
  
  return(as.numeric(u.agg))
}


utility.aggregate.max <- function(u,par=NA)
{
  # check input:
  
  ind <- which(!is.na(u))
  if ( length(ind) == 0 ) return(NA)
  
  # calculate aggregated value
  
  u.agg <- max(u[ind])
  
  return(as.numeric(u.agg))
}


utility.aggregate.mult <- function(u,par) 
{
  # check input:
  
  if ( length(u) != length(par) )
  {
    warning("Length of utilities/values and weights not equal: ",
            length(u)," ",length(par))
    return(NA)
  }
  ind <- which(!is.na(u))
  if ( length(ind) == 0 ) return(NA)
  if ( length(ind) == 1 )
  {
    return(as.numeric(u[ind]))
  }
  if ( sum( par < 0 | par > 1 ) > 0 )
  {
    warning("Parameter of multiplicative aggregation",
            "smaller than zero or larger than unity")
    return(NA)
  }
  
  # function used in uniroot to determine the scaling constant k:
  
  utility.aggregate.mult.root <- function(k,ki)
  {
    res <- 1
    for ( i in 1:length(ki) )
    {
      res <- res * ( 1 + k * ki[i] )
    }
    res <- 1 + k - res
    return(res)
  }
  
  # define numerical parameter:
  
  eps <- 1e-3   # maximum deviation of sum(par) from unity to use additive fcn 
  
  # rescale weights:
  
  s <- sum(par)   
  fact <- s/sum(par[ind])
  ki <- fact*par[ind]
  
  # calculate additive utility function if sum close to unity:
  
  if ( s > 1-eps & s < 1+eps )
  {
    return(utility.aggregate.add(u,par))
  }
  
  # calculate multiplicative utility function if sum not close to unity:
  
  # calculate k: 
  # (Keeney and Raiffa, Decisions with multiple objectives, 1976,
  # pp. 307, 347-348)
  
  if ( s < 1 )
  {
    lower <- 1
    i <- 0
    while ( utility.aggregate.mult.root(lower,ki) < 0 )
    {
      lower <- 0.1*lower
      i <- i+1
      if ( i > 20 )
      {
        warning("Problem solving equation for scaling constant")
        return(NA)
      }
    }
    upper <- 1
    i <- 0
    while ( utility.aggregate.mult.root(upper,ki) > 0 )
    {
      upper <- 10*upper
      i <- i+1
      if ( i > 20 )
      {
        warning("Problem solving equation for scaling constant")
        return(NA)
      }
    }
    k <- uniroot(utility.aggregate.mult.root,ki=ki,
                 lower=lower,upper=upper)$root
  }
  else  # s > 1
  {
    upper <- -0.1
    i <- 0
    while ( utility.aggregate.mult.root(upper,ki) < 0 )
    {
      upper <- 0.1*upper
      i <- i+1
      if ( i > 20 )
      {
        warning("Problem solving equation for scaling constant")
        return(NA)
      }
    }
    k <- uniroot(utility.aggregate.mult.root,ki=ki,
                 lower=-1,upper=upper)$root 
  }
  
  # evaluate multiplicative utility function:
  
  u.agg <- 1  
  for ( i in 1:length(ki) )
  {
    if ( !is.na(u[ind][i]) ) u.agg <- u.agg * (k*ki[i]*u[ind][i]+1) 
  }
  u.agg <- (u.agg - 1)/k
  
  # eliminate values out of range due to numerical inaccuracies:
  
  u.agg <- ifelse(u.agg < 0, 0, u.agg)
  u.agg <- ifelse(u.agg > 1, 1, u.agg)
  
  return(as.numeric(u.agg))
}


utility.aggregate.geo <- function(u,par) 
{
  # check input:
  
  if ( length(u) != length(par) )
  {
    warning("Length of utilities/values and weights not equal: ",
            length(u)," ",length(par))
    return(NA)
  }
  ind <- which(!is.na(u))
  if ( length(ind) == 0 ) return(NA)
  if ( sum( par < 0 ) > 0 )
  {
    warning("Parameter of geometric aggregation smaller than zero")
    return(NA)
  }
  
  # calculate aggregated value
  
  s <- sum(par[ind])
  if ( s <= 0 ) return(NA)
  u.agg <- 1
  for ( i in 1:length(ind) )
  {
    if ( par[ind][i]>0 ) u.agg <- u.agg*u[ind][i]^(par[ind][i]/s)
  }
  
  return(as.numeric(u.agg))
}


utility.aggregate.revgeo <- function(u,par) 
{
  return(1-utility.aggregate.geo(1-u,par))
}


utility.aggregate.geooff <- function(u,par)
{
  n <- length(u)
  
  # check input:
  
  if ( length(par) != n + 1)
  {
    warning("Length of parameter vector should be length of utilities/values (for weights) plus one (for offset): ",
            length(par)," ",n)
    return(NA)
  }
  u <- utility.aggregate.geo(u+par[n+1],par[1:n])-par[n+1]
  # correct for numerical errors due to differences of "large" numbers
  u <- ifelse(u>0,u,0)
  u <- ifelse(u<1,u,1)
  return(u) 
}


utility.aggregate.revgeooff <- function(u,par) 
{
  return(1-utility.aggregate.geooff(1-u,par))
}


utility.aggregate.cobbdouglas <- function(u,par) 
{
  return(utility.aggregate.geo(u,par))
}


utility.aggregate.harmo <- function(u,par) 
{
  # check input:
  
  if ( length(u) != length(par) )
  {
    warning("Length of utilities/values and weights not equal: ",
            length(u)," ",length(par))
    return(NA)
  }
  ind <- which(!is.na(u))
  if ( length(ind) == 0 ) return(NA)
  if ( sum( par < 0 ) > 0 )
  {
    warning("Parameter of harmonic aggregation smaller than zero")
    return(NA)
  }
  
  # calculate aggregated value
  
  s <- sum(par[ind])
  if ( s <= 0 ) return(NA)
  if ( sum(u==0) > 0 ) return(0)
  
  u.agg <- s / sum(par[ind]/u[ind])
  
  return(as.numeric(u.agg))
}


utility.aggregate.revharmo <- function(u,par) 
{
  return(1-utility.aggregate.harmo(1-u,par))
}


utility.aggregate.harmooff <- function(u,par)
{
  n <- length(u)
  
  # check input:
  
  if ( length(par) != n + 1)
  {
    warning("Length of parameter vector should be length of utilities/values (for weights) plus one (for offset): ",
            length(par)," ",n)
    return(NA)
  }
  return(utility.aggregate.harmo(u+par[n+1],par[1:n])-par[n+1])
}


utility.aggregate.revharmooff <- function(u,par) 
{
  return(1-utility.aggregate.harmooff(1-u,par))
}


utility.aggregate.mix <- function(u,par)  # par[i]: weight of u[i]
{                                         # par[n+j]: weight of technique j
  # check input:                         # (j = add, min, geo)
  
  n <- length(u)
  if ( n+3 != length(par) )
  {
    warning("Length of parameter vector must be equal to",
            "length of utilities/values plus three:",
            length(par),length(u))
    return(NA)
  }
  s <- sum(par[n+(1:3)])
  if ( s <= 0 | sum(par[n+(1:3)]<0) > 0 )
  {
    warning("Weights of aggregation techniques to average",
            "cannot be negative or not all of them equal to zero")
    return(NA)
  }
  
  u.add <- 0; if ( par[n+1] != 0 ) u.add <- utility.aggregate.add(u,par[1:n])
  u.min <- 0; if ( par[n+2] != 0 ) u.min <- utility.aggregate.min(u)
  u.geo <- 0; if ( par[n+3] != 0 ) u.geo <- utility.aggregate.geo(u,par[1:n])
  
  if ( is.na(u.add) | is.na(u.min) | is.na(u.geo) ) return(NA)
  u.agg <- (par[n+1]*u.add + par[n+2]*u.min + par[n+3]*u.geo)/s
  
  return(u.agg)
}


utility.aggregate.addmin <- function(u,par)
{
  n <- length(u)
  if ( length(par) != n+1 )
  {
    warning("Length of parameter vector should be length of utilities/values ",
            "(for weights) plus one (for weight between methods): ", 
            length(par), " ", n)
    return(NA)
  }
  return(   par[n+1]  * utility.aggregate.add(u,par[1:n]) +
              (1-par[n+1]) * utility.aggregate.min(u,NA))
}

