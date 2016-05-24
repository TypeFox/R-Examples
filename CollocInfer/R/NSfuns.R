make.NS <- function()
{

#  Last modified 25 April 2010

  #  -----------------------------------------------------------------------------

  NSfun = function(times, x, p, more)
  {
    #  Evaluate the right side of the ODE for fitting
    #  the North Shore data

    m2 = 0
    betabasis  = more$betabasis
    betanbasis = betabasis$nbasis
    m1 = m2 + 1
    m2 = m2 + betanbasis
    betacoef = p[m1:m2]
    betamat  = eval.basis(times, betabasis)
    betavec  = betamat %*% betacoef

    alphabasis  = more$alphabasis
    alphanbasis = alphabasis$nbasis
    m1 = m2 + 1
    m2 = m2 + alphanbasis
    alphacoef = p[m1:m2]
    alphamat  = eval.basis(times, alphabasis)
    alphavec  = alphamat %*% alphacoef

    rain = eval.fd(times, more$rainfd)

    f = -betavec*x + alphavec*rain

    return(f)

  }

  #  -----------------------------------------------------------------------------

  NSdfdx = function(times, x, p, more)
  {
    #  Evaluate the right side of the ODE for fitting
    #  the North Shore data

    #  Last modified 24 March 2010

    betabasis  = more$betabasis
    betanbasis = betabasis$nbasis
    betamat    = eval.basis(times, betabasis)
    betacoef   = p[1:betanbasis]
    betavec    = betamat %*% betacoef
    nobs = length(times)
    
    dfdx = array(0,c(nobs,1,1))
    dfdx[,1,1] = -betavec 

    return(dfdx)

  }

  #  -----------------------------------------------------------------------------

  NSdfdp = function(times, x, p, more)
  {
    #  Evaluate the first parameter derivative of the
    #  right side of the ODE for fitting
    #  the North Shore data

    n = length(times)

    betabasis   = more$betabasis
    betanbasis  = betabasis$nbasis
    phimat      = eval.basis(times, betabasis)
    dfdp1       = -phimat * matrix(x,n,betanbasis)

    alphabasis  = more$alphabasis
    alphanbasis = alphabasis$nbasis
    phimat      = eval.basis(times, alphabasis)
    rain        = eval.fd(times, more$rainfd)
    dfdp2       = phimat * matrix(rain,n,alphanbasis)

    nobs        = length(times)
    npar        = length(p)
    dfdp        = array(0,c(nobs,1,npar))
    dfdp[,1,]   = cbind(dfdp1,dfdp2)

    return(dfdp)

  }

  #  -----------------------------------------------------------------------------

  NSd2fdx2 = function(times, x, p, more)
  {
    #  Evaluate the right side of the ODE for fitting
    #  the North Shore data

     nobs   = length(times)
     d2fdx2 = array(0,c(nobs,1,1,1))
     
     return(d2fdx2)
     
  }

  #  -----------------------------------------------------------------------------

  NSd2fdxdp = function(times, x, p, more)
  {
    #  Evaluate the first parameter derivative of the
    #  right side of the ODE for fitting
    #  the North Shore data

    betabasis = more$betabasis
    phimat    = eval.basis(times, betabasis)

    alphanbasis = more$alphabasis$nbasis

    nobs = length(times)
    npar = length(p)
    d2fdxdp = array(0,c(nobs,1,1,npar))
    d2fdxdp[,1,1,] = cbind(-phimat, matrix(0, nobs, alphanbasis))

    return(d2fdxdp)

  }
  
  #  return the named list

  return(
      list(
          fn      = NSfun,
          dfdx    = NSdfdx,
          dfdp    = NSdfdp,
          d2fdx2  = NSd2fdx2,
          d2fdxdp = NSd2fdxdp
      )
    )

}



