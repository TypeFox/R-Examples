########## Tools
# simply coumpounded euribor rate
euriborfromprice <- function(t, T, ZC)
{
  # Brigo P. 7
  tau <- T-t
  return(as.vector((1/ZC - 1)/tau))
}
euriborfromprice <- cmpfun(euriborfromprice)


# simply coumpounded zero-coupon price
pricefromeuribor <- function(t, T, L)
{
  # Brigo P. 7
  tau <- T-t
  return(as.vector(1/(1 + L*tau)))
}
pricefromeuribor <- cmpfun(pricefromeuribor)

# simply coumpounded forward rate
fwdrate <- function(T_, S, ZC_T, ZC_S)
{
  # Brigo P. 12
  tau <- S-T_
  return((ZC_T/ZC_S - 1)/tau)
}
fwdrate <- cmpfun(fwdrate)

# Swap curve bootstrap
bootstrapswapcurve <- function(swaprate, maturity, typeres=c("rates", "prices"))
{  
  swaprate <- swaprate/100
  accrual <- diff(maturity)
  n <- length(swaprate)  
  ZC <- numeric(n)  
  ZC[1] <- 1/(1+maturity[1]*swaprate[1])
  
  for (i in seq_len(n)[-1])
  {
    ind <- seq_len(i-1)
    ZC[i] <- (1 - swaprate[i]*sum(accrual[ind]*ZC[ind]))/(1+swaprate[i]*accrual[i-1])
  }
  
  typeres <- match.arg(typeres)
  
  if (missing(typeres) || typeres == "prices")
  {return(ts(ZC))}
  
  if (typeres == "rates")
  {return((1-ZC)/(maturity*ZC))}  
}
bootstrapswapcurve <- cmpfun(bootstrapswapcurve)

########## Interpolation
# Smith-Wilson
tZC_SW <- function(yM = NULL, p = NULL, u, t, UFR, typeres=c("rates", "prices"), T_UFR = NULL)
{
  N <- length(u)  
  J <- length(t)
  
  fonctionsWilson <- function(t,u,alpha,UFR)
  {    
    N <- length(u)  
    J <- length(t)    
    u_mat <- matrix(rep.int(u,J), nrow=J, byrow=T)
    t_mat <- t(matrix(rep.int(t,N), nrow=N, byrow=T))
    min_u <- u_mat*(u_mat<=t_mat) + t_mat*(u_mat>t_mat)
    max_u <- u_mat + t_mat - min_u    
    return(exp(-UFR*(u_mat+t_mat))*(alpha*min_u - 0.5*exp(-alpha*max_u)*(exp(alpha*min_u)-exp(-alpha*min_u))))
  }
 
  if ((is.null(p) || missing(p)) && (!is.null(yM) || !missing(yM)))
  {
    p <- pricefromeuribor(t=0, T=u, L=yM)
  }
  
  alpha <- 0.1
  mu <- exp(-UFR*u)
  W <- fonctionsWilson(u,u,alpha,UFR)    
  Xhi <- solve(W)%*%(p-mu)
  W_interp <- fonctionsWilson(t,u,alpha,UFR)  
  P <- exp(-UFR*t) + W_interp %*%Xhi
  Fwd <- fwdrate(t[-J], t[-1], P[-J], P[-1])
  
  typeres <- match.arg(typeres)  
  if(max(t) > max(u))
  {
    if (missing(T_UFR)) stop("For the extrapolation T_UFR must be provided, check the output maturities")        
    alpha <- 0.1
    
    if (max(u)+T_UFR > max(t)) stop("Not enough extrapolation dates. Try a lower T_UFR")
    
    while(abs(Fwd[pmatch(max(u)+T_UFR, t)]- UFR) >= 0.0003)
    {
          alpha <- alpha+0.001
          mu <- exp(-UFR*u)
          W <- fonctionsWilson(u,u,alpha,UFR)    
          Xhi <- solve(W)%*%(p-mu)
          W_interp <- fonctionsWilson(t,u,alpha,UFR)  
          P <- exp(-UFR*t) + W_interp %*%Xhi      
          Fwd <- fwdrate(t[-J], t[-1], P[-J], P[-1])
    }    
  }
    
  if (typeres == "prices")
  {return (list(coefficients = list(alpha=alpha, Xhi=as.vector(Xhi)), 
                values = as.vector(P),
                fwd = Fwd))}
  
  if (typeres == "rates")
  {
    return (list(coefficients = list(alpha=alpha, Xhi=as.vector(Xhi)), values = euriborfromprice(t=0, T=t, ZC=P), fwd = Fwd))
  }
}
tZC_SW <- cmpfun(tZC_SW)

# Nelson-Siegel interpolation
tZC_NS <- function(yM = NULL, p = NULL, matsin, matsout, typeres=c("rates", "prices"))
{
#   require(randtoolbox)
#   require(mcGlobaloptim)
  
  N <- length(matsin)  
  J <- length(matsout)
  
  if ((!is.null(p) || !missing(p)) && (is.null(yM) || missing(yM)))
  {
    yM <- euriborfromprice(t=rep.int(0, N), T=matsin, ZC=p)
  }
  
  # generic objective function
  OF <- function(betaV)
  {     
    gam1 <- matsin/betaV[4]
    aux1 <- 1-exp(-gam1)     
    y <- betaV[1]+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)
    return(crossprod(y - yM))
  }
  OF <- cmpfun(OF)
  set.seed(123)
  betaV <- multiStartoptim(objectivefn = OF,
                           lower = c(0, -15, -30, 0),
                           upper = c (15, 30, 30, 3),
                           method = "nlminb",
                           nbtrials = 300, 
                           typerunif = "niederreiterlodisp")$par
  gam1 <- matsout/betaV[4]
  aux1 <- 1-exp(-gam1)
  L <- betaV[1]+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)
  aux2 <- 1 - aux1
  Fwd <- betaV[1] + betaV[2]*aux2 + betaV[3]*gam1*aux2
  
  if (typeres == "rates")
  {return (list(coefficients = betaV, 
                values = L, fwd = Fwd))}
  
  if (typeres == "prices")
  {
    return (list(coefficients = betaV, 
                 values = pricefromeuribor(t=0, T=matsout, L=L), fwd = Fwd))}
}
tZC_NS <- cmpfun(tZC_NS)

# Svensson interpolation
tZC_SV <- function(yM = NULL, p = NULL, matsin, matsout, typeres=c("rates", "prices"))
{
#   require(randtoolbox)
#   require(mcGlobaloptim)
  
  N <- length(matsin)  
  J <- length(matsout)
  
  if ((!is.null(p) || !missing(p)) && (is.null(yM) || missing(yM)))
  {
    yM <- euriborfromprice(t=0, T=matsin, ZC=p)
  }
  
    # generic objective function
    OF <- function(betaV)
    {   	
      gam1 <- matsin/betaV[5]
      gam2 <- matsin/betaV[6]
      aux1 <- 1-exp(-gam1)
      aux2 <- 1-exp(-gam2)      
      y <- betaV[1]+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)+
        betaV[4]*(aux2/gam2+aux2-1)
       
      return(crossprod(y - yM))
    }
    OF <- cmpfun(OF)
  set.seed(123)    
  betaV <- multiStartoptim(objectivefn = OF,
                         lower = c(0, -15, -30, -30 ,0 ,3),
                         upper = c (15, 30, 30, 30 ,3 ,6),
                         method = "nlminb",
                         nbtrials = 300, 
                         typerunif = "niederreiterlodisp")$par
  gam1 <- matsout/betaV[5]
  gam2 <- matsout/betaV[6]
  aux1 <- 1-exp(-gam1)
  aux2 <- 1-exp(-gam2)  
  L <- betaV[1]+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)+
      betaV[4]*(aux2/gam2+aux2-1)
  aux3 <- 1 - aux1
  aux4 <- 1 - aux2
  Fwd <- betaV[1] + betaV[2]*aux3 + betaV[3]*gam1*aux3 + betaV[4]*gam2*aux4
  
    if (typeres == "rates")
    {return (list(coefficients = betaV, 
                  values = L, fwd = Fwd))}
    
    if (typeres == "prices")
    {
      return (list(coefficients = betaV, 
                   values = pricefromeuribor(t=rep.int(0, J), T=matsout, L=L), fwd = Fwd))
    }
}
tZC_SV <- cmpfun(tZC_SV)

# Polynomial Hermite cubic spline interpolation
hermitecubicspline <- function(yM = NULL, p = NULL, matsin, matsout, typeres=c("rates", "prices"))
{
  u <- matsin
  t <- matsout
  
  if (max(t) > max(u)) stop("Only interpolation can be performed with cubic splines, check the output maturities")
    
  if (!is.null(yM) && !is.null(p)) stop("either yields OR prices must be provided")
  
  if (is.null(yM) && is.null(p)) stop("either yields or prices must be provided")
  
  if (!is.null(yM))
  {Sw <- yM}
  else{Sw <- p}
  
  n <- length(u)
  i <- seq_len(n)
  iup <- i[-1]
  idown <- i[-n]  
  Delta <- diff(u)
  delta <- diff(Sw)

  A <- delta/Delta    
  A_up <- A[iup]
  A_down <- A[idown] 
  
  Delta_up <- Delta[iup]
  Delta_down <- Delta[idown] 
  
  P <- (A_up*Delta_down + A_down*Delta_up)/(Delta_down + Delta_up)
  P <- c(2*A[1] - P[2], P[!is.na(P)], 0)
    
  m <- length(t)
  z <- rep.int(0, m)
  SWout <- rep.int(0, m)
  
  for (i in  seq_len(m))
  {          
    u_down <-  pmatch(floor(t[i]),u)
    u_up <-  pmatch(ceiling(t[i]),u)
    Sw_down <- Sw[u_down]
    Sw_up <- Sw[u_up]
    P_down <- P[u_down]
    P_up <- P[u_up]
    Delta_down <- Delta[u_down]
    Delta_up <- Delta[u_up]
            
    if (Sw_down != Sw_up) 
    {
      z <- (t[i] - u_down)/(u_up - u_down)                    
      SWout[i] <- Sw_down*((1-z)^3) + (3*Sw_down + Delta_down*P_down)*
        ((1-z)^2)*z + (3*Sw_up - Delta_down*P_up)*(1-z)*(z^2) + 
        Sw_up*(z^3)      
    }
    else 
    {
      SWout[i] <- Sw_down
    }    
  }
  
  if (typeres=="rates")
  {P <- pricefromeuribor(0, t, SWout)}
  else {P <- SWout}
  
  Fwd <- fwdrate(t[-m], t[-1], P[-m], P[-1])
  
  if ((!is.null(yM) && typeres == "rates") || (!is.null(p) && typeres == "prices"))
  {return (list(coefficients = NA, values = SWout, fwd = Fwd))}
  
  if (!is.null(yM) && typeres == "prices")
  {
    return (list(coefficients = NA, values = pricefromeuribor(0, t, SWout), fwd = Fwd))
  }

  if (!is.null(p) && typeres == "rates")
  {
    return (list(coefficients = NA, values = euriborfromprice(0, t, SWout), fwd = Fwd))
  }
}
hermitecubicspline <- cmpfun(hermitecubicspline)

# Yield curve interpolation
ycInterpolation <- function(yM = NULL, p = NULL, matsin, matsout, 
                            method=c("NS", "SV", "SW", "HCSPL"), typeres=c("rates", "prices"))
{
   if (missing(yM) && missing(p)) stop("market yield to maturities or prices must be provided")
  
  if (missing(matsin) || missing(matsout)) stop("Input and output maturities must be provided")
  
  if (max(matsout) > max(matsin)) warning("Try extrapolation function instead. Check the output maturities.")
  
  method <- match.arg(method)
  
  typeres <- match.arg(typeres) 
  
  if (is.null(yM) || missing(yM))
  {
    input_curve <- p
    # tZC_SW(UFR = UFR, T_UFR = T_UFR, u = matsin, p = input_curve, t = matsout, typeres = typeres)
    switch(method,
           NS = tZC_NS(p = input_curve, matsin = matsin, matsout = matsout, typeres = typeres),
           SV = tZC_SV(p = input_curve, matsin = matsin, matsout = matsout, typeres = typeres),
           SW = tZC_SW(UFR = 0.042,u = matsin, p = input_curve, t = matsout, typeres = typeres),
           HCSPL = hermitecubicspline(p=input_curve, matsin = matsin, matsout = matsout, typeres = typeres))}
  else {
        input_curve <- yM
        switch(method,
               NS = tZC_NS(yM = input_curve, matsin = matsin, matsout = matsout, typeres = typeres),
               SV = tZC_SV(yM = input_curve, matsin = matsin, matsout = matsout, typeres = typeres),
               SW = tZC_SW(UFR=0.042,u = matsin, yM = input_curve, t = matsout, typeres = typeres),
               HCSPL = hermitecubicspline(yM=input_curve, matsin = matsin, matsout = matsout, typeres = typeres))
  }
}
ycInterpolation <- cmpfun(ycInterpolation)

########## Extrapolation

# Nelson-Siegel
tZC_NS_extra <- function(yM = NULL, p = NULL, matsin, matsout, typeres=c("rates", "prices"), UFR)
{
#   require(randtoolbox)
#   require(mcGlobaloptim)
  
  N <- length(matsin)  
  J <- length(matsout)
  
  if ((!is.null(p) || !missing(p)) && (is.null(yM) || missing(yM)))
  {
    yM <- euriborfromprice(t=rep.int(0, N), T=matsin, ZC=p)
  }
  
  # generic objective function
  OF <- function(betaV)
  {     
    gam1 <- matsin/betaV[4]
    aux1 <- 1-exp(-gam1)     
    y <- UFR+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)
    return(crossprod(y - yM))
  }
  OF <- cmpfun(OF)
  set.seed(123)
  betaV <- multiStartoptim(objectivefn = OF,
                           lower = c(UFR, -15, -30, 0),
                           upper = c (UFR, 30, 30, 3),
                           method = "nlminb",
                           nbtrials = 300, 
                           typerunif = "niederreiterlodisp")$par
  gam1 <- matsout/betaV[4]
  aux1 <- 1-exp(-gam1)
  L <- UFR+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)
  aux2 <- 1 - aux1
  Fwd <- UFR + betaV[2]*aux2 + betaV[3]*gam1*aux2
  
  if (typeres == "rates")
  {return (list(coefficients = betaV, 
                values = L, fwd = Fwd))}
  
  if (typeres == "prices")
  {
    return (list(coefficients = betaV, 
                 values = pricefromeuribor(t=0, T=matsout, L=L), fwd = Fwd))
  }
}
tZC_NS_extra <- cmpfun(tZC_NS_extra)

# Svensson
tZC_SV_extra <- function(yM = NULL, p = NULL, matsin, matsout, typeres=c("rates", "prices"), UFR)
{  
#   require(randtoolbox)
#   require(mcGlobaloptim)
  
  N <- length(matsin)  
  J <- length(matsout)
  
  if ((!is.null(p) || !missing(p)) && (is.null(yM) || missing(yM)))
  {
    yM <- euriborfromprice(t=0, T=matsin, ZC=p)
  }
  
  # generic objective function
  OF <- function(betaV)
  {     
    gam1 <- matsin/betaV[5]
    gam2 <- matsin/betaV[6]
    aux1 <- 1-exp(-gam1)
    aux2 <- 1-exp(-gam2)      
    y <- UFR+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)+
      betaV[4]*(aux2/gam2+aux2-1)
    return(crossprod(y - yM))
  }
  OF <- cmpfun(OF)
  set.seed(123)
  betaV <- multiStartoptim(objectivefn = OF,
                           lower = c(UFR, -15, -30, -30 ,0 ,3),
                           upper = c (UFR, 30, 30, 30 ,3 ,6),
                           method = "nlminb",
                           nbtrials = 300, 
                           typerunif = "niederreiterlodisp")$par
  gam1 <- matsout/betaV[5]
  gam2 <- matsout/betaV[6]
  aux1 <- 1-exp(-gam1)
  aux2 <- 1-exp(-gam2)  
  L <- UFR+betaV[2]*(aux1/gam1)+betaV[3]*(aux1/gam1+aux1-1)+
    betaV[4]*(aux2/gam2+aux2-1)
  aux3 <- 1 - aux1
  aux4 <- 1 - aux2
  Fwd <- UFR + betaV[2]*aux3 + betaV[3]*gam1*aux3 + betaV[4]*gam2*aux4
  
  if (typeres == "rates")
  {return (list(coefficients = betaV, 
                values = L, fwd = Fwd))}
  
  if (typeres == "prices")
  {
    return (list(coefficients = betaV, 
                 values = pricefromeuribor(t=rep.int(0, J), T=matsout, L=L), fwd = Fwd))
  }
}
tZC_SV_extra <- cmpfun(tZC_SV_extra)

# Yield curve extrapolation wrapper
ycExtrapolation <- function(yM = NULL, p = NULL, matsin, matsout, 
                            method=c("NS", "SV", "SW"), typeres=c("rates", "prices"), UFR, T_UFR = NULL)
{  
  if (missing(yM) && missing(p)) stop("market yield to maturities or prices must be provided")  
  if (missing(matsin) || missing(matsout)) stop("Input and output maturities must be provided")  
  if (max(matsout) <= max(matsin)) warning("max(matsout) <= max(matsin) : try the interpolation function instead")  
  method <- match.arg(method)  
  typeres <- match.arg(typeres) 
        
  if (is.null(yM) || missing(yM))
  {
    input_curve <- p    
    switch(method,
           NS = tZC_NS_extra(p = input_curve, matsin = matsin, matsout = matsout, typeres = typeres, UFR = UFR),
           SV = tZC_SV_extra(p = input_curve, matsin = matsin, matsout = matsout, typeres = typeres, UFR = UFR),
           SW = tZC_SW(UFR = UFR, T_UFR = T_UFR, u = matsin, p = input_curve, t = matsout, typeres = typeres)) 
  }
  else {
        input_curve <- yM
        switch(method,
               NS = tZC_NS_extra(yM = input_curve, matsin = matsin, matsout = matsout, typeres = typeres, UFR = UFR),
               SV = tZC_SV_extra(yM = input_curve, matsin = matsin, matsout = matsout, typeres = typeres, UFR = UFR),
               SW = tZC_SW(UFR = UFR, T_UFR = T_UFR, u = matsin, yM = input_curve, t = matsout, typeres = typeres)) 
        }
  
}
ycExtrapolation <- cmpfun(ycExtrapolation)