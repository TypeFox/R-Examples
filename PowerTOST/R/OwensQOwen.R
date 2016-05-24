#-----------------------------------------------------------------------------
# Calculation of Owens Q-function by the algo given by Owen himself:
# repeated integration by parts
# 
# Author: dlabes Mar 2012
#-----------------------------------------------------------------------------

# Owen's T-function: 
# Calculates integral 0 to a of exp(-0.5*h^2*(1+x^2))/(1+x^2)/(2*pi)
# most simple implementation via integrate()
# TODO: consider eventually a numerical algo with high precision
# according to M. PATEFIELD, D. TANDY

OT_integrand <- function(x, h) exp(-0.5*h^2*(1+x^2))/(1+x^2)/(2*pi)

OwensT <- function(h, a)
{ 
  int <- integrate(OT_integrand, lower=0, upper=abs(a), h=h)$value
  # in case of a=Inf or -Inf the condition T(h,-a)=-T(h,a) is not maintained!
  int <- ifelse(a<0, -int, int)
  return(int) 
}

# Owen's Q-function
# Calculates Owen's Q-function via repeated integration by parts
# formulas as given in 
# Owen, D. B. (1965)
# "A Special Case of a Bivariate Non-central t-Distribution"
# Biometrika Vol. 52, p437-446.

OwensQOwen <- function(nu, t, delta, a=0, b)
{
  if (nu<1) stop("nu must be >=1!")
  if (a != 0) stop("Only a=0 implemented!")
  
  A   <- t/sqrt(nu)
  B   <- nu/(nu + t*t)
  upr <- nu-2           # upper index of integration by parts
  # the coefficients a(k)
  av  <- vector(mode="numeric", length=nu)
  for (k in seq_along(av)){
    if (k==1 | k==2) av[k] <- 1 else av[k] <- 1/((k-2)*av[k-1])
  }
  ll <- ifelse((upr-1)>0, upr-1, 0)
  L  <- vector(mode="numeric", length=ll)
  # k-1 of the formulas transformes to k here
  if(is.finite(b)){
    # all L[k] are NaN if b==Inf, since dnorm(Inf)==0 and 0*Inf=NaN
    # we use L[k]=0 instead
    for (k in seq_along(L)){
                       # gives NaN 
      if (k==1) L[1] <- 0.5*A*B*b*dnorm(b)*dnorm(A*b-delta)
        else L[k] <- av[k+3]*b*L[k-1]
    }
  } 
  # the series 0 to k is stored as 1 to k+1
  ll <- ifelse((upr+1)>0, upr+1, 0)
  H  <- vector(mode="numeric", length=ll)
  # k+1 of the formulas transformes to k here
  if(is.finite(b)){
    # since H[1]==0 we also get NaN here if b==Inf  
    for (k in seq_along(H)){
      if (k==1) H[1] <- - dnorm(b)*pnorm(A*b-delta)
        else H[k] <- av[k+1]*b*H[k-1]
    }
  }  
  M    <- vector(mode="numeric", length=ll)
  sB   <- sqrt(B)
  # k+1 in the formulas transformes to k here
  for (k in seq_along(M)){
    if (k==1) M[1] <- A*sB*dnorm(delta*sB)*( pnorm(delta*A*sB) - 
                                             pnorm((delta*A*B-b)/sB) )
    if (k==2) M[2] <- B*( delta*A*M[1] + A*dnorm(delta*sB)*
                        ( dnorm(delta*A*sB)-dnorm((delta*A*B-b)/sB)) )
    if (k>2)  M[k] <- ((k-2)/(k-1))*B*(av[k-1]*delta*A*M[k-1] + M[k-2]) - L[k-2]
  }
  
  sumt <- 0 
  if (2*(nu%/%2)!= nu){
    # odd values of nu
    # sum over odd indices 1, 3, ..., nu-2
    # they are stored in k+1 (0 as index is not allowed), 
    # i.e. at the even indices
    if (upr>=1){
      # for (k in seq(1, upr, by=2)) sumt <- sumt + M[k+1] + H[k+1]
      # vectorized form of a for loop
      k    <- seq(1, upr, by=2)
      sumt <- sum(M[k+1]) + sum(H[k+1])
    }
    # if b==Inf what then?
    # rewrite second argument in first OwensT call (A*b-delta)/b which 
    # gives NaN to A-delta/b which gives A if b==Inf
      qv <- pnorm(b) - 2*OwensT(b, A-delta/b) - 
                       2*OwensT(delta*sB, (delta*A*B-b)/B/delta) +
                       2*OwensT(delta*sB, A) - (delta>=0) + 2*sumt
  } else {
    # even values of nu
    # sum over even indices 0, 2, ..., nu-2
    # they are stored in k+1
    if (upr>=0){
      # vectorized form of the for loop
      k    <- seq(0, upr, by=2)
      sumt <- sum(M[k+1]) + sum(H[k+1])
    }
    qv <- pnorm(-delta) + sqrt(2*pi)*sumt
  }
  return(qv)
}

# ----------------------------------------------------------------------------
# raw power function using OwensQOwen
.power.TOST.Q0 <- function(alpha=0.05, ltheta1, ltheta2, diffm, se, n, df, bk=2)
{
  tval   <- qt(1 - alpha, df, lower.tail = TRUE)
  # if alpha>0.5 (very unusual) then b=R is negative if not using abs()
  # in the application of OwensQ the upper integration limit 
  # is lower then the lower integration limit!
  # SAS OwenQ gives missings if b or a are negative!
  
  delta1 <- (diffm-ltheta1)/(se*sqrt(bk/n))
  delta2 <- (diffm-ltheta2)/(se*sqrt(bk/n))
  # 0/0 -> NaN in case diffm=ltheta1 or diffm=ltheta2 and se=0!
  delta1[is.nan(delta1)] <- 0
  delta2[is.nan(delta2)] <- 0
    
  # R is infinite in case of alpha>=0.5
  R <- (delta1-delta2)*sqrt(df)/(2.*tval)
  R[R<0] <- Inf
  # to avoid numerical errors in OwensQ implementation
  if (min(df)>10000) { 
    # Joulious formula (57) or (67), normal approximation
    p1 <- pnorm( (abs(delta1)-tval), lower.tail = TRUE)
    p2 <- pnorm( (abs(delta2)-tval), lower.tail = TRUE)
    
    return(p1 + p2 -1.)
  }
  # attempt to vectorize (it vectorizes properly if diffm is a vector
  # OR se OR n,df are vectors) 
  nel <- length(delta1)
  dl <- length(tval)
  p1 <- c(1:nel)	
  p2 <- p1
#  print(df)
#  print(tval)
#  print(delta1)
#  print(delta2)
#  print(R)
  for (i in seq_along(delta1)) {
    if (dl>1) {ddf <- df[i]; ttt <- tval[i]} 
        else {ddf <- df[1]; ttt <- tval[1]}
    p1[i] <- OwensQOwen(ddf,  ttt, delta1[i], 0, R[i])
    p2[i] <- OwensQOwen(ddf, -ttt, delta2[i], 0, R[i])
  }
  return( p2-p1 )
}

# check the function(s)
#require(PowerTOST)
#se <- CV2se(0.075)
#.power.TOST.Q0(alpha=0.05, ltheta1=log(0.8), ltheta2=log(1.25), 
#               diffm=0, se=se, n=4, df=2, bk=2)
#PowerTOST:::.power.TOST(alpha=0.05, ltheta1=log(0.8), ltheta2=log(1.25), 
#                        diffm=0, se=se, n=4, df=2, bk=2)
#.power.TOST.Q0(alpha=0.05, ltheta1=log(0.8), ltheta2=log(1.25), 
#               diffm=0, se=se, n=6, df=4, bk=2)
#PowerTOST:::.power.TOST(alpha=0.05, ltheta1=log(0.8), ltheta2=log(1.25), 
#                        diffm=0, se=se, n=6, df=4, bk=2)
#           
#
#OwensQ(2,-2.919986,-4.213542,0,2.040712)
#OwensQ(2,2.919986,4.213542,0,2.040712)
#
#OwensQOwen(2, 2.919986, 4.213542,0,2.040712)
