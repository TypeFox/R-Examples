
###############################################################################
# Exact MISE for normal mixtures
###############################################################################

## nu, gamma.r, gamma.r2 written by Jose Chacon 10/2008

nu <- function(r, A)
{ ###Using the recursive formula provided in Kan (2008)
  ##if (!inverse) A <- solve(A)
  ei <- eigen(A)$values
  tr.vec <- numeric(r)
  for(p in 1:r)
    tr.vec[p] <- sum(ei^p)

  nu.val <- 1
  if (r>=1)
  {
    for(p in 1:r)
    {
      a <- sum(tr.vec[1:p]*rev(nu.val))/(2*p)
      nu.val <- c(nu.val,a)
    }
  } 
  return(factorial(r)*2^r*nu.val[r+1])   
}   

nu.rs <- function(r, s, A, B)
{
  if (s==0)
    return(nu(r=r, A=A))
  
  if (s >=1)
  {
    nu.val <- 0
    for (i in 0:r)
      for (j in 0:(s-1))
        nu.val <- nu.val + choose(r,i)*choose(s-1,j) *factorial(r+s-i-j-1)*2^(r+s-i-j-1)*tr(matrix.pow(A,r-i)%*%matrix.pow(B,s-j))*nu.rs(r=i,s=j, A=A, B=B)
       
    return(nu.val)
  }
}


## gamma functional for normal mixture MISE
gamma.r <- function(mu, Sigma, r)
{
  Sigmainv <- chol2inv(chol(Sigma))
  d <- ncol(Sigma)
  v <- 0
  for (j in 0:r)
    v <- v + (-1)^j*choose(2*r, 2*j)*OF(2*j)*nu.rs(r=r-j, s=j, A=Sigmainv%*%mu%*%t(mu)%*%Sigmainv, B=Sigmainv)
   
  v <- (-1)^r*v*drop(dmvnorm.deriv(x=rep(0,d), mu=mu, Sigma=Sigma,deriv.order=0)/OF(2*r))
  return(v)
}


## gamma functional for normal mixture AMISE 
gamma.r2 <- function(mu, Sigma, d, r, H)
{
  Sigmainv <- chol2inv(chol(Sigma))

  if (d==1)
    w <- vec(Kpow(Sigmainv %*% Sigmainv, r)) %x% vec(Kpow(Sigmainv %*% H %*% Sigmainv, 2))
  else
    w <- matrix(Sdrv(d=d,r=2*r+4, v=vec(Kpow(Sigmainv %*% Sigmainv, r)) %x% vec(Kpow(Sigmainv %*% H %*% Sigmainv, 2))), nrow=1)
  
  v <- rep(0,length=d^(2*r+4))
  for(j in 0:(r+2))
    v <- v+((-1)^j*OF(2*j)*choose(2*r+4, 2*j))*(Kpow(mu,2*r-2*j+4)%x%Kpow(vec(Sigma),j))
  
  gamr<-(-1)^r*dmvnorm(mu,mean=rep(0,d),sigma=Sigma)*sum(w %*% v)

  return(gamr)
}


###############################################################################
# Omega matrices (for exact MISE for normal mixtures)
#
# Parameters 
# mus - means
# Sigmas - variances
# k - number of mixture components
# a - subscript of Omega matrix
# H - bandwidth matrix
#
# Returns 
# Omega matrix
###############################################################################

omega <- function(mus, Sigmas, k, a, H, r)
{
  ## the (i,j) element of Omega matrix is
  ## dmvnorm(0, mu_i - mu_j, a*H + Sigma_i + Sigma_j)
  if (is.matrix(mus)) d <- ncol(mus)
    else d <- length(mus)

  if (k == 1)
    omega.mat <- gamma.r(mu=rep(0,d),Sigma=a*H + 2*Sigmas, r=r)  
  else
  {
    omega.mat <- matrix(0, nrow=k, ncol=k)
    for (i in 1:k)
    {
      Sigmai <- Sigmas[((i-1)*d+1):(i*d),]
      mui <- mus[i,]
      for (j in 1:k)
      {
        Sigmaj <- Sigmas[((j-1)*d+1):(j*d),]
        muj <- mus[j,]    
        omega.mat[i,j] <- gamma.r(mu=mui-muj, Sigma=a*H + Sigmai + Sigmaj, r=r) 
      }
    }
  }
  
  return(omega.mat)
}

omega.1d <- function(mus, sigmas, k, a, h, r)
{
  H <- h^2
  Sigmas <- sigmas^2
  
  if (k == 1)
    omega.mat <- gamma.r(mu=0, Sigma=as.matrix(a*H + 2*Sigmas), r=r)  
  else
  {   
    omega.mat <- matrix(0, nrow=k, ncol=k)
    for (i in 1:k)
    {
      Sigmai <- Sigmas[i]
      mui <- mus[i]
      for (j in 1:k)
      {
        Sigmaj <- Sigmas[j]
        muj <- mus[j]    
        omega.mat[i,j] <- gamma.r(mu=mui-muj, Sigma=as.matrix(a*H + Sigmai + Sigmaj),  r=r) 
      }
    }
  }
  
  return(omega.mat)
}


##############################################################################
# Exact MISE for normal mixtures
#
# Parameters
# mus - means
# Sigmas - variances
# Props - vector of proportions of each mixture component 
# H - bandwidth matrix
# samp - sample size
#
# Returns
# Exact MISE for normal mixtures
###############################################################################

mise.mixt <- function(H, mus, Sigmas, props, samp, h, sigmas, deriv.order=0)
{
  if (!(missing(h)))
    return(mise.mixt.1d(h=h, mus=mus, sigmas=sigmas, props=props, samp=samp, deriv.order=deriv.order)) 

  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus)
  k <- length(props)
  r <- deriv.order

  Hinv <- chol2inv(chol(H))
  ## formula is found in Wand & Jones (1993) and Chacon, Duong & Wand (2008)
  if (k == 1) 
  {
    mise <- 2^(-r)*nu(r,Hinv)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) + 
      (1-1/samp)*omega(mus, Sigmas, 1, 2, H, r) - 2*omega(mus, Sigmas, 1, 1, H, r) + omega(mus, Sigmas, 1, 0, H, r)
  }
  else
  {
    mise <- 2^(-r)*nu(r,Hinv)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) +
      props %*% ((1-1/samp)*omega(mus, Sigmas, k, 2, H, r) - 2*omega(mus, Sigmas, k, 1, H, r) +  omega(mus, Sigmas, k, 0, H, r)) %*% props
  }
  return(drop(mise)) 
}

mise.mixt.1d <- function(h, mus, sigmas, props, samp, deriv.order=0)
{  
  d <- 1
  k <- length(props)
  r <- deriv.order
  H <- as.matrix(h^2)
  Hinv <- chol2inv(chol(H))
  ## formula is found in Wand & Jones (1993) and Chacon, Duong & Wand (2008)
  if (k == 1) 
  {
    mise <- 2^(-r)*nu(r,Hinv)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) + 
        (1-1/samp)*omega.1d(mus, sigmas, 1, 2, h, r) -
        2*omega.1d(mus, sigmas, 1, 1, h, r) +
          omega.1d(mus, sigmas, 1, 0, h, r)
  }
  else
  {
    mise <- 2^(-r)*nu(r,Hinv)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) +
      props %*% ((1-1/samp)*omega.1d(mus, sigmas, k, 2, h, r) - 
                 2*omega.1d(mus, sigmas, k, 1, h, r) + 
                 omega.1d(mus, sigmas, k, 0, h, r)) %*% props
  }
  return(drop(mise)) 
}



 
###############################################################################
# Exact AMISE for bivariate normal mixtures
#
# Parameters
# mus - means
# Sigmas - variances
# props - mixing proportions 
# H - bandwidth matrix
# samp - sample size
#
# Returns   
# Exact AMISE for normal mixtures
###############################################################################

amise.mixt <- function(H, mus, Sigmas, props, samp, h, sigmas, deriv.order=0)
{
  if (!(missing(h)))
    return(amise.mixt.1d(h=h, mus=mus, sigmas=sigmas, props=props, samp=samp, deriv.order=deriv.order))

  r <- deriv.order
  if (is.vector(mus)) {d <- length(mus); mus <- t(matrix(mus))}
  else d <- ncol(mus)
  k <- length(props)
 
  if (k == 1)
  {
    Sigmasinv <- chol2inv(chol(Sigmas))
    omega.mat <- 2^(-d-r-2)*pi^(-d/2)*det(Sigmas)^(-1/2)*nu.rs(r=r, s=2, Sigmasinv, matrix.sqrt(Sigmasinv) %*% H %*% matrix.sqrt(Sigmasinv))
  }
  else
  {
    omega.mat <- matrix(0, nrow=k, ncol=k)
    for (i in 1:k)
    {
      Sigmai <- Sigmas[((i-1)*d+1):(i*d),]
      mui <- mus[i,]
      for (j in 1:k)
      {
        Sigmaj <- Sigmas[((j-1)*d+1):(j*d),]
        muj <- mus[j,]    
        omega.mat[i,j] <- gamma.r2(mu=mui-muj, Sigma= Sigmai + Sigmaj, d=d, r=r, H=H)
      }
    }
  }

  Hinv <- chol2inv(chol(H)) 
  if (k == 1)amise <- 2^(-r)*nu(r,Hinv)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) + omega.mat/4
  else amise <- 2^(-r)*nu(r,Hinv)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) + (props %*% omega.mat %*% props)/4
 
  return(drop(amise))
}

amise.mixt.1d <- function(h, mus, sigmas, props, samp, deriv.order=0)
{
  d <- 1
  r <- deriv.order
  k <- length(props)
  H <- as.matrix(h^2)
  
  if (k == 1)
    omega.mat <- gamma.r2(mu=rep(0,d),Sigma=as.matrix(2*sigmas^2), d=d, r=r, H=H)
  else
  {  
    omega.mat <- matrix(0, nrow=k, ncol=k)
    for (i in 1:k)
    {
      Sigmai <- as.matrix(sigmas[i]^2)
      mui <- mus[i]
      for (j in 1:k)
      {
        Sigmaj <- as.matrix(sigmas[j]^2)
        muj <- mus[j]    
        omega.mat[i,j] <- gamma.r2(mu=mui-muj, Sigma= Sigmai + Sigmaj, d=d, r=r, H=H)
      }
    }
  }

  Hinv <- h^(-2)
  
  if (k == 1) amise <- 2^(-r)*nu(r,Hinv)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) + omega.mat/4
  else amise <- 2^(-r)*nu(r,Hinv)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) + (props %*% omega.mat %*% props)/4
  
  return(drop(amise))
}



###############################################################################
# Lambda matrices (for exact AMISE for normal mixtures)
#
# Parameters 
# mus - means
# Sigmas - variances
# k - number of mixture components
# r - derivative (r1, r2)
#
# Returns 
# Lambda matrix
###############################################################################

lambda <- function(mus, Sigmas, k, r)
{
  ## the (i,j) element of Lambda matrix is d^r/ dx^r  dmvnorm(0, mu_i - mu_j,
  ## a*H + Sigma_i + Sigma_j)
    
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus)
  
  if (k == 1) 
    lambda.mat <- dmvnorm.deriv(deriv.order=r, x=rep(0, length(mus)), Sigma=2*Sigmas)
  else
  {   
    if (is.matrix(mus)) d <- ncol(mus)
    else d <- length(mus)
    lambda.mat <- matrix(0, nrow=k, ncol=k)
    for (i in 1:k)
    {
      Sigmai <- Sigmas[((i-1)*d+1) : (i*d),]
      mui <- mus[i,]
      for (j in 1:k)
      {
        Sigmaj <- Sigmas[((j-1)*d+1) : (j*d),]
        muj <- mus[j,]    
        lambda.mat[i,j] <- dmvnorm.deriv(deriv.order=r, x=mui-muj,Sigma=Sigmai+Sigmaj)
      }
    }
  }
  
  return(lambda.mat)
}


amise.mixt.2d <- function(H, mus, Sigmas, props, samp)
{  
  d <- ncol(Sigmas)
  k <- length(props)

  ## formula is found in Wand & Jones (1993)
  if (k == 1) 
  {
   amise <- 1/(samp * (4*pi)^(d/2) * sqrt(det(H))) +
       1/4 *(lambda(mus, Sigmas, k, r=c(4,0))*H[1,1]^2 +
           4*lambda(mus, Sigmas, k, r=c(3,1))*H[1,1]*H[1,2] +  
           2*lambda(mus, Sigmas, k, r=c(2,2))*(H[1,1]*H[2,2] + 2*H[1,2]^2) + 
           4*lambda(mus, Sigmas, k, r=c(1,3))*H[2,2]*H[1,2]+    
             lambda(mus, Sigmas, k, r=c(0,4))*H[2,2]^2) 
  }
  else
  {
    amise <- 1/(samp * (4*pi)^(d/2) * sqrt(det(H))) +
      1/4 * props %*% 
          (  lambda(mus, Sigmas, k, r=c(4,0))*H[1,1]^2 +
           4*lambda(mus, Sigmas, k, r=c(3,1))*H[1,1]*H[1,2] +  
           2*lambda(mus, Sigmas, k, r=c(2,2))*(H[1,1]*H[2,2] + 2*H[1,2]^2) + 
           4*lambda(mus, Sigmas, k, r=c(1,3))*H[2,2]*H[1,2]+    
             lambda(mus, Sigmas, k, r=c(0,4))*H[2,2]^2) %*% props
  }
  
  return(drop(amise)) 
}


###############################################################################
# Finds the bandwidth matrix that minimises the MISE for normal mixtures
#
# Parameters
# mus - means
# Sigmas - variances
# props - vector of proportions of each mixture component 
# Hstart - initial bandwidth matrix
# samp - sample size
# full - 1 minimise over full bandwidth matrices
#      - 0 minimise over diagonal bandwidth matrices
# 
# Returns
# H_MISE
###############################################################################

hmise.mixt <- function(mus, sigmas, props, samp, hstart, deriv.order=0)
{
  r <- deriv.order
  d <- 1
  
  if (missing(hstart))
  {
    x <- rnorm.mixt(n=1000, mus=mus, sigmas=sigmas)  
    hstart <- sqrt((4/(samp*(d+2*r+2)))^(2/(d+2*r+4)) * var(x))
  }
  mise.mixt.temp <- function(h)
  {  
    return(mise.mixt.1d(h=h, mus=mus, sigmas=sigmas, props=props, samp=samp, deriv.order=deriv.order))
  }

  result <- optimize(f=mise.mixt.temp, interval=c(0, 10*hstart))
  hmise <- result$minimum
  return(hmise)
}

Hmise.mixt <- function(mus, Sigmas, props, samp, Hstart, deriv.order=0)
{
  r <- deriv.order
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus) 
 
  ## use normal reference estimate as initial condition
  if (missing(Hstart))
  {
    x <- rmvnorm.mixt(10000, mus, Sigmas, props)
    Hstart <- (4/(samp*(d+2*r+2)))^(2/(d+2*r+4)) * var(x)
  }
  
  mise.mixt.temp <- function(vechH)
  {  
    H <- invvech(vechH) %*% invvech(vechH)
    return(mise.mixt(H=H, mus=mus, Sigmas=Sigmas, props=props, samp=samp, deriv.order=deriv.order))
  }

  Hstart <- vech(matrix.sqrt(Hstart))
  result <- nlm(p=Hstart, f=mise.mixt.temp)
  Hmise <- invvech(result$estimate) %*% invvech(result$estimate)
  ##result <- optim(Hstart, mise.mixt.temp, method="Nelder-Mead")
  ##Hmise <- invvech(result$par) %*% invvech(result$par) 
  
  return(Hmise)
}   

Hmise.mixt.diag <- function(mus, Sigmas, props, samp, Hstart, deriv.order=0)
{   
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus) 

  if (missing(Hstart))
  {
    x <- rmvnorm.mixt(10000, mus, Sigmas, props)
    Hstart <- (4/(samp*(d + 2)))^(2/(d + 4)) * var(x)
  }  

  mise.mixt.temp <- function(diagH)
  {  
    H <- diag(diagH) %*% diag(diagH)
    return(mise.mixt(H=H, mus=mus, Sigmas=Sigmas, props=props, samp=samp, deriv.order=deriv.order))
  }
  Hstart <- diag(matrix.sqrt(Hstart))
  result <- nlm(p=Hstart, f=mise.mixt.temp)
  Hmise <- diag(result$estimate) %*% diag(result$estimate)
  ##result <- optim(Hstart, mise.mixt.temp, method = "Nelder-Mead")
  ##Hmise <- diag(result$par) %*% diag(result$par) 
  
  return(Hmise)
}   



###############################################################################
## Finds bandwidth matrix that minimises the AMISE for normal mixtures - 2-dim
##
## Parameters
## mus - means
## Sigmas - variances
## props - vector of proportions of each mixture component 
## Hstart - initial bandwidth matrix
## samp - sample size
## 
## Returns
## Bandwidth matrix that minimises AMISE
###############################################################################

hamise.mixt <- function(mus, sigmas, props, samp, hstart, deriv.order=0)
{
  r <- deriv.order
  d <- 1
 
  if (missing(hstart))
  {
    x <- rnorm.mixt(n=1000, mus=mus, sigmas=sigmas)  
    hstart <- sqrt((4/(samp*(d+2*r+2)))^(2/(d+2*r+4)) * var(x))
  }
  amise.mixt.temp <- function(h)
  {  
    return(amise.mixt.1d(h=h, mus=mus, sigmas=sigmas, props=props, samp=samp, deriv.order=deriv.order))
  }

  result <- optimize(f=amise.mixt.temp, interval=c(0, 10*hstart))
  hamise <- result$minimum
  
  return(hamise)
}


Hamise.mixt <- function(mus, Sigmas, props, samp, Hstart, deriv.order=0)
{
  r <- deriv.order
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus) 

  ## use explicit formula for single normal
  if (length(props)==1)
  {
    Hamise <- (4/ (samp*(d+2*r+2)))^(2/(d+2*r+4)) * Sigmas
  }
  else
  {  
    ## use normal reference estimate as initial condition
    if (missing(Hstart)) 
    {
      x <- rmvnorm.mixt(10000, mus=mus, Sigmas=Sigmas, props=props)
      Hstart <- matrix.sqrt((4/ (samp*(d+2*r+2)))^(2/(d+2*r+4)) * var(x))
    }
    
    amise.mixt.temp <- function(vechH)
    {
      H <- invvech(vechH) %*% invvech(vechH)
      return(amise.mixt(H=H, mus=mus, Sigmas=Sigmas, props=props, samp=samp, deriv.order=deriv.order))
    }

    result <- nlm(p=vech(Hstart), f=amise.mixt.temp)
    Hamise <- invvech(result$estimate) %*% invvech(result$estimate)
  }
  
  return(Hamise)
}   

Hamise.mixt.diag <- function(mus, Sigmas, props, samp, Hstart, deriv.order=0)
{
  r <- deriv.order
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus) 

  ## use normal reference estimate as initial condition
  if (missing(Hstart)) 
  {
    x <- rmvnorm.mixt(10000, mus, Sigmas, props)
    Hstart <- matrix.sqrt((4/ (samp*(d+2*r+2)))^(2/(d+2*r+4)) * var(x))
  }
  
  amise.mixt.temp <- function(diagH)
  {
    H <- diag(diagH) %*% diag(diagH)
    return(amise.mixt(H=H, mus=mus, Sigmas=Sigmas, props=props, samp=samp, deriv.order=deriv.order))
  }
    
  result <- nlm(p=diag(Hstart), f=amise.mixt.temp)
  Hamise <- diag(result$estimate) %*% diag(result$estimate)
  
  return(Hamise)
}   
  


###############################################################################
# ISE for normal mixtures (fixed KDE)
# 
# Parameters
# x - data values
# H - bandwidth matrix
# mus - matrix of means (each row is a vector of means from each component
#       density)
# Sigmas - matrix of covariance matrices (every d rows is a covariance matrix 
#          from each component density) 
# props - mixing proportions
# lower - vector of lower end points of rectangle
# upper - vector of upper end points of rectangle
# gridsize - vector of number of grid points
# stepsize - vector of step sizes
# Returns
# ISE 
###############################################################################

ise.mixt <- function(x, H, mus, Sigmas, props, h, sigmas, deriv.order=0, binned=FALSE, bgridsize)
{
  if (!(missing(h)))
    return(ise.mixt.1d(x=x, h=h, mus=mus, sigmas=sigmas, props=props, deriv.order=deriv.order, binned=binned))
  
  if (is.vector(x)) x <- matrix(x, nrow=1)
  if (is.vector(mus)) mus <- matrix(mus, nrow=length(props))

  d <- ncol(x)
  n <- nrow(x)
  M <- length(props)
  r <- deriv.order
 
  ## formula is found in thesis

  vIdr <- vec(diag(d^r))
  ise1 <- 0
  ise2 <- 0
  ise3 <- 0
  
  if (binned)
  {
    ise1 <- dmvnorm.deriv.sum(x=x, Sigma=2*H, inc=1, deriv.order=2*r, binned=binned, bgridsize=bgridsize)
    for (j in 1:M)
    {
      Sigmaj <- Sigmas[((j - 1) * d + 1):(j * d), ]
      ise2 <- ise2 + props[j] * colSums(dmvnorm.deriv(x, mu=mus[j,],Sigma=H + Sigmaj, deriv.order=2*r))
      for (i in 1:M)
      {
        Sigmai <- Sigmas[((i - 1) * d + 1):(i * d), ]
        ise3 <- ise3 + props[i] * props[j] * dmvnorm.deriv(x=mus[i,],mu=mus[j,], Sigma = Sigmai + Sigmaj, deriv.order = 2*r)
      }
    }
    ise <- (-1)^r * sum(vIdr*(ise1/n^2 - 2 * ise2/n + ise3))
  }
  else
  {
    ise1 <- Qr(x=x, Sigma=2*H, inc=1, deriv.order=2*r)
    for (j in 1:M)
    {
      Sigmaj <- Sigmas[((j-1)*d + 1):(j*d), ]
      ise2 <- ise2 + props[j] * Qr(x=x, y=mus[j,], Sigma=H + Sigmaj, deriv.order=2*r, inc=1)
      for (i in 1:M)
      {
        Sigmai <- Sigmas[((i-1)*d + 1):(i*d), ]
        ise3 <- ise3 + props[i] * props[j] * Qr(x=mus[i,], y=mus[j,], Sigma=Sigmai + Sigmaj, deriv.order=2*r, inc=1)
      }
    }
    ise <- (-1)^r*(ise1 - 2*ise2 + ise3)
  }
 
  return(ise)
}

ise.mixt.1d <- function(x, h, mus, sigmas, props, deriv.order=0, binned=FALSE)
{  
  n <- length(x)
  M <- length(props)
  r <- deriv.order 
  ise1 <- 0
  ise2 <- 0
  ise3 <- 0
  
  ise1 <- dnorm.deriv.sum(x=x, sigma=sqrt(2)*h, inc=1, deriv.order=2*r, binned=binned)
  
  for (j in 1:M)
  {
    sigmaj <- sigmas[j]
    ise2 <- ise2 + sum(props[j]*dnorm.deriv(x=x, mu=mus[j], sigma=sqrt(h^2 + sigmaj^2), deriv.order=2*r))
    
    for (i in 1:M)
    {
      sigmai <- sigmas[i]
      ise3 <- ise3 + sum(props[i]*props[j]*dnorm.deriv(x=mus[i], mu=mus[j], sigma=sqrt(sigmai^2+sigmaj^2), deriv.order=2*r))
    }
  }  

  return ((-1)^r*(ise1/n^2 - 2*ise2/n + ise3))
}


Hise.mixt <- function(x, mus, Sigmas, props, Hstart, deriv.order=0)
{
  r <- deriv.order
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus) 
  samp <- nrow(x)
  
  ## use normal reference estimate as initial condition
  if (missing(Hstart))
  {
    xstart <- rmvnorm.mixt(10000, mus, Sigmas, props)
    Hstart <- (4/(samp*(d+2*r+2)))^(2/(d+2*r+4)) * var(xstart)
  }
  
  ise.mixt.temp <- function(vechH)
  {  
    H <- invvech(vechH) %*% invvech(vechH)
    return(ise.mixt(x=x, H=H, mus=mus, Sigmas=Sigmas, props=props, deriv.order=deriv.order))
  }

  Hstart <- vech(matrix.sqrt(Hstart))
  result <- nlm(p=Hstart, f=ise.mixt.temp)
  Hise <- invvech(result$estimate) %*% invvech(result$estimate)
  return(Hise)
}   

