
#############################################################################
## Kernel functional estimation
#############################################################################

kfe <- function(x, G, deriv.order, inc=1, binned=FALSE, bin.par, bgridsize, deriv.vec=TRUE, add.index=TRUE, verbose=FALSE)
{
  r <- deriv.order
  d <- ncol(x)
  psir <- dmvnorm.deriv.sum(x=x, Sigma=G, deriv.order=r, inc=inc, binned=binned, bin.par=bin.par, bgridsize=bgridsize, deriv.vec=deriv.vec, verbose=verbose, kfe=TRUE, add.index=FALSE)
 
 if (add.index)
 {
    ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=r, only.index=TRUE)
  
    if (deriv.vec) return(list(psir=psir, deriv.ind=ind.mat))
    else return(list(psir=psir, deriv.ind=unique(ind.mat)))
  }
  else return(psir=psir)
}


kfe.1d <- function(x, g, deriv.order, inc=1, binned=FALSE, bin.par)
{
  r <- deriv.order
  n <- length(x)
  psir <- dnorm.deriv.sum(x=x, sigma=g, deriv.order=r, inc=1, binned=binned, bin.par=bin.par, kfe=TRUE)
  if (inc==0)  psir <- (n^2*psir - n*dnorm.deriv(0, mu=0, sigma=g, deriv.order=r))/(n*(n-1))
  
  return(psir) 
}


kfe.scalar <- function(x, g, deriv.order, inc=1, binned=FALSE, bin.par, verbose=FALSE)
{
  r <- deriv.order
  d <- ncol(x)
  ##if (missing(bin.par) & binned) bin.par <- binning(x=x, H=g^2*diag(d))
  
  psir <- dmvnorm.deriv.scalar.sum(x=x, sigma=g, deriv.order=r, inc=inc, kfe=TRUE, binned=binned, bin.par=bin.par, verbose=verbose)
  return(psir)
}


###############################################################################
## Plug-in unconstrained bandwidth for KFE
##
## Returns
## Plug-in bandwidth
###############################################################################


hpi.kfe <- function(x, nstage=2, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0)
{
  n <- length(x)
  d <- 1
  r <- deriv.order
  k <- 2 ## kernel order
  Kr0 <- dnorm.deriv(x=0, mu=0, sigma=1, deriv.order=r)    
  mu2K <- 1  
  
  if (nstage==2)
  {
    psi4.hat <- psins.1d(r=r+k+2, sigma=sd(x))
    gamse2 <- (factorial(r+2)*Kr0/(mu2K*psi4.hat*n))^(1/(r+k+3))
    psi2.hat <- kfe.1d(x=x, g=gamse2, deriv.order=r+k, inc=1, binned=binned)
  }
  else 
  {
    psi2.hat <- psins.1d(r=r+k, sigma=sd(x))
  }

  ## formula for bias annihliating bandwidths from Wand & Jones (1995, p.70)
  gamse <- (factorial(r)*Kr0/(-mu2K*psi2.hat*n))^(1/(r+k+1)) 
  
  return(gamse)
}


Hpi.kfe <- function(x, nstage=2, pilot, pre="sphere", Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm")
{
  if (deriv.order!=0) stop("Currently only deriv.order=0 is implemented")
   
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  
 if (missing(pilot)) pilot <- "dscalar"
  if(!is.matrix(x)) x <- as.matrix(x)
  pilot1 <- match.arg(pilot, c("dunconstr", "dscalar"))  
  pre1 <- match.arg(pre, c("scale", "sphere"))
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))
  
  if (pre1=="scale")
  {
    x.star <- pre.scale(x)
    S12 <- diag(sqrt(diag(var(x))))
    Sinv12 <- chol2inv(chol(S12))
  }
  else if (pre1=="sphere")
  {
    x.star <- pre.sphere(x)
    S12 <- matrix.sqrt(var(x))
    Sinv12 <- chol2inv(chol(S12))
  }

  D2K0 <- t(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=2))
  K0 <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=0)

  if (nstage==2)
  {  
    ## stage 1
    if (pilot1=="dscalar")
    {
      psi4.ns <- psins(r=r+4, Sigma=var(x.star), deriv.vec=TRUE)
      ##h2 <- gdscalar(x=x.star, d=d, r=r-2, n=n, verbose=verbose, nstage=nstage, scv=FALSE)     
      ## gdscalar not used because it's an implicit computation without
      ## symmetriser matrices and is slower than direct computation with 
      ## symmetriser matrices. 

      A1 <- drop(t(D2K0) %*% D2K0)
      A2 <- drop(t(D2K0) %*% t(vec(diag(d)) %x% diag(d^2)) %*% psi4.ns) 
      A3 <- drop(t(psi4.ns) %*% (vec(diag(d)) %*% t(vec(diag(d))) %x% diag(d^2)) %*% psi4.ns)

      ## Special case from Chacon & Duong (2009): bias minimisation
      ##h2 <- ((4*d+8)*A1/(-d*A2 + sqrt(d^2*A2^2 + (8*d+16)*A1*A3)))^(1/(d+4))*n^(-1/(d+4))
      ## Special from Chacon & Duong (2009): bias annihilation
      h2 <- (-A1/(2*A2*n))^(1/(d+4))
      H2 <- h2^2*diag(d)   
      psi2.hat <- kfe(x=x.star, G=H2, deriv.order=r+2, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
    }
    else if (pilot1=="dunconstr")
    {
      psi4.ns <- psins(r=r+4, Sigma=var(x), deriv.vec=TRUE)
 
      amse2.temp <- function(vechH)
      { 
        H <- invvech(vechH) %*% invvech(vechH)
        Hinv <- chol2inv(chol(H))
        Hinv12 <- matrix.sqrt(Hinv)
        amse2.temp <- 1/(det(H)^(1/2)*n)*((Hinv12 %x% Hinv12) %*% D2K0) + 1/2* t(vec(H) %x% diag(d^2)) %*% psi4.ns
        return(sum((amse2.temp)^2)) 
      }
      
      Hstart2 <- matrix.sqrt(Gns(r=r+2, n=n, Sigma=var(x)))
     
 
      if (optim.fun1=="nlm")
      {
         result <- nlm(p=vech(Hstart2), f=amse2.temp, print.level=2*as.numeric(verbose))    
         H2 <- invvech(result$estimate) %*% invvech(result$estimate)
      }
      else if (optim.fun1=="optim")
      {
         result <- optim(vech(Hstart2), amse2.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
         H2 <- invvech(result$par) %*% invvech(result$par)
      }
 
      psi2.hat <- kfe(x=x, G=H2, deriv.order=r+2, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
    }
  }
  else
  {
    if (pilot1=="dscalar") psi2.hat <- psins(r=r+2, Sigma=var(x.star), deriv.vec=TRUE)
    else if (pilot1=="dunconstr") psi2.hat <- psins(r=r+2, Sigma=var(x), deriv.vec=TRUE)    
  }

  if (pilot1=="dscalar") {if (missing(Hstart)) Hstart <- Gns(r=r, n=n, Sigma=var(x.star))}
  else if (pilot1=="dunconstr") {if (missing(Hstart)) Hstart <- Gns(r=r, n=n, Sigma=var(x))}
  
  ## stage 2
  amse.temp <- function(vechH)
  { 
    H <- invvech(vechH) %*% invvech(vechH)
    amse.val <- 1/(det(H)^(1/2)*n)*K0 + 1/2* t(vec(H)) %*% psi2.hat
    return(sum((amse.val^2))) 
  }
  
  Hstart <- matrix.sqrt(Hstart)
  if (optim.fun1=="nlm")
  {
    result <- nlm(p=vech(Hstart), f=amse.temp, print.level=2*as.numeric(verbose)) 
    H <- invvech(result$estimate) %*% invvech(result$estimate)
    amise.star <- result$minimum
  }
  else if (optim.fun1=="optim")
  {
    result <- optim(vech(Hstart), amse.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- invvech(result$par) %*% invvech(result$par)
    amise.star <- result$value
  }
  
  ## back-transform
  if (pilot1=="dscalar") H <- S12 %*% H %*% S12 

  if (!amise) return(H)
  else return(list(H=H, PI=amise.star))
}

###############################################################################
## Plug-in diagonal bandwidth for KFE
##
## Returns
## Plug-in bandwidth
###############################################################################


 
Hpi.diag.kfe <- function(x, nstage=2, pilot, pre="scale", Hstart, binned=FALSE, bgridsize, amise=FALSE,  deriv.order=0, verbose=FALSE, optim.fun="nlm")
{
  if (deriv.order!=0) stop("Currently only dervi.order=0 is implemented")

  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  D2K0 <- t(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=2))
  K0 <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=0)

  if (missing(pilot)) pilot <- "dscalar"
  pilot1 <- match.arg(pilot, c("dunconstr", "dscalar"))  
  pre1 <- match.arg(pre, c("scale", "sphere"))
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))

  if (pre1=="sphere") stop("Using pre-sphering won't give diagonal bandwidth matrix")
  if (pilot1=="dunconstr") stop("Unconstrained pilot selectors are not suitable for Hpi.diag.kfe")
  
  if (pre1=="scale")
  {
    x.star <- pre.scale(x)
    S12 <- diag(sqrt(diag(var(x))))
    Sinv12 <- chol2inv(chol(S12))
  }

  if (nstage==2)
  {  
    ## stage 1
    psi4.ns <- psins(r=r+4, Sigma=var(x.star), deriv.vec=TRUE)
   
    if (pilot1=="dscalar")
    {
      ## h2 is on pre-transformed data scale
      A1 <- drop(t(D2K0) %*% D2K0)
      A2 <- drop(t(D2K0) %*% t(vec(diag(d)) %x% diag(d^2)) %*% psi4.ns) 
      A3 <- drop(t(psi4.ns) %*% (vec(diag(d)) %*% t(vec(diag(d))) %x% diag(d^2)) %*% psi4.ns)
      h2 <- ((4*d+8)*A1/(-d*A2 + sqrt(d^2*A2^2 + (8*d+16)*A1*A3)))^(1/(d+4))*n^(-1/(d+4))
      H2 <- h2^2*diag(d)
    }  
    psi2.hat <- kfe(x=x.star, G=H2, deriv.order=r+2, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose) 
  }
  else
    psi2.hat <- psins(r=r+2, Sigma=var(x.star), deriv.vec=TRUE)
 
  ## stage 2
  amse.temp <- function(diagH)
  { 
    H <- diag(diagH) %*% diag(diagH)
    amse.val <- 1/(det(H)^(1/2)*n)*K0 + 1/2* t(vec(H)) %*% psi2.hat
    return(sum((amse.val^2))) 
  }

  if (missing(Hstart)) Hstart <- Hns(x=x.star, deriv.order=r)
  Hstart <- matrix.sqrt(Hstart)
  
  if (optim.fun1=="nlm")
  {
    result <- nlm(p=diag(Hstart), f=amse.temp, print.level=2*as.numeric(verbose))    
    H <- diag(result$estimate) %*% diag(result$estimate)
    amise.star <- result$minimum
  }
  else if (optim.fun1=="optim")
  {
    result <- optim(diag(Hstart), amse.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- diag(result$par) %*% diag(result$par)
    amise.star <- result$value
  }
  ## back-transform
  if (pilot1=="dscalar") H <- S12 %*% H %*% S12 
  if (!amise) return(H)
  else return(list(H=H, PI=amise.star))
}

