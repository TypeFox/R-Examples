###############################################################################
## Estimate g_AMSE pilot bandwidths for even orders - 2-dim
##
## Parameters
## r - (r1, r2) partial derivative
## n - sample size
## psi1 - psi_(r + (2,0))
## psi2 - psi_(r + (0,2))
##
## Returns
## g_AMSE pilot bandwidths for even orders
###############################################################################

gamse.even.2d <- function(r, n, psi1, psi2)
{
  d <- 2
  num <- -2 * dmvnorm.deriv(x=c(0,0), deriv.order=r, Sigma=diag(2), deriv.vec=FALSE)
  den <- (psi1 + psi2) * n   
  g.amse <- (num/den)^(1/(2 + d + sum(r)))
  
  return(g.amse)
}

###############################################################################
## Estimate g_AMSE pilot bandwidths for odd orders - 2-dim
##
## Parameters
##
## r - (r1, r2) partial derivative
## n - sample size
## psi1 - psi_(r + (2,0))
## psi2 - psi_(r + (0,2))
## psi00 - psi_(0,0)
## RK - R(K^(r))
## 
## Returns
## g_AMSE pilot bandwidths for odd orders
###############################################################################

gamse.odd.2d <- function(r, n, psi1, psi2, psi00, RK)
{  
  d <- 2
  num <- 2 * psi00 * (2 * sum(r) + d) * RK
  den <- (psi1 + psi2)^2 * n^2
  g.amse <- (num/den)^(1/(2*sum(r) + d + 4))
  
  return(g.amse)
}


###############################################################################
## Estimate g_SAMSE pilot bandwidth - 2- to 6-dim 
##
## Parameters
## Sigma.star - scaled variance matrix
## n - sample size
##
## Returns
## g_SAMSE pilot bandwidth
###############################################################################

gsamse <- function(Sigma.star, n, modr, nstage=1, psihat=NULL)
{
  d <- ncol(Sigma.star)
  K <- numeric(); psi <- numeric()

  ## 4th order g_SAMSE

  K <- dmvnorm.deriv(x=rep(0,d), deriv.order=modr, Sigma=diag(d), add.index=TRUE, deriv.vec=FALSE)
  K <- K$deriv[apply(K$deriv.ind, 1, is.even)]
 
  if (modr==4)
  {
    derivt4 <- dmvnorm.deriv(x=rep(0,d), deriv.order=4, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
    derivt6 <- dmvnorm.deriv(x=rep(0,d), deriv.order=6, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
  
    for (i in 1:nrow(derivt4))
    {
      r <- derivt4[i,]
      if (is.even(r))
      {
        A3psi <- 0
        for (j in 1:d)
        {
          if (nstage==1)
          {
            A3psi <- A3psi + psins(r=r+2*elem(j,d), Sigma=Sigma.star)
            ##A3psi <- A3psi + psins6[which.mat(r=r+2*elem(j,d), mat=derivt6)]
          }
          else if (nstage==2)
            A3psi <- A3psi + psihat[which.mat(r=r+2*elem(j,d), mat=derivt6)]
        }
        psi <- c(psi, A3psi)    
      }
    }
  }
  ## 6th order g_SAMSE
  else if (modr==6)
  {
    derivt6 <- dmvnorm.deriv(x=rep(0,d), deriv.order=6, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
  
    for (i in 1:nrow(derivt6))
    {
      r <- derivt6[i,]        
      if (is.even(r))
      {
        A3psi <- 0
        for (j in 1:d)
          A3psi <- A3psi + psins(r=r+2*elem(j,d), Sigma=Sigma.star)
        psi <- c(psi, A3psi)
      }
    }
  }

  ## see thesis for formula
  A1 <- sum(K^2)
  A2 <- sum(K * psi)  
  A3 <- sum(psi^2)
  B1 <- (2*modr + 2*d)*A1
  B2 <- (modr + d - 2)*A2
  B3 <- A3
  gamma1 <- (-B2 + sqrt(B2^2 + 4*B1*B3)) / (2*B1)
  g.samse <- (gamma1 * n)^(-1/(modr + d + 2))
  
  return (g.samse)      
}

##############################################################################
## Scalar pilot selector for derivatives r>0 from Chacon & Duong (2011)
## Generalisation of gsamse for r>0
##############################################################################

gdscalar <- function(x, d, r, n, verbose, binned=FALSE, nstage=1, scv=FALSE)
{
  if (scv) cf <- c(2^(-d), 2^(-d/2+1), 4)
  else cf <- c(1,1,1)
  if (nstage==1)
  {
    G2r4 <- Gns(r=2*r+4,n=n,Sigma=var(x))
    g2r4 <- sqrt(G2r4[1,1])
  }
  else if (nstage==2)
  {
    G2r6.NR <- Gns(r=2*r+6,n=n,Sigma=var(x))
    g2r6.nr <- prod(sqrt(diag(G2r6.NR)))^(1/d) 
    L0 <- dmvnorm.mixt(x=rep(0,d), mus=rep(0,d), Sigmas=diag(d), props=1)
    if (binned)
      eta2r6 <- drop(kfe(x=x, G=g2r6.nr^2*diag(d), inc=1, binned=binned, deriv.order=2*r+6, add.index=FALSE, verbose=verbose) %*% vec(diag(d^(r+3))))
    else
      eta2r6 <- Qr(x=x, deriv.order=2*r+6, Sigma=g2r6.nr^2*diag(d), inc=1)

    A1 <- cf[1]*(2*d+4*r+8)*L0^2*OF(2*r+4)*nu(r=r+2, A=diag(d))
    A2 <- cf[2]*(-1)^(r+2)*(d+2*r+2)*L0*OF(2*r+4)*eta2r6
    A3 <- cf[3]*eta2r6^2
    
    g2r4 <- (2*A1/((-A2+ sqrt(A2^2 +4*A1*A3))*n))^(1/(d+2*r+6))
  }
  return(g2r4)
}

##############################################################################
## Unconstrained pilot selector for derivatives r>0 from Chacon & Duong (2011)
## Generalisation of Gunconstr for r>0
##############################################################################

Gdunconstr <- function(x, d, r, n, nstage=1, verbose, binned=FALSE, scv=FALSE, optim.fun="nlm")
{
  if (scv) cf <- c(2^(-d/2), 2)
  else cf <- c(1,1)
  S <- var(x)
  
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim")) 
  if (nstage==1)
  {
    G2r4 <- Gns(r=2*r+4,n=n,Sigma=S)
  }
  else if (nstage==2)
  {
    G2r4.NR <- Gns(r=2*r+4,n=n,Sigma=S)
    G2r6.NR <- Gns(r=2*r+6,n=n,Sigma=S)
    vecPsi2r6 <- kfe(x=x, G=G2r6.NR, binned=binned, deriv.order=2*r+6, deriv.vec=TRUE, add.index=FALSE, verbose=verbose)   
    
    dls <- (0:(d^2-1))*d^(2*r+4) 

    AB2 <- function(vechG)
    {
      G <- invvech(vechG) %*% invvech(vechG)
      Ginv <- chol2inv(chol(G))
    
      ## direct computation
      ##v1 <- n^(-1)*det(Ginv12)*drop(Kpow(A=Ginv12,pow=2*r+4)%*%D2r4phi0) 
      v1 <- n^(-1)*dmvnorm(rep(0,d),rep(0,d),sigma=G)*(-1)^(r+2)*OF(2*r+4)*Sdrv(d=d,r=2*r+4,v=Kpow(vec(Ginv),r+2))
      
      ##v2 <- (1/2)*drop((t(vec(G))%x%Id2r4) %*% vecPsi2r6)
      v2 <- numeric(d^(2*r+4))
        for(k in 1:d^(2*r+4)){v2[k]<-(1/2)*sum(vec(G)*vecPsi2r6[dls+k])}
      
      AB <- cf[1]*v1 + cf[2]*v2
      AB2.val <- sum(AB^2)
      return(AB2.val)
    }
    
    Gstart <- matrix.sqrt(G2r4.NR)
   
    if (optim.fun1=="nlm")
    {
      result <- nlm(p=vech(Gstart), f=AB2, print.level=2*as.logical(verbose))
      G2r4 <- result$estimate
    } 
    else if (optim.fun=="optim")
    { 
      result <- optim(vech(matrix.sqrt(Gstart)), AB2, method="BFGS", control=list(trace=as.numeric(verbose)))
      G2r4 <- result$par
    }   
    G2r4 <- invvech(G2r4)%*%invvech(G2r4)
  }
  return(G2r4)
}

##############################################################################
## Estimate psi functionals using 1-stage plug-in 
##
## Parameters
## x.star - pre-transformed data points
## pilot - "amse" = different AMSE pilot bandwidths
##       - "samse" = optimal SAMSE pilot bandwidth
##
## Returns
## estimated psi functionals
###############################################################################

psifun1 <- function(x.star, pilot="samse", binned, bin.par, deriv.order=0, verbose=FALSE)
{
  d <- ncol(x.star)
  r <- deriv.order
  S.star <- var(x.star)
  n <- nrow(x.star)
  
  ## pilots are based on (2r+4)-th order derivatives
  ## compute 1 pilot for SAMSE
  if (pilot=="samse")
  {
    g.star <- gsamse(S.star, n, 4)
    psihat.star <- kfe(x=x.star, G=g.star^2*diag(d), deriv.order=4, deriv.vec=TRUE, binned=binned, add.index=TRUE, verbose=verbose)
  }
  ## compute 5 different pilots for AMSE
  else if ((pilot=="amse") & (d==2))
  {
    derivt4 <- dmvnorm.deriv(x=rep(0,d), deriv.order=4, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
    derivt4.vec <- dmvnorm.deriv(x=rep(0,d), deriv.order=4, add.index=TRUE, deriv.vec=TRUE, only.index=TRUE)

    RK31 <- 15/(64*pi)
    psi00 <- psins(r=c(0,0), Sigma=S.star) 
    psihat.star <- vector()
    g.star <- vector()
    
    for (k in 1:nrow(derivt4))
    { 
      r <- derivt4[k,]
      psi1 <- psins(r=r + 2*elem(1, 2), Sigma=S.star)
      psi2 <- psins(r=r + 2*elem(2, 2), Sigma=S.star)

      ## odd order
      if (prod(r) == 3)g.star[k] <- gamse.odd.2d(r=4, n, psi1, psi2, psi00, RK31)
      ## even order
      else  g.star[k] <- gamse.even.2d(r=4, n, psi1, psi2)[k]
      psihat.star[k] <- kfe.scalar(x=x.star, deriv.order=r, g=g.star[k], binned=binned, bin.par=bin.par)
    }

    ## create replicated form of psihat
    psihat.star.vec <- rep(0, nrow(derivt4.vec))
    for (k in 1:nrow(derivt4.vec))
      psihat.star.vec[k] <- psihat.star[which.mat(r=derivt4.vec[k,], mat=derivt4)]

    psihat.star <- list(psir=psihat.star.vec, deriv.ind=derivt4.vec)
  }
  
  return(psihat.star)
}


###############################################################################
# Estimate psi functionals using 2-stage plug-in 
#
# Parameters
# x - pre-transformed data points
# pilot - "amse" - different AMSE pilot
#       - "samse" - SAMSE pilot
# Returns
# estimated psi functionals
###############################################################################

psifun2 <- function(x.star, pilot="samse", binned, bin.par, deriv.order=0, verbose=FALSE)
{ 
  d <- ncol(x.star)
  r <- deriv.order
  S.star <- var(x.star)
  n <- nrow(x.star)

  ## pilots are based on (2r+4)-th order derivatives
  ## compute 1 pilot for SAMSE    
  if (pilot=="samse")
  {
    g6.star <- gsamse(S.star, n=n, modr=6)
    psihat6.star <- kfe(x=x.star, G=g6.star^2*diag(d), deriv.order=6, deriv.vec=TRUE, binned=binned, add.index=FALSE, verbose=verbose)
    g.star <- gsamse(S.star, n=n, modr=4, nstage=2, psihat=psihat6.star)
    psihat.star <- kfe(x=x.star, G=g.star^2*diag(d), deriv.order=4, deriv.vec=TRUE, binned=binned, add.index=TRUE, verbose=verbose)
  }
  ## compute different pilots for AMSE
  else if ((pilot=="amse") & (d==2))
  {
    derivt4 <- dmvnorm.deriv(x=rep(0,d), deriv.order=4, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
    derivt4.vec <- dmvnorm.deriv(x=rep(0,d), deriv.order=4, add.index=TRUE, deriv.vec=TRUE, only.index=TRUE)
    derivt6 <- dmvnorm.deriv(x=rep(0,d), deriv.order=6, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
    
    RK31 <- 15/(64*pi)
    RK51 <- 945/(256*pi)
    RK33 <- 225/(256*pi)
    psi00 <- psins(r=rep(0,d), Sigma=S.star) 

    psihat6.star <- vector()
    g6.star <- vector()
    psihat.star <- vector()
    g.star <- vector()
    
    for (k in 1:nrow(derivt6))
    {
      r <- derivt6[k,]
      psi1 <- psins(r=r + 2*elem(1, 2), Sigma=S.star)
      psi2 <- psins(r=r + 2*elem(2, 2), Sigma=S.star)
      if (prod(r) == 5)
        g6.star[k] <- gamse.odd.2d(r=6, n, psi1, psi2, psi00, RK51)
      else if (prod(r) == 9)
        g6.star[k] <- gamse.odd.2d(r=6, n, psi1, psi2, psi00, RK33) 
      else  
        g6.star[k] <- gamse.even.2d(r=6, n, psi1, psi2)[k]
      
      psihat6.star[k] <- kfe.scalar(x=x.star, deriv.order=r, g=g6.star[k], binned=binned, bin.par=bin.par) 
    }
    
    ## pilots are based on 4th order derivatives using 6th order psi functionals
    ## computed above 'psihat6.star'
    
    for (k in 1:nrow(derivt4))
    {
      r <- derivt4[k,]
      psi1 <- psihat6.star[7 - (r + 2*elem(1,2))[1]]
      psi2 <- psihat6.star[7 - (r + 2*elem(2,2))[1]]
      
      if (prod(r) == 3)
        g.star[k] <- gamse.odd.2d(r=4, n, psi1, psi2, psi00, RK31)
      else
        g.star[k] <- gamse.even.2d(r=4, n, psi1, psi2)[k]

      psihat.star[k] <- kfe.scalar(x=x.star, deriv.order=r, g=g.star[k],  binned=binned, bin.par=bin.par) 
    }

    ## create replicated form of psihat
    psihat.star.vec <- rep(0, nrow(derivt4.vec))
    for (k in 1:nrow(derivt4.vec))
      psihat.star.vec[k] <- psihat.star[which.mat(r=derivt4.vec[k,], mat=derivt4)]

    psihat.star <- list(psir=psihat.star.vec, deriv.ind=derivt4.vec)
  }
  
  return(psihat.star)
}


#############################################################################
## Estimate psi functionals for 6-variate data using 1-stage plug-in 
## with unconstrained pilot
##
## Parameters
## x - data points
## Sd4, Sd6 - symmetrizer matrices of order 4 and 6
##
## Returns
## estimated psi functionals
#############################################################################

psifun1.unconstr <- function(x, binned, bgridsize, deriv.order=0, verbose=FALSE)
{
  n <- nrow(x)
  r <- deriv.order
  S <- var(x)
 
  ## stage 1 of plug-in
  G2r4 <- Gns(r=2*r+4,n=n,Sigma=S) 

  vecPsi2r4 <- kfe(x=x, G=G2r4, deriv.order=2*r+4, binned=binned, bgridsize=bgridsize, deriv.vec=TRUE, add.index=FALSE, verbose=verbose) 
  return (vecPsi2r4)
}


#############################################################################
## Estimate psi functionals for 6-variate data using 2-stage plug-in 
## with unconstrained pilot
##
## Parameters
## x - data points
## Sd4, Sd6 - symmetrizer matrices of order 4 and 6
##
## Returns
## estimated psi functionals
############################################################################

psifun2.unconstr <- function(x, rel.tol=10^-10, binned, bgridsize, deriv.order=0, verbose=FALSE, optim.fun="nlm")
{
  d <- ncol(x)
  n <- nrow(x)
  S <- var(x)
  r <- deriv.order
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))
  
  ## stage 1 of plug-in
  G2r6 <- Gns(r=2*r+6,n=n,Sigma=S) 
  vecPsi2r6 <- kfe(x=x, G=G2r6, binned=binned, bgridsize=bgridsize, deriv.order=2*r+6, deriv.vec=TRUE, add.index=FALSE, verbose=verbose)

  ## asymptotic squared bias for r = 4 for MSE-optimal G
  D2r4phi0 <- DrL0(d=d, r=2*r+4)
  Id2r4 <- diag(d^(2*r+4))
  
  AB2<-function(vechG){
    rr <- 2*r+4
    G <- invvech(vechG)%*%invvech(vechG)
    G12 <- matrix.sqrt(G)
    Ginv12 <- chol2inv(chol(G12))
    AB <- n^(-1)*det(Ginv12)*(Kpow(A=Ginv12,pow=rr)%*%D2r4phi0)+(1/2)*(t(vec(G))%x%Id2r4) %*% vecPsi2r6
    return (sum(AB^2))
  }
      
  Gstart <- Gns(r=2*r+4,n=n,Sigma=S) 
  Gstart <- matrix.sqrt(Gstart)
  
  if (optim.fun1=="nlm")
  {
    res <- nlm(p=vech(Gstart), f=AB2, print.level=2*as.logical(verbose))    
    G2r4 <- res$estimate
  }
  else if (optim.fun1=="optim")
  {
    res <- optim(vech(Gstart), AB2, control=list(reltol=rel.tol, trace=as.numeric(verbose)))
    G2r4 <- res$par
  }
  G2r4 <- invvech(G2r4)%*%invvech(G2r4)

  ## stage 2 of plug-in
  vecPsi2r4 <- kfe(x=x, G=G2r4, binned=binned, bgridsize=bgridsize, deriv.order=2*r+4, deriv.vec=TRUE, add.index=FALSE, verbose=verbose)
  
  return (vecPsi2r4)
}



#############################################################################
# Plug-in bandwidth selectors
#############################################################################

    
############################################################################
## Computes plug-in full bandwidth matrix - 2 to 6 dim
##
## Parameters
## x - data points
## Hstart - initial value for minimisation
## nstage - number of plug-in stages (1 or 2)
## pilot - "amse" - different AMSE pilot
##       - "samse" - SAMSE pilot
##       - "unconstr" - unconstrained pilot
## pre - "scale" - pre-scaled data
##     - "sphere"- pre-sphered data 
##
## Returns
## Plug-in full bandwidth matrix
###############################################################################

hpi <- function(x, nstage=2, binned=TRUE, bgridsize, deriv.order=0)
{
  ## 1-d selector is taken from KernSmooth's dpik
  d <- 1 
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  if (deriv.order==0) h <- dpik(x=x, level=nstage, gridsize=bgridsize)
  else
  {
    n <- length(x)
    d <- 1
    r <- deriv.order  
    K2r4 <- dnorm.deriv(x=0, mu=0, sigma=1, deriv.order=2*r+4)
    K2r6 <- dnorm.deriv(x=0, mu=0, sigma=1, deriv.order=2*r+6) 
    m2 <- 1  
    mr <- psins.1d(r=2*r, sigma=1)

    ## formula for bias annihilating bandwidths from Wand & Jones (1995, p.70)
    if (nstage==2)
    {
      psi2r8.hat <- psins.1d(r=2*r+8, sigma=sd(x))
      gamse2r6 <- (2*K2r6/(-m2*psi2r8.hat*n))^(1/(2*r+9)) 
      psi2r6.hat <- kfe.1d(x=x, g=gamse2r6, deriv.order=2*r+6, inc=1, binned=binned)
      gamse2r4 <- (2*K2r4/(-m2*psi2r6.hat*n))^(1/(2*r+7))
      psi2r4.hat <- kfe.1d(x=x, g=gamse2r4, deriv.order=2*r+4, inc=1, binned=binned)
    }
    else 
    {
      psi2r6.hat <- psins.1d(r=2*r+6, sigma=sd(x))
      gamse2r4 <- (2*K2r4/(-m2*psi2r6.hat*n))^(1/(2*r+7))
      psi2r4.hat <- kfe.1d(x=x, g=gamse2r4, deriv.order=2*r+4, inc=1, binned=binned)
    }
    
    ## formula form Wand & Jones (1995, p.49)
    h <- ((2*r+1)*mr/(m2^2*psi2r4.hat*n))^(1/(2*r+5))
    }
  return(h)
}

Hpi <- function(x, nstage=2, pilot, pre="sphere", Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm")
{
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
 
  if (!is.matrix(x)) x <- as.matrix(x)
  if (missing(pilot)) {if (d==2 & r==0) pilot <- "samse" else pilot <- "dscalar"}
  pilot1 <- match.arg(pilot, c("amse", "samse", "unconstr", "dunconstr", "dscalar"))
  pre1 <- match.arg(pre, c("scale", "sphere"))
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))
 
  if (pilot1=="amse" & (d>2 | r>0)) stop("amse pilot selectors not defined for d>2 and/or r>0")
  if ((pilot1=="samse" | pilot1=="unconstr") & r>0) stop("dscalar or dunconstr pilot selectors are better for derivatives r>0")
  if (pilot1=="unconstr" & d>=6) stop("Unconstrained pilots are not implemented for d>6")
  
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
  
  Idr <- diag(d^r)
  RKr <- nu(r=r, diag(d))*2^(-d-r)*pi^(-d/2)

  if (d > 4) binned <- FALSE
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  
  if (pilot1=="unconstr")
  {
    ## psi4.mat is on data scale
    if (nstage==1)
      psi.fun <- (-1)^r*psifun1.unconstr(x=x, binned=binned, bgridsize=bgridsize, deriv.order=r, verbose=verbose)
    else if (nstage==2)
      psi.fun <- psifun2.unconstr(x=x, binned=binned, bgridsize=bgridsize, deriv.order=r, verbose=verbose)
    psi2r4.mat <- (-1)^r*invvec(psi.fun)   
    
    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) Hstart <- Hns(x=x, deriv.order=r)
  }
  else if (pilot1=="dunconstr")
  {
    ## G2r4 is on data scale
    G2r4 <- Gdunconstr(x=x, d=d, r=r, n=n, nstage=nstage, verbose=verbose, binned=binned, optim.fun=optim.fun)
    vecPsi2r4 <- kfe(x=x, G=G2r4, binned=binned, deriv.order=2*r+4, deriv.vec=TRUE, add.index=FALSE, verbose=verbose)
    if (missing(Hstart)) Hstart <-  Hns(x=x, deriv.order=r)
  }
  else if (pilot1=="dscalar")
  {
    ## g2r4 is on pre-transformed data scale
    g2r4 <- gdscalar(x=x.star, r=r, n=n, d=d, verbose=verbose, nstage=nstage, binned=binned)

    G2r4 <- g2r4^2 * diag(d)
    vecPsi2r4 <- kfe(x=x.star, G=G2r4, binned=binned, deriv.order=2*r+4, deriv.vec=TRUE, add.index=FALSE, verbose=verbose)
    if (missing(Hstart)) Hstart <-  Hns(x=x.star, deriv.order=r)
  }
  else
  {
    if (binned)
    {
      H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RKr)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
      bin.par.star <- binning(x=x.star, bgridsize=bgridsize, H=H.max) 
    }
      
    ## psi4.mat is on pre-transformed data scale
    if (nstage==1)
      psi.fun <- psifun1(x.star, pilot=pilot, binned=binned, bin.par=bin.par.star, deriv.order=r, verbose=verbose)$psir
    else if (nstage==2)
      psi.fun <- psifun2(x.star, pilot=pilot, binned=binned, bin.par=bin.par.star, deriv.order=r, verbose=verbose)$psir
    psi2r4.mat <- invvec(psi.fun)

    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) Hstart <- Hns(x=x.star, deriv.order=r)
    else Hstart <- Sinv12 %*% Hstart %*% Sinv12
  }

  ## PI is estimate of AMISE
  pi.temp <- function(vechH)
  { 
    H <- invvech(vechH) %*% invvech(vechH)
    Hinv <- chol2inv(chol(H))
    IdrvH <- Idr%x%vec(H)
    int.var <- 1/(det(H)^(1/2)*n)*nur(r=r, A=Hinv, mu=rep(0,d), Sigma=diag(d))*2^(-d-r)*pi^(-d/2)
    
    if (pilot1=="dunconstr" | pilot1=="dscalar") 
    {
      pi.val <- int.var + (-1)^r*1/4*vecPsi2r4 %*% (vec(diag(d^r) %x% vec(H) %x% vec(H)))
    }
    else
      pi.val <- int.var + (-1)^r*1/4* sum(diag(t(IdrvH) %*% psi2r4.mat %*% IdrvH))
    pi.val <- drop(pi.val)
    return(pi.val)
  }

  Hstart <- matrix.sqrt(Hstart)
  if (optim.fun1=="nlm")
  {
    result <- nlm(p=vech(Hstart), f=pi.temp, print.level=2*as.numeric(verbose))    
    H <- invvech(result$estimate) %*% invvech(result$estimate)
    amise.star <- result$minimum
  }
  else if (optim.fun1=="optim")
  {
    result <- optim(vech(Hstart), pi.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- invvech(result$par) %*% invvech(result$par)
    amise.star <- result$value
  }
  if (!(pilot1 %in% c("dunconstr","unconstr")))  H <- S12 %*% H %*% S12   ## back-transform
  
  if (!amise) return(H)
  else return(list(H = H, PI.star=amise.star))
}     



###############################################################################
## Computes plug-in diagonal bandwidth matrix for 2 to 6-dim
##
## Parameters
## x - data points
## nstage - number of plug-in stages (1 or 2)
## pre - "scale" - pre-scaled data
##     - "sphere"- pre-sphered data 
##
## Returns
## Plug-in diagonal bandwidth matrix
###############################################################################

Hpi.diag <- function(x, nstage=2, pilot, pre="scale", Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm")
{
  if(!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  RK <- (4*pi)^(-d/2)
 
  if (missing(pilot)) {if (d==2 & r==0) pilot <- "samse" else pilot <- "dscalar"}
  pilot1 <- match.arg(pilot, c("amse", "samse", "unconstr", "dunconstr", "dscalar"))
  pre1 <- match.arg(pre, c("scale", "sphere"))
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))
 
  if (pre1=="sphere") stop("Using pre-sphering won't give diagonal bandwidth matrix")
  if (pilot1=="amse" & (d>2 | r>0)) stop("samse pilot selectors are better for higher dimensions and/or deriv.order>0")
  if (pilot1=="samse" & r>0) stop("dscalar or dunconstr pilot selectors are better for derivatives r>0")
  if (pilot1=="unconstr" | pilot1=="dunconstr") stop("Unconstrained pilot selectors are not suitable for Hpi.diag")
  
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
  
  if (d > 4) binned <- FALSE
  if (binned)
  {
    if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    bin.par <- binning(x=x.star, bgridsize=bgridsize, H=H.max) 
  }
 
  Idr <- diag(d^r)  
  if (pilot1=="amse" | pilot1=="samse")
  {
    if (nstage==1)
      psi.fun <- psifun1(x.star, pilot=pilot, binned=binned, bin.par=bin.par, deriv.order=r, verbose=verbose)$psir
    else if (nstage==2)
      psi.fun <- psifun2(x.star, pilot=pilot, binned=binned, bin.par=bin.par, deriv.order=r, verbose=verbose)$psir
    psi2r4.mat <- invvec(psi.fun)
  }
  else if (pilot1=="dscalar")
  {
    g2r4 <- gdscalar(x=x.star, r=r, n=n, d=d, verbose=verbose, nstage=nstage, binned=binned)
    G2r4 <- g2r4^2 * diag(d)
    vecPsi2r4 <- kfe(x=x.star, G=G2r4, binned=binned, deriv.order=2*r+4, deriv.vec=TRUE, add.index=FALSE, verbose=verbose)
  }
  
  if (d==2 & r==0 & (pilot1=="amse" | pilot1=="samse"))
  {
    ## diagonal bandwidth matrix for 2-dim has exact formula 
    psi40 <- psi.fun[1]
    psi22 <- psi.fun[6]
    psi04 <- psi.fun[16]
    s1 <- sd(x[,1])
    s2 <- sd(x[,2])
    h1 <- (psi04^(3/4)*RK/(psi40^(3/4)*(sqrt(psi40*psi04)+psi22)*n))^(1/6)
    h2 <- (psi40/psi04)^(1/4) * h1
    H <- diag(c(s1^2*h1^2, s2^2*h2^2))
    psimat4.D <- invvech(c(psi40, psi22, psi04))
    amise.star <- drop(n^(-1)*RK*(h1*h2)^(-1) + 1/4*c(h1,h2)^2 %*% psimat4.D %*% c(h1,h2)^2)
  }
  else
  {  
    ## PI is estimate of AMISE
    pi.temp <- function(diagH)
    { 
      H <- diag(diagH) %*% diag(diagH)
      Hinv <- chol2inv(chol(H))
      IdrvH <- Idr%x%vec(H)
      int.var <- 1/(det(H)^(1/2)*n)*nu(r=r, Hinv)*2^(-d-r)*pi^(-d/2)

      if (pilot1=="dscalar")
        pi.val <- int.var + (-1)^r*1/4*vecPsi2r4 %*% (vec(diag(d^r) %x% vec(H) %x% vec(H)))
      else
        pi.val <- int.var + (-1)^r*1/4* sum(diag(t(IdrvH) %*% psi2r4.mat %*% IdrvH))
      
      return(drop(pi.val))
    }
    
    ## use normal reference bandwidth as initial condition
    if (missing(Hstart)) Hstart <- Hns(x=x.star, deriv.order=r)
    else Hstart <- Sinv12 %*% Hstart %*% Sinv12
    Hstart <- matrix.sqrt(Hstart)
    
    if (optim.fun1=="nlm")
    {
      result <- nlm(p=diag(Hstart), f=pi.temp, print.level=2*as.numeric(verbose))    
      H <- diag(result$estimate^2)
      amise.star <- result$minimum
    }
    else if (optim.fun1=="optim")
    {  
      result <- optim(diag(Hstart), pi.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
      H <- diag(result$par) %*% diag(result$par)
      amise.star <- result$value
    }
    
    H <- S12 %*% H %*% S12  ## back-transform
  }

  if (!amise) return(H)
  else return(list(H = H, PI.star=amise.star))
}


###############################################################################
## Cross-validation bandwidth selectors
###############################################################################

###############################################################################
## Computes the least squares cross validation LSCV function for 2 to 6 dim
## 
## Parameters
## x - data values
## H - bandwidth matrix
##
## Returns
## LSCV(H)
###############################################################################

lscv.1d <- function(x, h, binned, bin.par, deriv.order=0)
{
  r <- deriv.order
  lscv1 <- kfe.1d(x=x, g=sqrt(2)*h, inc=1, binned=binned, bin.par=bin.par, deriv.order=2*r) 
  lscv2 <- kfe.1d(x=x, g=h, inc=0, binned=binned, bin.par=bin.par, deriv.order=2*r) 
  return((-1)^r*(lscv1 - 2*lscv2))     
}

lscv.mat <- function(x, H, binned=FALSE, bin.par, bgridsize, deriv.order=0)
{
  r <- deriv.order
  d <- ncol(x)
  n <- nrow(x)

  if (!binned)
  {
    lscv1 <- Qr(x=x, deriv.order=2*r, Sigma=2*H, inc=1)
    lscv2 <- Qr(x=x, deriv.order=2*r, Sigma=H, inc=0)
    lscv <- drop(lscv1 - 2*lscv2)
  }
  else
  {
    lscv1 <- kfe(x=x, G=2*H, inc=1, binned=binned, bin.par=bin.par, bgridsize=bgridsize, deriv.order=2*r, add.index=FALSE)
    lscv2 <- kfe(x=x, G=H, inc=0, binned=binned, bin.par=bin.par, bgridsize=bgridsize, deriv.order=2*r, add.index=FALSE)
     lscv <- (-1)^2*drop((lscv1 - 2*lscv2) %*% vec(diag(d^r)))
  }

  return(lscv)  
}


   
###############################################################################
## Finds the bandwidth matrix that minimises LSCV for 2 to 6 dim
## 
## Parameters
## x - data values
## Hstart - initial bandwidth matrix
##
## Returns
## H_LSCV
###############################################################################

hlscv <- function(x, binned=TRUE, bgridsize, amise=FALSE, deriv.order=0)
{
  if (any(duplicated(x)))
    warning("Data contain duplicated values: LSCV is not well-behaved in this case")
  n <- length(x)
  d <- 1
  r <- deriv.order
  hnorm <- sqrt((4/(n*(d + 2)))^(2/(d + 4)) * var(x))

  if (binned)
  {
    if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
    bin.par <- binning(x, bgridsize=bgridsize, h=hnorm)
    lscv.1d.temp <- function(h) { return(lscv.1d(x=x, h=h, binned=binned, bin.par=bin.par, deriv.order=r)) }
  }
  else
  {
    if (r>0) stop("Unbinned hlscv not yet implemented for deriv.order>0") 
    difs<-x%*%t(rep(1,n))-rep(1,n)%*%t(x)
    difs<-difs[lower.tri(difs)]  
    edifs<-exp(-difs^2/2)
    RK<-1/(2*sqrt(pi))

    lscv.1d.temp <- function(h)
    {
      lscv1 <- (1-1/n)*sum(edifs^(1/(2*h^2)))/(h*sqrt(2)*sqrt(2*pi))
      lscv2 <- 2*sum(edifs^(1/h^2))/(h*sqrt(2*pi))
      return(RK/(n*h)+2*(lscv1-lscv2)/(n^2-n))
    }    
  }
  opt <- optimise(f=lscv.1d.temp, interval=c(0.2*hnorm, 5*hnorm, tol=.Machine$double.eps))

  if (!amise) return(opt$minimum)
  else return(list(h=opt$minimum, LSCV=opt$objective))
}


Hlscv <- function(x, Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm", trunc)
{
  if (any(duplicated(x))) warning("Data contain duplicated values: LSCV is not well-behaved in this case")
  if (!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order 
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))
  
  ## use normal reference selector as initial condn
  Hnorm <- Hns(x=x, deriv.order=r)
  if (missing(Hstart)) Hstart <- Hnorm
  Hstart <- matrix.sqrt(Hstart)
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  if (d>4) binned <- FALSE
  if (binned) bin.par <- binning(x=x, H=Hnorm, bgridsize=bgridsize)
  if (missing(trunc)) {if (deriv.order==0) trunc <- 1e10 else trunc <- 4}

  lscv.init <- lscv.mat(x=x, H=Hnorm, binned=binned, bin.par=bin.par, deriv.order=r)
  lscv.mat.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    lscv <- lscv.mat(x=x, H=H, binned=binned, bin.par=bin.par, deriv.order=r)
    if (det(H) < 1/trunc*det(Hnorm) | det(H) > trunc*det(Hnorm) | abs(lscv) > trunc*abs(lscv.init)) lscv <- lscv.init 
    return(lscv)  
  }
   
  if (optim.fun1=="nlm")
  {
    result <- nlm(p=vech(Hstart), f=lscv.mat.temp, print.level=2*as.numeric(verbose))   
    H <- invvech(result$estimate) %*% invvech(result$estimate)
    amise.opt <- result$minimum
  }
  else if (optim.fun1=="optim")
  {
    result <- optim(vech(Hstart), lscv.mat.temp, control=list(trace=as.numeric(verbose)), method="Nelder-Mead")     
    H <- invvech(result$par) %*% invvech(result$par)
    amise.opt <- result$value
  }

  if (!amise) return(H)
  else return(list(H=H, LSCV=amise.opt))
}

###############################################################################
## Finds the diagonal bandwidth matrix that minimises LSCV for 2 to 6 dim
## 
## Parameters
## x - data values
## Hstart - initial bandwidth matrix
##
## Returns
## H_LSCV,diag
###############################################################################

Hlscv.diag <- function(x, Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm", trunc)
{
  if (any(duplicated(x))) warning("Data contain duplicated values: LSCV is not well-behaved in this case")
  if (!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))
  
  Hnorm <- Hns(x=x, deriv.order=r)
  if (missing(Hstart)) Hstart <- Hnorm
  if (d>4) binned <- FALSE

  ## don't truncate optimisation for deriv.order==0
  if (missing(trunc)) {if (deriv.order==0) trunc <- 1e10 else trunc <- 4}
  
  ## linear binning
  if (binned)
  {
    if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
    bin.par <- binning(x=x, bgridsize=bgridsize, H=Hnorm)
  }
  
  lscv.init <- lscv.mat(x=x, H=Hnorm, binned=binned, bin.par=bin.par, deriv.order=r)

  lscv.mat.temp <- function(diagH)
  {
    H <- diag(diagH^2)
    lscv <- lscv.mat(x=x, H=H, binned=binned, bin.par=bin.par, deriv.order=r)
    if (det(H) < 1/trunc*det(Hnorm) | det(H) > trunc*det(Hnorm) | abs(lscv) > trunc*abs(lscv.init)) lscv <- lscv.init 
    return(lscv)  
  }

  Hstart <- matrix.sqrt(Hstart)
  if (optim.fun1=="nlm")
  {
    result <- nlm(p=diag(Hstart), f=lscv.mat.temp, print.level=2*as.numeric(verbose))
    H <- diag(result$estimate^2)
    amise.opt <- result$minimum
  }
  else if (optim.fun1=="optim")
  {
    result <- optim(diag(Hstart), lscv.mat.temp, method="Nelder-Mead", control=list(trace=as.numeric(verbose)))
    H <- diag(result$par^2)
    amise.opt <- result$value
  }
  
  if (!amise) return(H)
  else return(list(H=H, LSCV=amise.opt))
}

hucv <- function(...) { hlscv(...) }
Hucv <- function(...) { Hlscv(...) }
Hucv.diag <- function(...) { Hlscv.diag(...) }

###############################################################################
## Computes the biased cross validation BCV function for 2-dim
## 
## Parameters
## x - data values
## H1, H2 - bandwidth matrices
##
## Returns
## BCV(H)
###############################################################################

bcv.mat <- function(x, H1, H2, binned=FALSE)
{
  n <- nrow(x)
  d <- 2

  psi <- kfe(x, G=H2, deriv.order=4, add.index=TRUE, deriv.vec=TRUE, inc=0, binned=binned)
  psi40 <- psi$psir[1]
  psi31 <- psi$psir[2]
  psi22 <- psi$psir[4]
  psi13 <- psi$psir[8]
  psi04 <- psi$psir[16]
  
  coeff <- c(1, 2, 1, 2, 4, 2, 1, 2, 1)
  psi.fun <- c(psi40, psi31, psi22, psi31, psi22, psi13, psi22, psi13,psi04)/(n*(n-1))
  psi4.mat <- matrix(coeff * psi.fun, ncol=3, nrow=3)
  
  RK <- (4*pi)^(-d/2) 
  bcv <- drop(n^(-1)*det(H1)^(-1/2)*RK + 1/4*t(vech(H1)) %*% psi4.mat %*% vech(H1))
  
  return(list(bcv=bcv, psimat=psi4.mat))
}



###############################################################################
## Find the bandwidth matrix that minimises the BCV for 2-dim
## 
## Parameters
## x - data values
## whichbcv - 1 = BCV1
##          - 2 = BCV2 
## Hstart - initial bandwidth matrix
##
## Returns
## H_BCV
###############################################################################

Hbcv <- function(x, whichbcv=1, Hstart, binned=FALSE, amise=FALSE, verbose=FALSE)
{
  n <- nrow(x)
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)
  if(!is.matrix(x)) x <- as.matrix(x)
  
  ## use normal reference b/w matrix for bounds
  k <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*n*gamma((d+8)/2)*(d+2)))^(2/(d+4))
  Hmax <- k * abs(var(x))
  up.bound <- Hmax
  if (missing(Hstart)) Hstart <- matrix.sqrt(0.9*Hmax)

  bcv.mat.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    bcv <- bcv.mat(x=x, H1=H, H2=whichbcv*H, binned=binned)$bcv
    return(bcv)  
  }

  result <- optim(vech(Hstart), bcv.mat.temp, method="L-BFGS-B", upper=vech(matrix.sqrt(up.bound)), lower=-vech(matrix.sqrt(up.bound)), control=list(trace=as.numeric(verbose)))

  H <- invvech(result$par) %*% invvech(result$par)
  amise.opt <- result$value
  if (!amise) return(H)
  else return(list(H = H, BCV=amise.opt))
}

###############################################################################
## Find the diagonal bandwidth matrix that minimises the BCV for 2-dim
## 
## Parameters
## x - data values
## whichbcv - 1 = BCV1
##          - 2 = BCV2
## Hstart - initial bandwidth matrix
##
## Returns
## H_BCV, diag
###############################################################################

Hbcv.diag <- function(x, whichbcv=1, Hstart, binned=FALSE, amise=FALSE, verbose=FALSE)
{
  n <- nrow(x)
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)
  if(!is.matrix(x)) x <- as.matrix(x)

  ## use maximally smoothed b/w matrix for bounds
  k <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*n*gamma((d+8)/2)*(d+2)))^(2/(d+4))
  Hmax <- k * abs(var(x))
  up.bound <- diag(Hmax)
  
  if (missing(Hstart)) Hstart <- 0.9*matrix.sqrt(Hmax)

  bcv.mat.temp <- function(diagH)
  {
    H <- diag(diagH) %*% diag(diagH) 
    return(bcv.mat(x, H, whichbcv*H, binned=binned)$bcv)
  } 
 
  result <- optim(diag(Hstart), bcv.mat.temp, method="L-BFGS-B", upper=sqrt(up.bound), control=list(trace=as.numeric(verbose)))
  H <- diag(result$par) %*% diag(result$par)
  amise.opt <- result$value
  if (!amise) return(H)
  else return(list(H = H, BCV=amise.opt))
}

###############################################################################
## Estimate scalar g_AMSE pilot bandwidth for SCV for 2 to 6 dim
##
## Parameters
## Sigma.star - scaled/ sphered variance matrix
## Hamise - (estimate) of H_AMISE 
## n - sample size
##
## Returns
## g_AMSE pilot bandwidth
###############################################################################

Theta6.elem <- function(d)
{
  Theta6.mat <- list()
  Theta6.mat[[d]] <- list()
  for (i in 1:d)
    Theta6.mat[[i]] <- list()
  
  for (i in 1:d)
    for (j in 1:d)
    {  
      temp <- numeric()
      for (k in 1:d)     
        for (ell in 1:d)    
          temp <- rbind(temp, elem(i,d)+2*elem(k,d)+2*elem(ell,d)+elem(j,d))
      
      Theta6.mat[[i]][[j]] <- temp
    }
  
  return(Theta6.mat)
}

gamse.scv <- function(x.star, d, Sigma.star, Hamise, n, binned=FALSE, bin.par, bgridsize, verbose=FALSE, nstage=1, Theta6=FALSE)
{
  if (nstage==0)
  {
    psihat6.star <- psins(r=6, Sigma=Sigma.star, deriv.vec=TRUE) 
  }
  else if (nstage==1)
  {  
    g6.star <- gsamse(Sigma.star, n, 6) 
    G6.star <- g6.star^2 * diag(d)
    if (Theta6) psihat6.star <- kfe(x=x.star, bin.par=bin.par, deriv.order=6, G=G6.star, deriv.vec=FALSE, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
    else psihat6.star <- kfe(x=x.star, bin.par=bin.par, deriv.order=6, G=G6.star, deriv.vec=TRUE, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
  }

  if (Theta6)
  {
    derivt6 <- dmvnorm.deriv(x=rep(0,d), deriv.order=6, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
    Theta6.mat <- matrix(0, ncol=d, nrow=d)
    Theta6.mat.ind <- Theta6.elem(d)
    for (i in 1:d)
      for (j in 1:d)
      {
        temp <- Theta6.mat.ind[[i]][[j]]
        temp.sum <- 0
        for (k in 1:nrow(temp))
          temp.sum <- temp.sum + psihat6.star[which.mat(temp[k,], derivt6)]
        Theta6.mat[i,j] <- temp.sum 
      }
    
    eye3 <- diag(d)
    D4 <- dupl(d)$d
    trHamise <- tr(Hamise) 

    ## required constants - see thesis
    Cmu1 <- 1/2*t(D4) %*% vec(Theta6.mat %*% Hamise)
    Cmu2 <- 1/8*(4*pi)^(-d/2) * (2*t(D4)%*% vec(Hamise) + trHamise * t(D4) %*% vec(eye3))
    num <- 2 * (d+4) * sum(Cmu2*Cmu2)
    den <- -(d+2) * sum(Cmu1*Cmu2) + sqrt((d+2)^2 * sum(Cmu1*Cmu2)^2 + 8*(d+4)*sum(Cmu1*Cmu1) * sum(Cmu2*Cmu2))
    gamse <- (num/(den*n))^(1/(d+6)) 
  }
  else
  {  
    ## updated constants using Chacon & Duong (2010) notation
    Cmu1Cmu1 <- drop(1/4*psihat6.star %*% (Hamise %x% diag(d^4) %x% Hamise) %*% psihat6.star)
    Cmu1Cmu2 <- 3/4*(4*pi)^(-d/2)*drop(vec(Hamise %x% diag(d) %x% Hamise) %*% psihat6.star)
    Cmu2Cmu2 <- 1/64*(4*pi)^(-d)*(4*tr(Hamise%*%Hamise) + (d+8)*tr(Hamise)^2)
    num <- 2 * (d+4) * Cmu2Cmu2
    den <- -(d+2) * Cmu1Cmu2 + sqrt((d+2)^2 * Cmu1Cmu2^2 + 8*(d+4)*Cmu1Cmu1 * Cmu2Cmu2)
    gamse <- (num/(den*n))^(1/(d+6)) 
  }
  return(gamse)
}


###############################################################################
## Estimate unconstrained G_AMSE pilot bandwidth for SCV for 2 to 6 dim
## (J.E. Chacon)
##
## Returns
## G_AMSE pilot bandwidth
###############################################################################

Gunconstr.scv <- function(x, binned=FALSE, bin.par, bgridsize, rel.tol=10^-10, verbose=FALSE, nstage=1, optim.fun="nlm")
{
  d <- ncol(x)
  n <- nrow(x)
  S <- var(x)
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))
  
  ## stage 1 of plug-in
  if (nstage==1)
  {  
    G6 <- (2^(d/2+5)/((d+6)*n))^(2/(d+8))*S
    psihat6 <- kfe(x=x, deriv.order=6, G=G6, deriv.vec=TRUE, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
  }
  else if (nstage==0)
  {
    psihat6 <- psins(r=6, Sigma=S, deriv.vec=TRUE) 
  }
  
  ## constants for normal reference
  D4phi0 <- DrL0(d=d,r=4)
  Id4 <- diag(d^4)

  ## asymptotic squared bias for r = 4
  AB2<-function(vechG){
    G <- invvech(vechG)%*%invvech(vechG)
    G12 <- matrix.sqrt(G)
    Ginv12 <- chol2inv(chol(G12))
    AB <- n^(-1)*det(Ginv12)*(Kpow(A=Ginv12,pow=4)%*%D4phi0)*2^(-(d+4)/2) + (t(vec(G))%x%Id4)%*%psihat6
    return (sum(AB^2))
  }

  Hstart <- matrix.sqrt(Hns(x=x, deriv.order=1))
  if (optim.fun1=="nlm")
  {
     result <- nlm(p=vech(Hstart), f=AB2, print.level=2*as.logical(verbose))    
     G4 <- result$estimate
  }
  else if (optim.fun1=="optim")    
  { 
    result <- optim(vech(Hstart), AB2, method="BFGS", control=list(trace=as.numeric(verbose)))
    G4 <- result$par
  }   
  G4 <- invvech(G4)%*%invvech(G4)
  return(G4) 
}


###############################################################################
## Computes the smoothed cross validation function for 2 to 6 dim
## 
## Parameters
## x - data values
## H - bandwidth matrix
## G - pilot bandwidth matrix
##
## Returns
## SCV(H)
###############################################################################


scv.1d <- function(x, h, g, binned=TRUE, bin.par, inc=1, deriv.order=0)
{
  r <- deriv.order
  if (!missing(x)) n <- length(x)
  if (!missing(bin.par)) n <- sum(bin.par$counts)
  scv1 <- kfe.1d(x=x, deriv.order=2*r, bin.par=bin.par, g=sqrt(2*h^2+2*g^2), binned=binned, inc=inc)
  scv2 <- kfe.1d(x=x, deriv.order=2*r, bin.par=bin.par, g=sqrt(h^2+2*g^2), binned=binned, inc=inc)
  scv3 <- kfe.1d(x=x, deriv.order=2*r, bin.par=bin.par, g=sqrt(2*g^2), binned=binned, inc=inc)

  bias2 <- (-1)^r*(scv1 - 2*scv2 + scv3)
  if (bias2 < 0) bias2 <- 0
  scv <- (n*h)^(-1)*(4*pi)^(-1/2)*2^(-r)*OF(2*r) + bias2

  return(scv)
}

scv.mat <- function(x, H, G, binned=FALSE, bin.par, bgridsize, verbose=FALSE, deriv.order=0)
{
  d <- ncol(x)
  n <- nrow(x)
  r <- deriv.order
  vId <- vec(diag(d))
  Hinv <- chol2inv(chol(H))

  if (!binned)
  {
    scv1 <- Qr(x=x, deriv.order=2*r, Sigma=2*H+2*G, inc=1)
    scv2 <- Qr(x=x, deriv.order=2*r, Sigma=H+2*G, inc=1)
    scv3 <- Qr(x=x, deriv.order=2*r, Sigma=2*G, inc=1)
    bias2 <- (-1)^r*(scv1 - 2*scv2 + scv3)
    if (bias2 < 0) bias2 <- 0
  }
  else
  {
    scv1 <- kfe(x=x, G=2*H + 2*G, deriv.order=2*r, inc=1, binned=binned, bin.par=bin.par, bgridsize=bgridsize, verbose=verbose, add.index=FALSE)
    scv2 <- kfe(x=x, G=H + 2*G, deriv.order=2*r, inc=1, binned=binned, bin.par=bin.par, bgridsize=bgridsize, verbose=verbose, add.index=FALSE)
    scv3 <- kfe(x=x, G=2*G, deriv.order=2*r, inc=1, binned=binned, bin.par=bin.par, bgridsize=bgridsize, verbose=verbose, add.index=FALSE)
    
    bias2 <- drop((-1)^r*Kpow(vId,r) %*% (scv1 - 2*scv2 + scv3))
    if (bias2 < 0) bias2 <- 0
  }
  scvmat <- 1/(det(H)^(1/2)*n)*nur(r=r, A=Hinv, mu=rep(0,d), Sigma=diag(d), type="direct")*2^(-d-r)*pi^(-d/2) + bias2
  return (scvmat)
}


###############################################################################
# Find the bandwidth that minimises the SCV for 1 to 6 dim
# 
# Parameters
# x - data values
# pre - "scale" - pre-scaled data
#     - "sphere"- pre-sphered data
# Hstart - initial bandwidth matrix
#
# Returns
# H_SCV
###############################################################################

hscv <- function(x, nstage=2, binned=TRUE, bgridsize, plot=FALSE)
{
  sigma <- sd(x)
  n <- length(x)
  d <- 1
  hnorm <- sqrt((4/(n*(d + 2)))^(2/(d + 4)) * var(x))
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  hmin <- 0.1*hnorm
  hmax <- 2*hnorm

  bin.par <- binning(x=x, bgridsize=bgridsize, h=hnorm)
  if (nstage==1)
  {
    psihat6 <- psins.1d(r=6, sigma=sigma)
    psihat10 <- psins.1d(r=10, sigma=sigma)
  }
  else if (nstage==2)
  {
    g1 <- (2/(7*n))^(1/9)*2^(1/2)*sigma
    g2 <- (2/(11*n))^(1/13)*2^(1/2)*sigma

    psihat6 <- kfe.1d(x=x, bin.par=bin.par, binned=binned, deriv.order=6, g=g1, inc=1)
    psihat10 <- kfe.1d(x=x, bin.par=bin.par, binned=binned, deriv.order=10, g=g2, inc=1)
  }

  g3 <- (-6/((2*pi)^(1/2)*psihat6*n))^(1/7) 
  g4 <- (-210/((2*pi)^(1/2)*psihat10*n))^(1/11)
  psihat4 <- kfe.1d(x=x, bin.par=bin.par, binned=binned, deriv.order=4, g=g3, inc=1)
  psihat8 <- kfe.1d(x=x, bin.par=bin.par, binned=binned, deriv.order=8, g=g4, inc=1)

  C <- (441/(64*pi))^(1/18) * (4*pi)^(-1/5) * psihat4^(-2/5) * psihat8^(-1/9)
  
  scv.1d.temp <- function(h)
  {
    return(scv.1d(x=x, bin.par=bin.par, h=h, g=C*n^(-23/45)*h^(-2), binned=binned, inc=1))
  }

  if (plot)
  {  
    hseq <- seq(hmin,hmax, length=400)
    hscv.seq <- rep(0, length=length(hseq))
    for (i in 1:length(hseq))
      hscv.seq[i] <- scv.1d.temp(hseq[i])
    plot(hseq, hscv.seq, type="l", xlab="h", ylab="SCV(h)")
  }
  
  opt <- optimise(f=scv.1d.temp, interval=c(hmin, hmax))$minimum
  if (n >= 1e5) warning("hscv is not always stable for large samples")
  
  return(opt)
}


Hscv <- function(x, nstage=2, pre="sphere", pilot, Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm")
{
  n <- nrow(x)
  d <- ncol(x)
  ##S <- var(x)
  r <- deriv.order

  if(!is.matrix(x)) x <- as.matrix(x)
  if (missing(pilot)) {if (d==2 & r==0) pilot <- "samse" else pilot <- "dscalar"}
  pilot1 <- match.arg(pilot, c("amse", "samse", "unconstr", "dunconstr", "dscalar"))
  pre1 <- match.arg(pre, c("scale", "sphere"))
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))
  
  if (pilot1=="amse" & (d>2 | r>0)) stop("amse pilot selectors not defined for d>2 and/or r>0")
  if ((pilot1=="samse" | pilot1=="unconstr") & r>0) stop("dscalar or dunconstr pilot selectors are better for deriv.order>0")
    
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
    
  RK <- (4*pi)^(-d/2)

  if (d > 4) binned <- FALSE
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  if (d>=4 & nstage==2) bgridsize <- rep(11,d)
  
  if (binned)
  {
    if (pilot1=="unconstr" | pilot1=="dunconstr")
        H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x)
    else
        H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    if (default.bflag(d=d, n=n))
        bin.par <- binning(x=x.star, bgridsize=bgridsize, H=matrix.sqrt(H.max))
  }
  if (pilot1=="unconstr")
  {
    ## Gu pilot matrix is on data scale
    Gu <- Gunconstr.scv(x=x, binned=binned, bgridsize=bgridsize, verbose=verbose, nstage=nstage-1, optim.fun=optim.fun)
    if (missing(Hstart)) Hstart <- Hns(x=x, deriv.order=r)
  }
  else if (pilot1=="dunconstr")
  {
    ## Gu pilot matrix is on data scale
    Gu <- Gdunconstr(x=x, d=d, r=r, n=n, nstage=nstage, verbose=verbose, binned=binned, scv=TRUE,  optim.fun=optim.fun)
  
    if (missing(Hstart)) Hstart <- Hns(x=x, deriv.order=r)
  }
  else if (pilot1=="dscalar")
  {
    ## Gs is on pre-transformed data scale
    g2r4 <- gdscalar(x=x.star, d=d, r=r, n=n, nstage=nstage, verbose=verbose, scv=TRUE, binned=binned)
    Gs <- g2r4^2*diag(d)
    if (missing(Hstart)) Hstart <-Hns(x=x.star, deriv.order=r)
  }
  else
  {
    ## Gs is on transformed data scale    
    Hamise <- Hpi(x=x.star, nstage=1, deriv.order=r, pilot=pilot, pre="sphere", binned=TRUE, bgridsize=bgridsize, verbose=verbose, optim.fun=optim.fun) 
    if (any(is.na(Hamise)))
    {
      warning("Pilot bandwidth matrix is NA - replaced with maximally smoothed")
      Hamise <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    }
    gs <- gamse.scv(x.star=x.star, d=d, Sigma.star=var(x.star), Hamise=Hamise, n=n, binned=binned, bgridsize=bgridsize, verbose=verbose, nstage=nstage-1)
    Gs <- gs^2*diag(d)

    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) Hstart <- Hns(x=x.star, deriv.order=r) 
    else Hstart <- Sinv12 %*% Hstart %*% Sinv12
  }

  ## SCV is estimate of AMISE
  
  scv.mat.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    if (pilot1=="samse" | pilot1=="amse" | pilot1=="dscalar"){ Gpilot <- Gs; xx <- x.star }
    else if (pilot1=="unconstr" | pilot1=="dunconstr") { Gpilot <- Gu; xx <- x }

    if (default.bflag(d=d, n=n))
        scvm <- scv.mat(x=xx, H=H, G=Gpilot, binned=binned, bin.par=bin.par, verbose=FALSE, deriv.order=r)
    else
        scvm <- scv.mat(x=xx, H=H, G=Gpilot, binned=binned, verbose=FALSE, deriv.order=r)
    return(scvm)
  }

  Hstart <- matrix.sqrt(Hstart)
  if (optim.fun1=="nlm")
  {
    result <- nlm(p=vech(Hstart), f=scv.mat.temp, print.level=2*as.numeric(verbose))    
    H <- invvech(result$estimate) %*% invvech(result$estimate)
    amise.star <- result$minimum
  }
  else if (optim.fun=="optim")
  {
    result <- optim(vech(Hstart), scv.mat.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- invvech(result$par) %*% invvech(result$par)
    amise.star <- result$value
  }
  if (!(pilot1 %in% c("dunconstr","unconstr")))  H <- S12 %*% H %*% S12   ## back-transform

  if (!amise) return(H)
  else return(list(H = H, SCV.star=amise.star))
}


Hscv.diag <- function(x, nstage=2, pre="scale", pilot, Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm")
{
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  RK <- (4*pi)^(-d/2)
  
  if(!is.matrix(x)) x <- as.matrix(x)
  if (missing(pilot)) {if (d==2 & r==0) pilot <- "samse" else pilot <- "dscalar"}
  pilot1 <- match.arg(pilot, c("amse", "samse", "unconstr", "dunconstr", "dscalar"))
  pre1 <- match.arg(pre, c("scale", "sphere"))
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))
  
  if (pilot1=="amse" & (d>2 | r>0)) stop("samse pilot selectors are better for higher dimensions and/or deriv.order>0")
  if (pilot1=="samse" & r>0) stop("dscalar pilot selectors are better for deriv.order>0")
  if (pilot1=="unconstr" | pilot1=="dunconstr") stop("Unconstrained pilot selectors are not suitable for Hscv.diag")

  if (pre1=="sphere") stop("Using pre-sphering doesn't give a diagonal bandwidth matrix")

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

  if (d > 4) binned <- FALSE
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  if (binned)
  {
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    bin.par.star <- binning(x=x.star, bgridsize=bgridsize, H=H.max) 
  }
  
  if (pilot1=="dscalar")
  {
    ## Gs is on pre-transformed data scale
    g2r4 <- gdscalar(x=x.star, r=r, n=n, d=d, verbose=verbose, nstage=nstage, scv=TRUE, binned=binned)
    Gs <- g2r4^2*diag(d)
    if (missing(Hstart)) Hstart <- Hns(x=x.star, deriv.order=r)
  }
  else
  {
    ## Gs is on transformed data scale
      Hamise <- Hpi(x=x.star, nstage=1, pilot=pilot, pre="sphere", binned=binned, bgridsize=bgridsize, verbose=verbose, optim.fun=optim.fun) 
      if (any(is.na(Hamise)))
      {
        warning("Pilot bandwidth matrix is NA - replaced with maximally smoothed")
        Hamise <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
      }
      gs <- gamse.scv(x.star=x.star, d=d, Sigma.star=var(x.star), Hamise=Hamise, n=n, binned=binned, bgridsize=bgridsize, verbose=verbose, nstage=nstage-1)
    
    Gs <- gs^2*diag(d)

    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) Hstart <- Hns(x=x.star, deriv.order=r)
    else Hstart <- Sinv12 %*% Hstart %*% Sinv12
  }

  scv.mat.temp <- function(diagH)
  {
    H <- diag(diagH) %*% diag(diagH)
    if (default.bflag(d=d, n=n)) scvm <- scv.mat(x.star, H, Gs, binned=binned, verbose=FALSE, bin.par=bin.par.star, deriv.order=r)
    else scvm <- scv.mat(x.star, H, Gs, binned=binned, verbose=FALSE, deriv.order=r)
    return(scvm)
  }

  Hstart <- matrix.sqrt(Hstart)
  if (optim.fun1=="nlm")
  {
    result <- nlm(p=diag(Hstart), f=scv.mat.temp, print.level=2*as.numeric(verbose))    
    H <- diag(result$estimate) %*% diag(result$estimate)
    amise.star <- result$minimum
  }
  else if (optim.fun1=="optim")
  {  
    result <- optim(diag(Hstart), scv.mat.temp, method="Nelder-Mead", control=list(trace=as.numeric(verbose)))
    H <- diag(result$par) %*% diag(result$par)
    amise.star <- result$value
  }
  ## back-transform
  H <- S12 %*% H %*% S12

  if (!amise) return(H)
  else return(list(H = H, SCV.star=amise.star))
}

##############################################################################
## Normal scale selector H_ns for kernel density derivate estimators
##############################################################################

Hns <- function(x, deriv.order=0)
{
  if (is.vector(x)){ n<-1; d <- length(x)} 
  else { n <- nrow(x); d <- ncol(x)}
  r <- deriv.order
  H <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x)
  return(H)
}

hns <- function(x, deriv.order=0)
{
  n <- length(x)
  d <- 1
  r <- deriv.order
  h <- (4/(n*(d+2*r+2)))^(1/(d+2*r+4))*sd(x)
  return(h)
}

#######################################################################
## Normal scale G_ns for kernel functional estimators
#######################################################################

Gns <- function(r,n,Sigma)
{
  d <- ncol(Sigma)
  G <- (2/((n*(d+r))))^(2/(d+r+2))*2*Sigma
  return(G)
}
