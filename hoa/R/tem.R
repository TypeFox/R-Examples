## file tem/R/tem.R, v 0.1-2 2009/10/01
##
##  Copyright (C) 2007-2015 Anthony C. Davison 
##
##  This file is part of the "hoa" package for R.  This program 
##  is free software; you can redistribute it and/or modify it under the 
##  terms of the GNU General Public License as published by the Free 
##  Software Foundation; either version 2 of the License, or (at your 
##  option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
##  MA 02111-1307 USA or look up the web page 
##  http://www.gnu.org/copyleft/gpl.html.
##
##  Please send any comments, suggestions or errors found to:
##  Anthony C. Davison, IMA-FSB-EPFL, Station 8, CH-1015 Lausanne, 
##  Switzerland.  Email: Anthony.Davison@epfl.ch.  
##  Web: http://statwww.epfl.ch/davison/davison.html
##  OR
##  Alex-Antoine Fortin, alex@fortin.bio.

tem <- function(psi=NULL, nlogL, phi, make.V, th.init, data, tol=10^(-5), n.psi=50)
{
#psi should be NULL or a scalar
  if (length(psi) == 2 | (is.null(psi) & n.psi==2)) stop("psi must be NULL or not a vector of length 2")
#
# define some ancillary functions
#
even <- function(n) (2*floor(n/2)==n)
#
d.phi <- function(x, V, data, tol=tol)
  {  # numerical differentiation of function phi at x
    f0 <- phi(x, V, data)
    p <- length(f0)
    out <- matrix(NA,p,p)
    for (i in 1:p) 
  {  x1 <- x
     x1[i] <- x1[i]+tol
     out[i,] <- (phi(x1, V, data)-f0)/tol }
# dxd matrix with columns for elements of phi and rows for elements of theta
# first row corresponds to psi, rest to lambda
    out
  }
#
  nlogL.full <- function(th, data)  nlogL(th[1], th[-1], data)
  nlogL.rest <- function(psi, lam, data)  nlogL(psi, lam, data)
#
# full maximum likelihood fit 
#
  suppressWarnings(full <- nlm(nlogL.full, jitter(th.init), hessian=TRUE, data=data))
  if(full$code>2) stop("Full fit failed: try different initial values")
  th.full <- full$estimate
  th.se <- sqrt(diag(solve(full$hessian)))
  psi.se <- th.se[1]
  L.full <- -full$minimum
  out <- NULL
  out$normal <- c(th.full[1],psi.se)
  out$th.hat <- th.full
  out$th.hat.se <- th.se
# 
# compute matrix of Vs for later computation of local parametrization
#
  V <- make.V(th.full, data)
# 
# choose a range of psi or use those given
#
  if(is.null(psi)) psi <- seq(from=th.full[1]-5*psi.se,to=th.full[1]+5*psi.se,length=n.psi)
#
# restricted maximum likelihood fits
#
  n.psi <- length(psi)
  K <- ceiling(n.psi/2)
  if(even(n.psi)) K <- K+1
    L.rest <- J.rest <- psi
  out$th.rest <- matrix(th.full,n.psi,length(th.full),byrow=TRUE)
  if (n.psi == 1)
  {
    suppressWarnings(rest <- nlm(nlogL.rest, jitter(out$th.rest[,-1]), hessian=TRUE, psi=psi, data=data))
    if(rest$code>2) stop("Restricted fit failed: try different initial values")
    out$th.rest <- c(psi,rest$estimate)
    dim(out$th.rest) <- dim(matrix(0, nrow=n.psi, ncol=length(out$th.rest)))
    L.rest <- -rest$minimum
    J.rest <- det(rest$hessian)
  }
  else 
{ 
  for(j in (K-1):1)
  {
    suppressWarnings(rest <- nlm(nlogL.rest, jitter(out$th.rest[j+1,-1]), hessian=TRUE, psi=psi[j], data=data))
    out$th.rest[j,] <- c(psi[j],rest$estimate)
    L.rest[j] <- -rest$minimum
    J.rest[j] <- det(rest$hessian)
  }
  if(even(n.psi))
  { 
    suppressWarnings(rest <- nlm(nlogL.rest, jitter(th.full[-1]), hessian=TRUE, psi=psi[K], data=data))
    out$th.rest[K,] <- c(psi[K],rest$estimate)
    L.rest[K] <- -rest$minimum
    J.rest[K] <- det(rest$hessian)
  }
  else 
  {
    suppressWarnings(rest <- nlm(nlogL.rest, jitter(out$th.rest[K-1,-1]), hessian=TRUE, psi=psi[K], data=data))
    out$th.rest[K,] <- c(psi[K],rest$estimate)
    L.rest[K] <- -rest$minimum
    J.rest[K] <- det(rest$hessian)
  }
  for(j in (K+1):n.psi)
  {
    suppressWarnings(rest <- nlm(nlogL.rest, jitter(out$th.rest[j-1,-1]), hessian=TRUE, psi=psi[j], data=data))
    out$th.rest[j,] <- c(psi[j],rest$estimate)
    L.rest[j] <- -rest$minimum
    J.rest[j] <- det(rest$hessian)
  }
}
#
# compute r
#
  out$r <- sign(out$normal[1]-out$th.rest[,1])*sqrt(2*(L.full-L.rest))
#
# compute q, first doing computations at full MLE
#
  dphi.dth.full <- d.phi(th.full, V, data, tol)
  D.bot <- det(dphi.dth.full)
  j.th.th <- det(full$hessian)
#
# prepare output
#
  out$q <- out$psi <- psi
#
  for (j in 1:n.psi)
{
  dphi.dth.rest <- d.phi(out$th.rest[j,], V, data, tol)
  dphi.dth.rest[1,] <- phi(out$th.hat, V, data) - phi(out$th.rest[j,],V,data)
  D.top <- det(dphi.dth.rest)
  j.lam.lam <- J.rest[j]
  out$q[j] <- (D.top/D.bot)*sqrt(j.th.th/j.lam.lam)
}
#
# output
#
  out$rstar <- out$r + log(out$q/out$r)/out$r
  class(out) <- "fr"
  out
}

fraser.reid <- tem

## Summary S3Methods for fr objects
lik.ci <- function(object, conf=c(0.975,0.025), ...)
{  # read likelihood confidence interval off from fraser-reid object
  fr <- object
  fit.r <- smooth.spline(fr$r, fr$psi)
  fit.rstar <- smooth.spline(fr$rstar, fr$psi)
  p.r <- predict(fit.r, qnorm(conf))
  p.rstar <- predict(fit.rstar, qnorm(conf))
  mle.hoa <- predict(fit.rstar, 0)
  se.hoa <- predict(fit.rstar, 0, deriv=1)
  mle.hoa <- c(mle.hoa$y, -se.hoa$y)
  z.lims <- fr$normal[1]-qnorm(conf)*fr$normal[2]
  r.lims <- p.r$y
  rstar.lims <- p.rstar$y
  
  pointEst.z <- fr$normal[1]
  pointEst.r <- predict(fit.r, 0.5)$y
  pointEst.rstar <- predict(fit.rstar, 0.5)$y
  
  cat("Statistics for the parameter of interest psi\n")
  cat("Point estimate using: \n")
  cat("Wald statistic, z           :",pointEst.z,"\n")
  cat("Likelihood root, r          :",pointEst.r,"\n")
  cat("Modified likelihood root, r*:",pointEst.rstar,"\n\n")
  cat("Confidence intervals, levels:", 1-conf,"\n")
  cat("Wald statistic, z           :", z.lims,"\n")
  cat("Likelihood root, r          :", r.lims,"\n")
  cat("Modified likelihood root, r*:", rstar.lims,"\n")
  invisible(list(mle=fr$normal, mle.hoa=mle.hoa, pointEst.z=pointEst.z,
                 pointEst.r=pointEst.r, pointEst.rstar=pointEst.rstar,
                 z.lims=z.lims, r.lims=r.lims, rstar.lims=rstar.lims))
}

summary.fr <- lik.ci

## Plot S3Methods for fr objects
plot.fr <- function(x, psi=NULL, all=FALSE, xl="Interest parameter", ...)
{ # plot a fraser-reid object
  
  old.pars <- par(no.readonly = TRUE)
  if(all) par(mfrow=c(2,2), pty="s") else par(mfrow=c(1,2), pty="s") 
  
  # top left: plot of pivot as a function of psi
  
  fr <- x
  plot(fr$psi,fr$r,type="l",xlab=xl,ylab="Value of pivot",ylim=c(-4, 4), 
       panel.first=abline(h=qnorm(c(0.005,0.025,0.05,0.5,0.95,0.975,0.995)),col="grey",lwd=0.7))
  lines(fr$psi,(fr$normal[1]-fr$psi)/fr$normal[2],col="green")
  lines(fr$psi,fr$q,col="red")
  lines(fr$psi,fr$r)
  lines(fr$psi,fr$rstar,col="blue")
  #  abline(h=qnorm(c(0.005,0.025,0.05,0,0.95,0.975,0.995)),lty=2)
  legend(sum(fr$normal),4,c("Wald pivot","Likelihood root","Modified root","q(psi)"),
         lty=c(1,1,1,1),col=c("green","black","blue","red"),bty="n")
  
  # top right: log likelihood (and adjusted version, I think?) as a function of psi
  
  plot(fr$psi,-fr$r^2/2,type="l",xlab=xl,ylab="Log likelihood",ylim=c(-8, 0),
       panel.first=abline(h=-qchisq(c(0.95,0.99),df=1)/2,col="grey"),lwd=0.7)
  lines(fr$psi,-fr$rstar^2/2,col="blue")
  legend(fr$normal[1],-7,c("Profile log likelihood","Modified profile log likelihood"),
         lty=c(1,1),col=c("black","blue"),bty="n")  
  corr <- log(fr$q/fr$r)/fr$r
  
  # optional: add diagnostic panels
  
  if (all) {  
    
    # lower left: plot of Phi(pivot) as a function of psi
    
    plot(fr$psi,pnorm(fr$r),type="l",xlab=xl,ylab="Significance function",ylim=c(0,1),
         panel.first=abline(h=c(0.025,0.05,0.5,0.95,0.975),col="grey",lwd=0.7))
    lines(fr$psi,pnorm(fr$q),col="red")
    lines(fr$psi,pnorm(fr$rstar),col="blue")
    legend(sum(fr$normal),0.9,c("Likelihood root","Modified root","q(psi)"),
           lty=c(1,1,1),col=c("black","blue","red"),bty="n")
    
    # lower right: log(q/r)/r as a function of r (should be smooth)
    
    plot(fr$r,corr,type="l",xlab="Likelihood root r",ylab="Correction log(q/r)/r",
         panel.first={ abline(h=0,col="grey"); abline(v=0,col="grey")})
    points(fr$r,log(fr$q/fr$r)/fr$r)
    #  lines(fr$r,fitted(lm(corr~poly(fr$r,7))),col="green")
    #  lines(smooth.spline(fr$r,corr))
    
  }
  
  par(old.pars)
  
  # computation of significance probability if psi is specified
  
  if (!is.null(psi)) 
  {
    fit.r <- smooth.spline(fr$psi, fr$r)
    fit.rstar <- smooth.spline(fr$psi, fr$rstar)
    p.r <- predict(fit.r, psi)
    p.rstar <- predict(fit.rstar, psi)
    z <- c((fr$normal[1]-psi)/fr$normal[2],p.r$y,p.rstar$y)
    cat("z, r, r*   :",z,"\n")
    cat("Sig. levels:",pnorm(z),"\n") 
  }
}

