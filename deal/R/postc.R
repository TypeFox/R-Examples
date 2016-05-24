## postc.R --- 
## Author          : Claus Dethlefsen
## Created On      : Tue Mar 12 06:52:02 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Fri Apr 20 09:25:29 2007
## Update Count    : 144
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Susanne Gammelgaard Bøttcher, Claus Dethlefsen
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################

postc <- function(mu,tau,rho,phi,y,z,timetrace=FALSE) {
    ## Posterior for continuous node with continuous parents
    ## written as a for-loop in R (slow)
    if (timetrace) {t1 <- proc.time();cat("[postc ")}
    
    
    loglik <- 0
    for (i in 1:length(y)) {

        ## likelihood
        logscale <- log(phi) + log( 1 + t(z[i,])%*%solve(tau)%*%z[i,])
        logk     <- lgamma( (rho+1)/2 ) - lgamma(rho/2) - 0.5*(logscale  +  log(pi)) 
        mscore   <- logk - 0.5*(rho+1)*log(1 + ((y[i] - z[i,]%*%mu)^2)/exp(logscale))

        loglik <- loglik + mscore

## update
        oldtau <- tau
        oldmu <- mu
        tau <- tau + z[i,]%*%t(z[i,])
        mu <- solve(tau)%*%(oldtau%*%mu+z[i,]*y[i])
        rho<- rho + 1
        phi<- phi + (y[i]-t(z[i,])%*%mu)*y[i] + t(oldmu-mu)%*%oldtau%*%oldmu
    }

  if (timetrace) {
    t2 <- proc.time()
    cat((t2-t1)[1],"]")
  }

    list(mu=mu,tau=tau,rho=rho,phi=phi,loglik=loglik)
}


post <- function(mu,tau,rho,phi,y,z,timetrace=FALSE) {
    ## Posterior for continuous node with continuous parents
    ## written as matrix notation in R
    if (timetrace) {t1 <- proc.time();cat("[post ")}
    
    mu.n  <- solve(tau+t(z)%*%z)%*%(tau%*%mu+t(z)%*%y)
    tau.n <- tau + t(z)%*%z
    rho.n <- rho + length(y)
    phi.n <- phi + t(y - z%*%mu.n)%*%y + t(mu - mu.n)%*%tau%*%mu

    loglik <- 0
    s <- as.numeric(phi)/rho*(diag(nrow(z))+ z%*%solve(tau)%*%t(z))
    k <- lgamma( (rho + length(y))/2 ) - lgamma(rho/2)-0.5*log(det(rho*s*pi))
    ind <- log( 1 + (mahalanobis(y,center=z%*%mu,cov=s,inverted=FALSE))/rho)
    loglik <- as.numeric(k) - (rho+length(y))/2 * ind

        
    if (timetrace) {
        t2 <- proc.time()
        cat((t2-t1)[1],"]")
    }

    list(mu=mu.n,tau=tau.n,rho=rho.n,phi=phi.n,loglik=loglik)
}

## postM <- function(mu,tau,rho,phi,y,z,timetrace=FALSE) {
##     ## Posterior for continuous node with continuous parents
##     ## written as Matrix notation in R (needs Matrix)
##     if (timetrace) {t1 <- proc.time();cat("[postM ")}
    
##     z <- as.Matrix(z)
##     mu.n  <- solve(as.Matrix(tau+t(z)%*%z))%*%(tau%*%mu+t(z)%*%y)
##     tau.n <- tau + t(z)%*%z
##     rho.n <- rho + length(y)
##     phi.n <- phi + t(y - z%*%mu.n)%*%y + t(mu - mu.n)%*%tau%*%mu


##     loglik <- 0
##     s <- as.numeric(phi)/rho*(diag(nrow(z))+ z%*%solve(tau)%*%t(z))
##     k <- lgamma( (rho + length(y))/2 ) - lgamma(rho/2)-0.5*log(det(rho*s*pi))
##     ind <- log( 1 + (mahalanobis(y,center=z%*%mu,cov=s,inverted=FALSE))/rho)
##     loglik <- as.numeric(k) - (rho+length(y))/2 * ind

        
##     if (timetrace) {
##         t2 <- proc.time()
##         cat((t2-t1)[1],"]")
##     }

##     list(mu=mu.n,tau=tau.n,rho=rho.n,phi=phi.n,loglik=loglik)
## }


postcc <- function(mu,tau,rho,phi,y,z,timetrace=FALSE) {
    ## Posterior for continuous node with x parents
    ## written as for-loop in C (fast)
    if (timetrace) {t1 <- proc.time();cat("[postcc ")}
    
    
    ## call to C
    res <- .C("postc",
              mu =as.double(c(mu)),
              tau=as.double(t(tau)),
              rho=as.double(rho),
              phi=as.double(phi),
              loglik=as.double(0),
              as.double(y),
              as.double(t(z)),
              as.integer(length(y)),
              as.integer(ncol(z)),
              PACKAGE="deal"
              )
    if (timetrace) {
        t2 <- proc.time()
        cat((t2-t1)[1],"]")
    }
    list(mu=res$mu,tau=matrix(res$tau,ncol(z),ncol(z)),rho=res$rho,phi=res$phi,loglik=res$loglik)
}
    
