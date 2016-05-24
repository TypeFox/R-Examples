## postc0.R --- 
## Author          : Claus Dethlefsen
## Created On      : Tue Mar 12 06:52:02 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Thu Jul 24 15:12:23 2003
## Update Count    : 100
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


post0 <- function(mu,tau,rho,phi,y,timetrace=FALSE) {
    ## Posterior for continuous node with 0 parents
    if (timetrace) {t1 <- proc.time();cat("[post0 ")}
    
    mu.n  <- (tau*mu+sum(y))/(tau+length(y))
    tau.n <- tau + length(y)
    rho.n <- rho + length(y)    
    phi.n <- phi + (y - mu.n)%*%y + (mu - mu.n)*tau*mu

    s <- as.numeric(phi)/rho*(diag(length(y)) + matrix(1/tau,length(y),length(y)))
    k <- lgamma( (rho + length(y))/2 ) - lgamma(rho/2)-0.5*log(det(rho*s*pi))
    ind <- log( 1 + (mahalanobis(y,center=mu,cov=s,inverted=FALSE))/rho)
    loglik <- k - (rho+length(y))/2 * ind
    
    if (timetrace) {
        t2 <- proc.time()
        cat((t2-t1)[1],"]")
    }
    
    list(mu=mu.n,tau=tau.n,rho=rho.n,phi=phi.n,loglik=loglik)
}


postc0c <- function(mu,tau,rho,phi,y,timetrace=FALSE) {
    ## Posterior for continuous node with 0 parents
    if (timetrace) {t1 <- proc.time();cat("[postc0 ")}
    
    
    ## call to C
    res <- .C("postc0",
              mu =as.double(mu),
              tau=as.double(tau),
              rho=as.double(rho),
              phi=as.double(phi),
              loglik=as.double(0),
              as.double(y),
              as.integer(length(y)),
              PACKAGE="deal"
              )
    
    
    if (timetrace) {
        t2 <- proc.time()
        cat((t2-t1)[1],"]")
    }
    
    list(mu=res$mu,tau=res$tau,rho=res$rho,phi=res$phi,loglik=res$loglik)
}
    
