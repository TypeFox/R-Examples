## conditional.R
## Author          : Claus Dethlefsen
## Created On      : Sun Dec 02 14:18:04 2001
## Last Modified By: Claus Dethlefsen
## Last Modified On: Tue Jul 22 15:31:42 2003
## Update Count    : 291
## Status          : Unknown, Use with caution!
#######################################################################
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

conditional.cont <- function(A,mu,nu,rho,phi) {
    ## Conditional distribution for continuous node with index A
    ## The master parameters mu, nu, rho and phi
    ## See Bottcher (2002) for details.
    
    B <- A ## renaming due to compatibility
    
    ## calculate conditional probabilities
    ## p. 14 in Bottcher
    ##
    A <- setdiff(1:ncol(phi),B)
    if (length(A)<1) A <- TRUE
    
    rho.BlA <- rho + length(A)
    phi.AA.inv <- solve(phi[A,A])
    phi.tmp <- phi[B,A]%*%phi.AA.inv
    phi.BlA <- phi[B,B] - phi.tmp%*%phi[A,B]
    mu.BlA  <- c(mu[B] - phi.tmp%*%mu[A], phi.tmp)
    tau.BlA.inv.11 <- 1/nu + t(mu[A])%*%phi.AA.inv%*%mu[A]
    tau.BlA.inv.22 <- phi.AA.inv
    tau.BlA.inv.12 <- -t(mu[A]%*%phi.AA.inv)
    
    tau.inv <- rbind(cbind(tau.BlA.inv.11,t(tau.BlA.inv.12)),
                     cbind(tau.BlA.inv.12,tau.BlA.inv.22)
                     )

    tau <- solve(tau.inv)

    list(tau=tau,phi=phi.BlA,mu=mu.BlA,rho=rho.BlA)
}

conditional.disc <- function(A,master) {
    list(list(alpha=apply(master,A,sum)))
}

conditional <- function(A,master,nw) {
    ## From node index A and given the master prior, calculate the
    ## conditional of A given the parents. (In nw, we use parents,
    ## discrete and continuous)
    
    ## A is always 1-dimensional
    
    family <- sort(c(nw$nodes[[A]]$idx,nw$nodes[[A]]$parents))

    ## didx and cidx are used as indices for A in the master
    didx    <- match(A,intersect(family,nw$discrete))
    didx    <- didx[!is.na(didx)]
    cidx    <- match(A,intersect(family,nw$continuous))
    cidx    <- cidx[!is.na(cidx)]
    
    if (nw$nodes[[A]]$type=="continuous") {
        cond <- list()
        
        if (!is.list(master$phi)) {
            cond[1] <- list(conditional.cont(cidx,
                                             master$mu,
                                             master$nu,
                                             master$rho,
                                             master$phi
                                             ))
        }
        else {
            for (i in 1:length(master$phi)) {
                
                cond[i] <- list(conditional.cont(cidx,
                                                 master$mu[i,],
                                                 master$nu[i],
                                                 master$rho[i],
                                                 master$phi[[i]]
                                                 ))
            }
        }
    }
    else if (nw$nodes[[A]]$type=="discrete") {
        
        cond <- list(list(alpha=master$alpha)) 
    }
    else stop("Wrong node type in conditional\n")
    
    cond
}


