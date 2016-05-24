## master.R
## Author          : Claus Dethlefsen
## Created On      : Thu Nov 29 21:28:29 2001
## Last Modified By: Claus Dethlefsen
## Last Modified On: Wed Jul 23 19:22:41 2003
## Update Count    : 299
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

localmaster <- function(family,nw,prior=jointprior(nw)) {
  
    ## family: indices of a subset of nodes in the network 'nw'
    ## jointprior: jointprior(nw,N)
    ##
    ## Returns: the joint local master prior for the family
    
    
    listsum <- function(liste,idx=1:nrow(liste[[1]])) {
        ## sum elements of list containing a matrix as each element
        ## narrow down to liste[[i]][idx,idx] (always made to be a matrix)
        
        res <- matrix(0,
                      nrow(as.matrix(liste[[1]][idx,idx])),
                      ncol(as.matrix(liste[[1]][idx,idx])))
        
        for (i in 1:length(liste)) 
            res <- res + as.matrix(liste[[i]][idx,idx])
        
        res
    }
    
    
    
    ## determine indices of discrete and cont. nodes
    didx <- match(family,nw$discrete)
    didx <- didx[!is.na(didx)]
    cidx <- match(family,nw$continuous)
    cidx <- cidx[!is.na(cidx)]
    
    ## initialize
    alpha <- NA
    nu    <- NA
    rho   <- NA
    mu    <- NA
    phi   <- NA
    
    if (!length(cidx)>=1) { ## no cont. nodes
        alpha <- apply(prior$jointalpha,didx,sum)
    }
    else if(!length(didx)>=1) { ## no disc. nodes
        
        nu <- sum(prior$jointnu)
        rho<- sum(prior$jointrho)
        
        M <- as.matrix(prior$jointmu[,cidx]*c(prior$jointnu))
        if (nrow(prior$jointmu)==1)
            dim(M) <- c(1,length(prior$jointmu[,cidx]))
        
        mu <- apply( M ,2,sum )/nu
        
        ss <- matrix(0,length(cidx),length(cidx))
        for (i in 1:nrow(prior$jointmu)) {
            thismu <- as.matrix(prior$jointmu[i,cidx])
            mumean <- as.matrix(mu)
            
            ss <- ss+prior$jointnu[i]*(thismu-mumean)%*%t(thismu-mumean)
        }
        
        phi<- listsum(prior$jointphi,cidx)+ss
    }
    
    else { ## mixed
        nu    <- apply(prior$jointnu   ,didx, sum)
        rho   <- apply(prior$jointrho  ,didx, sum)
        nconfig <- length(nu) # number of configs.
        mu    <- matrix(0,nconfig,length(cidx))
        phi    <- list()
        for (i in 1:nconfig) phi[[i]] <- matrix(0,length(cidx),length(cidx))
        
        ## find dimension from  levels of discrete nodes
        D <- c()
        for (i in 1:length(didx)) {
            D <- c(D,nw$nodes[[nw$discrete[didx[i]]]]$levels)
        }
        jmu <- prior$jointmu

        for (i in 1:nrow(jmu)) {
            ## the corresp. configuration of the disc. variables in the
            ## joint distribution
            idx <- findex(i,dim(prior$jointalpha),config=FALSE)
            y   <- findex(matrix(idx[didx],1),D,config=TRUE)
            mu[y,] <- mu[y,] + jmu[i,cidx]*prior$jointnu[i]
            phi[[y]][,] <- phi[[y]][,] +
                prior$jointphi[[i]][cidx,cidx]
        }
        for (i in 1:nrow(mu)) 
            mu[i,] <- mu[i,]/nu[i]
      
        ## adjust phi with sum(nu_j(mu_j-mean(mu))(mu_j-mean(mu))^t)
        for (i in 1:nrow(jmu)) {
            idx <- findex(i,dim(prior$jointalpha),config=FALSE)
            y   <- findex(matrix(idx[didx],1),D,config=TRUE)
            phi[[y]] <- phi[[y]] +
                prior$jointnu[i]*(jmu[i,cidx]-mu[y,])%*%t(jmu[i,cidx]-mu[y,])
            rownames(phi[[y]]) <- colnames(phi[[y]])
        }
        colnames(mu) <- colnames(phi[[1]])
    }
  
    list(alpha=alpha,
         nu=nu,
         rho=rho,
         mu=mu,
         phi=phi)
}



