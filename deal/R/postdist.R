## postnw.R --- 
## Author          : Claus Dethlefsen
## Created On      : Sat Sep 28 17:15:47 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Thu Jul 24 15:21:36 2003
## Update Count    : 17
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

postdist <- function(nw) {
    ## calculate means of parameters and overwrite the prob attributes
    ## of the nodes
    
    nw$nodes <- lapply(nw$nodes,postdist.node,nw)
    nw
}

postdist.node <- function(nd,nw,vtype="mode") {
    ## calc. local prob from post.parameters (in cond.posterior)
    if (nd$type=="discrete") {
        if (length(nd$parents)>0) {
            a <- nd$condposterior[[1]]$alpha
            npa <- length(dim(a))
            as<- apply(a,2:npa,sum)
            bs<- sweep(a,2:npa,as,"/")
            nd$prob <- bs
        }
        else {
                    nd$prob <- nd$condposterior[[1]]$alpha/
            sum(nd$condposterior[[1]]$alpha)
                }
    }
    if (nd$type=="continuous") {
        dpar <- intersect(nd$parents,nw$discrete)
        cpar <- intersect(nd$parents,nw$continuous)
        Dim <- c()                
        for (i in dpar) {
            Dim <- c(Dim, nw$nodes[[i]]$levels)
        }
        TD <- prod(Dim)

        res <- matrix(NA,nrow=0,ncol=(2+length(cpar)))
        for (i in 1:TD) {
            cp <- nd$condposterior[[i]]
            mu <- cp$mu

            if (vtype=="mean") {
                ## mean
                s2 <- cp$phi/(cp$rho-2)
            }
            if (vtype=="mode") {
                ## mode
                s2 <- cp$phi/(cp$rho+2)
            }

            res <- rbind(res,c(s2,mu))
        }
        nd$prob <- res
    }
        
    
    nd
}
