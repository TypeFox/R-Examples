## cycletest.R
## Author          : Claus Dethlefsen
## Created On      : Fri Dec 21 14:04:58 2001
## Last Modified By: Claus Dethlefsen
## Last Modified On: Sun Sep 15 08:05:24 2002
## Update Count    : 59
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

cycletest <- function(nw) {
    ## Does nw contain a cycle?
    ##
    ## Algorithm:
    ##  if nw$n == 1   return(F)
    ##  else res <- findleaf(nw)  ## res=0 if no leaf else idx of leaf
    ##       if res==0 return(T)
    ##       else nw <- (nw with node idx deleted)
    ##            cycletest(nw)
    ##
    ## Uses: findleaf
    ## and network attributes:
    ##      n,nodes
    ## Used by: networkfamily, drawnetwork, autosearch,
    ##          addrandomarrow, turnrandomarrow
    
    if (nw$n == 1) { #cat("only one node\n");
        return(FALSE)}
    else {
        res <- findleaf(nw)
        if (res == 0) {
            ## cat("No leaf found\n");
            return(TRUE)
        }
        else {
            ## cat("deleting node: ",nw$nodes[[res]]$name,"\n")
            nw$nodes <- nw$nodes[-res]
            nw$n     <- nw$n - 1
            ## should update cont and disc, but I won't. Just be careful
            ## how you use the procedure!
            cycletest(nw)
        }
    }
}

findleaf <- function(nw) {
    ## find a node not being a parent to any other node
    ##
    ## Uses network attributes: n, nodes
    ##   and  node  attributes: idx, parents
    ##
    ## Used by: cycletest
    jump <- FALSE
    for (i in 1:nw$n) {    ## for each node
        for (j in 1:nw$n) {  ## testing i against i (hmm)
            ## cat("Is",nw$nodes[[i]]$name,"parent to",nw$nodes[[j]]$name,"?")
            
            ## is i a parent to j?
            ## Here, it is necessary to use 'idx', since we have been
            ## deleting nodes. Thus the indices are no longer 1:n
            res <- match(nw$nodes[[i]]$idx, nw$nodes[[j]]$parents)
            if (!is.na(res)) { ## i is not a leaf
                jump <- TRUE
                break ## next i
            }
        }
        if (!jump)  return(i)
        jump <- FALSE
    }
    ## did not find any
    res <- 0
    res
}

