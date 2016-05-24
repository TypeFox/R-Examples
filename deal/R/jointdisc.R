## jointdisc.R --- 
## Author          : Claus Dethlefsen
## Created On      : Wed Mar 06 12:52:57 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Mon Dec 15 12:05:40 2008
## Update Count    : 31
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

jointdisc <- function(nw,timetrace=FALSE) {
    ## From the discrete part of nw, the joint distribution is
    ## determined from the local distributions in nodes$prob.
    ##
    ## If eg. A|B,C, B|C, C are given, the joint distribution of A,B,C
    ## is returned
    ##
    
    if (timetrace) {t1 <- proc.time();cat("[jointdisc ")}
    
    ## First, determine the discrete nodes and their dimensions
    
    Dim <- c()
    lablist <- list()
    for (i in nw$discrete) {
        Dim <- c(Dim, nw$nodes[[i]]$levels)
        lablist <- c(lablist,list(nw$nodes[[i]]$levelnames))
    }
    
    ## Dim is the dimension of the returned joint distribution
    jointprob <- array(1,Dim)
    dimnames(jointprob) <- lablist
    
    ## for each node, multiply jointprob by the local distribution
    ## (blown up appropriately).
    
    for (nid in nw$discrete) {
        node    <- nw$nodes[[nid]] 
        Pn      <- node$prob        ## the local distribution
        parents <- node$parents     ## the parents, 
        if (nw$nd>0)    dparents<- sort(intersect(parents,nw$discrete))
        else dparents <- c()

        idx <- c(node$idx, dparents) ## sequence in Pn
        pidx<- 1:length(idx)         ## corresponding idx
        jidx<- 1:nw$nd               ## idx in jointprior
    
        ## dimension of c(node,parents)
        nDim <- c(node$levels)
        for (i in dparents) 
            nDim <- c(nDim,nw$nodes[[i]]$levels)
        
        ## blow up
        ## first, permute Dim appropriately
        ivek <- c(pidx,setdiff(jidx,pidx))
#        ivek <- c(idx,setdiff(jidx,idx)) # changed 25/6-2007 due to
        ## Jean-Baptiste DENIS 

        jDim <- Dim[ivek]
        bigPn <- array(Pn,jDim)
        ## permute indices appropriately
        permvek <- match(1:nw$nd,ivek)
        bigPn <- aperm(bigPn, permvek)
        
        jointprob <- jointprob * bigPn
    } ## for
    
    if (timetrace) {
        t2 <- proc.time()
        cat((t2-t1)[1],"]")
    }
    jointprob
} ## function discjoint
  
