## makesimprob.R
## Author          : Claus Dethlefsen
## Created On      : Tue Feb 26 13:25:44 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Thu Dec 07 08:59:15 2006
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

makesimprob <- function(nw,
                        s2=function(idx,cf) {
                            cf <- as.vector(cf)
                            xs <- (1:length(cf))
                            log(xs%*%cf+1)
                        },
                        m0=function(idx,cf) {
                            cf <- as.vector(cf)
                            xs <- (1:length(cf))^2
                            .69*(xs%*%cf)
                        },
                        m1=function(idx,cf) {
                            cf <- as.vector(cf)
                            xs <- (1:length(cf))*10
                            idx*(cf%*%xs)
                        }) {
    ## sets up (and asks user for) probablities to simulate from
    ##
    ## Idea: let s2 and m1 depend on the node-index and on j
    ## Perhaps passing functions as arguments?
    ##
    ## Discrete variables are organised as follows
    ## The table always has the node itself as the first one. The
    ## remaining (conditioning) are sorted according to their index. We
    ## let the probabilities be equal.
    
    for (nid in 1:nw$n) {
        
        node <- nw$nodes[[nid]]
        parents <- node$parents
        if (nw$nd>0)    dparents<- sort(intersect(parents,nw$discrete))
        else dparents <- c()
        if (nw$nc>0)    cparents<- sort(intersect(parents,nw$continuous))
        
        if (length(dparents)>0) {
            Dim <- c()
            dnames <- list(node$levelnames)
            for (i in dparents) {
                Dim <- c(Dim,nw$nodes[[i]]$levels)
                dnames <- c(dnames,list(nw$nodes[[i]]$levelnames))
            }
            TD <- prod(Dim)
            
            ## create labels
            lvek <- c()
            for (i in 1:TD) {
                cf <- findex( i, Dim, FALSE)
                label <- ""
                for (j in 1:ncol(cf)) {
                    label <- paste(label, nw$nodes[[dparents[j]]]$levelnames[cf[1,j]]
                                   ,sep=":")
                }
                lvek <- c(lvek,label)
            }
        }
        else {
            dnames <- list(node$levelnames)
            TD  <- 1
            Dim <- c()
        }
    
        if (node$type=="continuous") {
            M <- matrix(NA,TD,1+1+length(cparents))
            
            if (length(dparents)>0) rownames(M) <- lvek
            
            colnames(M) <- c("s2","m0",names(nw$nodes[cparents]))
            
            for (it in 1:nrow(M)) {
                ifelse(TD>1,cf <- findex( it, Dim, FALSE), cf <- 1)        
                M[it,1] <- s2(nid,cf)
                M[it,2] <- m0(nid,cf)
                if (length(cparents)>0) {
                    for (itt in 1:length(cparents))
                        M[it,3:(2+itt)] <- m1(nid,cf)
                }
            }
            
            nw$nodes[[nid]]$simprob <- M
        }
        else if (node$type=="discrete") {

            Dim <- c(node$levels,Dim)
            simtab <- array(1/prod(Dim),dim=Dim)
            dimnames(simtab) <- dnames
            if (length(node$parents)>0)
                simtab <- prop.table(simtab,2:(length(node$parents)+1))
            
            nw$nodes[[nid]]$simprob <- simtab
        }
        else stop("makesimprob: Type is wrong")
    }
    
    nw
}
