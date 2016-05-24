## addarrows.R
## Author          : Claus Dethlefsen
## Created On      : Fri Nov 02 21:02:07 2001
## Last Modified By: Claus Dethlefsen
## Last Modified On: Mon Jan 12 14:45:43 2004
## Update Count    : 197
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


addarrows <- function(nw, node, data, prior,trylist=vector("list",size(nw))) {
    ## Create all possible networks with arrows to/from node from/to
    ## nodes with smaller index.
    ##
    ## data:    dataframe
    ## prior:   jointprior
    ## returns  a list of networks (nwl) that have been learned
    ## trylist: a list of networks wherefrom some learning may be reused
    ##
    ## Used by: networkfamily
    ## Uses:    insert
    
    nwl  <- list(nw) # working network list
    
    for (i in 1:(node-1)) {
        for (j in 1:length(nwl)) {
            
            newnet <- insert(nwl[[j]],node,i,data,prior,trylist=trylist)
            if (length(newnet$nw) > 0 ) {   # Prevent NULL networks
                nwl <- c(nwl, list(newnet$nw))
                trylist <- newnet$trylist
            }
            newnet <- insert(nwl[[j]],i,node,data,prior,trylist=trylist)
            if (length(newnet$nw) > 0) {     # Prevent NULL networks
                nwl <- c(nwl, list(newnet$nw))
                trylist <- newnet$trylist
            }
        }
    }
    nwl <- nwl[-1]
    class(nwl) <- "networkfamily"
    list(nw=nwl,trylist=trylist)
}


insert <- function(nw,j,i,df,prior,nocalc=FALSE,
                   trylist=vector("list",size(nw))) {
    ## insert one arrow from node j to node i in network nw
    ## df: dataframe
    ## prior: jointprior
    ## nocalc: if F, relearn the net; else do not relearn
    ## trylist: a list of networks wherefrom some learning may be reused
    
    ## If arrow is illegal, returns a NULL network. Otherwise, returns a
    ## network with the arrow added (and relearned, if nocalc=F)
    
    ## Used by: addarrows, drawnetwork, addarrow, turnarrow
    ## Uses: learn(.network) if nocalc=F
    ## network attributes: nodes[[]]$type, nodes[[]]$parents,
    ##                     nw$banlist, nodes[[]]$tvar
    
    ## examines if the arrow is legal (no continuous parents for discrete
    ## node), is not banned.

    if (i==j) {
        ##        cat("Arrow (",i,"<-",j,") illegal\n")
        return(list(nw=NULL,trylist=trylist))  # RETURNS a NULL network
    }
    
    if (nw$nodes[[i]]$type=="discrete" &
        nw$nodes[[j]]$type=="continuous")
    {
        ##      cat("Arrow (",i,"<-",j,") illegal\n")
        return(list(nw=NULL,trylist=trylist))  # RETURNS a NULL network
    }
    else if (!is.na(match(j,nw$nodes[[i]]$parents))) {
        ##      cat("Arrow (",i,"<-",j,") already exists\n")
        return(list(nw=NULL,trylist=trylist))  # RETURNS a NULL network
    }
    else if (!is.na(match(i,nw$nodes[[j]]$parents))) {
        ##      cat("Arrow (",j,"<-",i,") already exists\n")
        return(list(nw=NULL,trylist=trylist))  # RETURNS a NULL network
    }
    else if (!is.null(nw$banlist)) {
        if (nrow(nw$banlist)>0) {
            idx <- (1:nrow(nw$banlist))[nw$banlist[,1]==j]
            if (length(idx)>0) 
                if (!is.na(match(i,nw$banlist[idx,2]))) {
                    ##    cat("Arrow (",j,"<-",i,") banned\n")
                    return(list(nw=NULL,trylist=trylist))  
                                        # RETURNS a NULL network
                }
        }
    }

    ## update parents
    nw$nodes[[i]]$parents <- sort(c(nw$nodes[[i]]$parents,j))
    if (!nocalc) {
        nw <- learn(nw,df,prior,i,trylist=trylist)
        trylist <- nw$trylist
        nw <- nw$nw
    }
    list(nw=nw,trylist=trylist)
}

remover <- function(nw,j,i,df,prior,nocalc=FALSE,
                    trylist=vector("list",size(nw))) {
    ## remove one arrow from node j to node i in network nw
    ## df: dataframe
    ## prior: jointprior
    ## nocalc: if F, relearn the net; else do not relearn
    ## trylist: a list of networks wherefrom some learning may be reused
    
    ## Used by: drawnetwork
    ## Uses: learn(.network) if nocalc=F
    ## network attributes: nodes[[]]$parents
    
    if (i==j) {
        ##    cat("Arrow (",i,"<-",j,") illegal\n")
        return(list(nw=NULL,trylist=trylist))  # RETURNS a NULL network
    }
    
    ## check if there *is* an arrow from i to j.
    parents <- nw$nodes[[i]]$parents
    if (!length(intersect(parents,j))>0) {
        cat("There's no arrow there!\n")
        return(list(nw=NULL,trylist=trylist))  # RETURNS a NULL network
    }
    else { 
        ## update parents
        nw$nodes[[i]]$parents <- setdiff(nw$nodes[[i]]$parents,j)
    }
    if (!nocalc) { nw <- learn(nw,df,prior,i,trylist=trylist)
                   trylist <- nw$trylist
                   nw <- nw$nw
               }
    list(nw=nw,trylist=trylist)
}
