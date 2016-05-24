## networkfamily.R
## Author          : Claus Dethlefsen
## Created On      : Tue Oct 30 16:43:05 2001
## Last Modified By: Claus Dethlefsen
## Last Modified On: Wed Jul 28 09:39:43 2004
## Update Count    : 429
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

networkfamily <- function(data,nw=network(data),prior=jointprior(nw),
                          trylist=vector("list",size(nw)),timetrace=TRUE) {
    ## Creator class for networkfamily
    ##
    ## Generates all possible networks with the restriction that
    ## discrete nodes cannot have continuous parents. (see insert)
    
    ## Uses: numbermixed, addarrows, learn.network, cycletest
    ## and attributes of nw: nd,nc
    
    ## Value: 
    ## networklist: A list of network-objects
    ## trylist: an updated trylist

    if (timetrace) {t1 <- proc.time();cat("[networkfamily ")}

    nw <- learn(nw,data,prior,trylist=trylist)
    trylist <- nw$trylist
    nw <- nw$nw
    
    ndiscrete  <- nw$nd
    ncontinuous<- nw$nc
    
    cat("Creating all (",
        numbermixed(ndiscrete,ncontinuous),
        " minus restrictions) networks with ",ndiscrete," discrete and ",
        ncontinuous," continuous nodes\n",sep="")
    
    nwl     <- list() # network list
    n       <- ndiscrete + ncontinuous
    
    nwl     <- list(nw)  # current network list
    for (node in 2:n) {
        for (idx in 1:length(nwl)) {
            nws <- addarrows(nwl[[idx]],node,data,prior,trylist=trylist)
            trylist <- nws$trylist
            nwl <- c(nwl,nws$nw)
        }
    }
    
    cat("Created",length(nwl),"networks, ")
    
    if (ndiscrete>2|ncontinuous>2) {
        cat("removing cycles...\n")
        
        nwlres <- nwl[!unlist(lapply(nwl,cycletest))]
        
        cat(length(nwl)-length(nwlres),"cycles removed, ending up with",length(nwlres),"networks\n")
    }
    else nwlres <- nwl
    class(nwlres) <- "networkfamily"
    
    if (timetrace) {
        t2 <- proc.time()
        cat((t2-t1)[1],"]\n")
    }
    
    list(nw=nwlres,trylist=trylist)
}



plot.networkfamily <- function(x,
                               layout=rep(min(1+floor(sqrt(length(x))),5),2),
                               cexscale=5,arrowlength=0.1,
                               sscale=7,...) {
    nwf <- x
    par(mfrow=layout)
    for (i in 1:length(nwf)) {
        par(mar=c(0,0,0,0))
        plot(nwf[[i]],cexscale=cexscale,arrowlength=arrowlength,sscale=sscale,showban=FALSE,...)
    }
    par(mfrow=c(1,1))
}

nwfsort <- function(nwf) {
    ## sort according to network score, and add relative scores
  
    n <- length(nwf)
    ## first, create a vector with the indices and scores
    tab <- rep(NA,n)
    for (i in 1:n)
        tab[i] <- nwf[[i]]$score
    
    ## then find the sort list of indices
    sl <- sort.list(-tab)
    
  relscore <- exp(tab - tab[sl[1]]) 
  
    ## create the sorted family
    nwf <- nwf[sl]
    for (i in 1:n)
        nwf[[i]]$relscore <- relscore[sl[i]]
    class(nwf) <- "networkfamily"
    nwf
}

print.networkfamily <- function(x,...) {
    
    nwf <- nwfsort(x) ## ensure they are sorted
    
    g <- function(x) x$name
    nw <- nwf[[1]]
    
    cat("Discrete:  ")
    if (nw$nd>0) {
        nn <- nw$discrete[1]
        cat(nw$nodes[[nn]]$name,"(",nw$nodes[[nn]]$levels,")",sep="")
        if (nw$nd>1) {
            for (i in nw$discrete[-1])
                cat(",",nw$nodes[[i]]$name,"(",nw$nodes[[i]]$levels,")",sep="")
        }
        cat("\n")
    }
    else cat("\n")
    
    
    cat("Continuous:")
    if (nw$nc>0) {
        nn <- nw$continuous[1]
        cat(nw$nodes[[nn]]$name,sep="")
        if (nw$nc>1) {
            for (i in nw$continuous[-1])
                cat(",",nw$nodes[[i]]$name,sep="")
        }
        cat("\n")
    }
    else cat("\n")
    
    cat("  log(Score)\t|Relscore\t|Network\n")
    printline()
    for (i in 1:length(nwf)) {
        nw <- nwf[[i]]
        cat(i,". ",nw$score,"\t",nw$relscore,sep="")
        if (i==1) cat("\t")
        cat("\t",sep="")
        for (j in 1:nw$n) {
            nd <- nw$nodes[[j]]
            cat("[",nd$name,sep="")
            if (length(nd$parents)>0) {
                cat("|",
                unlist(lapply(nw$nodes[nd$parents],g)),
                sep="")
            }
            cat("]")
        }
        cat("\n")
        
    } ## for
    invisible(nwf)
}

