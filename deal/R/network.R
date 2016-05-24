## network.R
## Author          : Claus Dethlefsen
## Created On      : Fri Nov 02 21:20:16 2001
## Last Modified By: Claus Dethlefsen
## Last Modified On: Fri Jan 09 10:00:47 2004
## Update Count    : 319
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

network <- function(df,specifygraph=FALSE,inspectprob=FALSE,
                    doprob=TRUE,
                    yr=c(0,350),xr=yr) {
    ## creator for class 'network'
    ## df is a dataframe with one column per variable and one row per
    ## observation. Discrete variables are factors. We assume complete
    ## data, that is no NA's and at least one observation for each
    ## configuration of the factors.
    ##
    ## We create a 'trivial' network, which is a network without any arrows.
    
    if (length(dim(df))<1) stop("Can't handle networks with one node, sorry\n")
    
    nw   <- list()
    nw$n <- ncol(df)  ## df must have at least 2 columns...
    
    nw$discrete <- c()
    nw$continuous <- c()
    
    nw$nodes <- list()
    unit <- 2*pi/nw$n
    xc <- mean(xr)
    yc <- mean(yr)
    for (i in 1:nw$n) {
        pos <- c(cos( unit*i+pi/4),sin(unit*i+pi/4))*xc*.8 + c(xc,yc)
        ## create one node per column
        if (is.factor(df[,i])) {
            ## the node is discrete
            nw$nodes[[i]] <- node(i,c(),"discrete",
                                  names(df)[i],
                                  length(levels(df[,i])),
                                  levels(df[,i]),
                                  position=pos
                                  )
            nw$discrete <- c(nw$discrete,i)
        }
        else {
            ## the node is continuous
            nw$nodes[[i]] <- node(i,c(),"continuous",
                                  names(df)[i],
                                  position=pos
                                  )
            nw$continuous <- c(nw$continuous,i)
            
        }
    }
    
    nw$nd <- length(nw$discrete)
    nw$nc <- length(nw$continuous)
    stopifnot(nw$nd+nw$nc==nw$n) # invariant
    
    names(nw$nodes) <- names(df)
    
    class(nw) <- "network"
    
    if (specifygraph) {
        nw <- drawnetwork(nw,nocalc=TRUE)$nw
    }
    
    if (doprob) 
        nw <- prob(x=nw,df=df)
    
    if (inspectprob) nw <- inspectprob(nw)
    
    nw
}




print.network <- function(x,filename=NA,condposterior=FALSE,
                          condprior=FALSE,...) {
    nw <- x
    str <- paste("## ",nw$n,"(",nw$nd,"discrete+",nw$nc,") nodes;score=",
                 nw$score,";relscore=",nw$relscore,"\n")
    if (is.na(filename)) cat(str)
    else cat(str,file=filename)
    
    for (i in 1:nw$n)
        print(nw$nodes[[i]],filename=filename,condposterior,condprior)
    invisible(nw)
}

plot.network <- function(x,arrowlength=.25,
                         notext=FALSE,sscale=7,showban=TRUE,
                         yr=c(0,350),xr=yr
                         ,unitscale=20,cexscale=8,...) {
    
    nw <- x
    
    plot(0,0,xlim=xr,
         ylim=yr,type="n",
         axes=FALSE,xlab="",ylab="",...)
    
    unit <- 2*pi/nw$n
    xc <- mean(xr) # center coordinates
    yc <- mean(yr) # 
    
    ## show nodes
    for (i in 1:nw$n) 
        plot(nw$nodes[[i]],
             cexscale=cexscale,notext=notext,...)
    
    ## show score and relscore
    if (length(nw$score)>0 && !notext) {
        
        string <- paste("Score:",format(nw$score,2))
        if (length(nw$relscore)>0)
            string <- paste(string,"\n","Relscore:",format(nw$relscore,2))
        
        text(xc,0.97*yr[2],string)
    }
    
    ## show banlist
    if (showban) {
        if (!is.null(nw$banlist))
            if (nrow(nw$banlist)>0) {
                bl <- nw$banlist
                for (i in 1:nrow(bl)) {
                    from <- bl[i,2]
                    to   <- bl[i,1]
                    x  <- nw$nodes[[from]]$position 
                    y  <- nw$nodes[[to]]$position 
                    u <- (x - y) / sqrt(sum( (x-y)^2 )) 
                    
                    x <- x - u*unitscale
                    y <- y + u*unitscale
                    arrows( y[1],y[2],x[1],x[2],length=arrowlength,col="red",lty=2)
                } ## for
            } ## if (nrow...)
    } ## if (showban)
    
    
    ##< show arrows
    
    for (i in 1:nw$n) {
        ni <- nw$nodes[[i]]    # node i
        if (length(ni$parents)>0) {
            for (j in 1:length(ni$parents)) {
                x  <- ni$position # coords of ni
                pj <- ni$parents[j]  # parent j (index)
                y  <- nw$nodes[[pj]]$position # coords of pj
                
                u <- (x - y) / sqrt(sum( (x-y)^2 )) # unit vector from y to x
                
                x <- x - u*unitscale
                y <- y + u*unitscale
                
                arrows( y[1],y[2],x[1],x[2],length=arrowlength,...)
            }
        }
    }
    
}

score <- function(x,...) {
    UseMethod("score")
}

score.network <- function(x,...) {
    return(x$score)
}

score.node <- function(x,...) {
    return(x$loglik)
}


prob <- function(x,df,...) {
    UseMethod("prob")
}

prob.network <- function(x,df,...) {
    ## calculate initial probability
    x$nodes <- lapply(x$nodes,prob,df,x)
    x
}


banlist <- function(x) { x$banlist }

"banlist<-" <- function(x,value) {x$banlist <- value;x}

getnetwork <- function(x) x$nw

gettrylist <- function(x) x$trylist

gettable <- function(x) x$table

size <- function(x) x$n


