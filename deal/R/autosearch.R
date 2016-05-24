## autosearch.R
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 11 10:54:00 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Thu Dec 04 12:43:15 2008
## Update Count    : 307
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

autosearch <- function(initnw,data,prior=jointprior(network(data)),maxiter=50,
           trylist= vector("list",size(initnw)),trace=TRUE,
           timetrace=TRUE,showban=FALSE,removecycles=FALSE)
{
    ## Greedy search
    
    ## initnw: initial network with conditionals calculated
    ##
    ## output: networklist: a sorted list of all tried networks.
    
    ## used by: heuristic.
    ## uses: addarrow,removearrow,turnarrow,nwfsort,cycletest
    ##       initnw$score
    
    ## Algorithm:
    ## Create all networks with one arrow added (addarrow)
    ## Create all networks with one arrow turned (turnarrow)
    ## Create all networks with one arrow removed (removearrow)
    ## Calculated scores for all networks
    ## Choose the non-cyclic network that increases the score the most,
    ## or stop. 
    
    if (timetrace) {t1 <- proc.time();cat("[Autosearch ")
                    tadd <- 0
                    trem <- 0
                    ttur <- 0
                    tsor <- 0
                    tcho <- 0
                }
    
    
    nw <- initnw

    model <- modelstring(initnw)
    score <- initnw$score
    
    slut <- FALSE
    it   <- 0
    hiscore <- initnw$score

    while (!slut & it < maxiter) {
        it <- it + 1
        
        if (timetrace) {s1 <- proc.time()[1]}
#         cat("adding arrows\n")
        thisnwl.add <- addarrow(nw,data,prior,trylist=trylist)
        trylist     <- thisnwl.add$trylist
        thisnwl.add <- thisnwl.add$nw
        if (timetrace) {s2 <- proc.time()[1];
                        tadd <- tadd+s2-s1
                    }
#         cat("removing arrows\n")
        thisnwl.rem <- removearrow(nw,data,prior,trylist=trylist)
        trylist <- thisnwl.rem$trylist
        thisnwl.rem <- thisnwl.rem$nw
        if (timetrace) {s3 <- proc.time()[1];
                        trem <- trem+s3-s2
                    }
#         cat("turning arrows\n")
        thisnwl.tur <- turnarrow(nw,data,prior,trylist=trylist)
        trylist <- thisnwl.tur$trylist
        thisnwl.tur <- thisnwl.tur$nw
        if (timetrace) {s4 <- proc.time()[1];
                        ttur <- ttur+s4-s3
                    }
        thisnwl <- c(thisnwl.add,thisnwl.rem,thisnwl.tur)
        class(thisnwl) <- "networkfamily"

        thisnwl <- nwfsort(thisnwl)
        if (timetrace) {s5 <- proc.time()[1];
                        tsor <- tsor+s5-s4
                    }
        
        
        ## remove cycles and then choose the best
        if (removecycles)
        {
            thisnwl <- thisnwl[!unlist(lapply(thisnwl,cycletest))]
            nwcand <- thisnwl[[1]]
        ## what if all of them contains cycles? They do not.
        }
        else
        {
            ## choose the 'best' and then check for cycle.
            kk <- 1
            while (TRUE) {
                nwcand <- thisnwl[[kk]]
                kk <- kk + 1
                if (!cycletest(nwcand)) break
                if (timetrace) cat(".")
            }
        }
        if (timetrace) {s6 <- proc.time()[1];
                        tcho <- tcho+s6-s5
                    }
        
        model <- c(model,unlist(lapply(thisnwl,modelstring)))
        score <- c(score,unlist(lapply(thisnwl,function(x) x$score)))
        
        
        if (nwcand$score > hiscore) {
            hiscore <- nwcand$score
            nw <- nwcand
            if (trace) {plot(nw,showban=showban)
                    }
            cat("(",it,") ",hiscore," ",modelstring(nw),"\n",sep="")
        }
        else
        {
            slut <- TRUE
        }
        
    } ## end while
    
    if (timetrace) {
        t2 <- proc.time()
        total <- (t2-t1)[1]
        cat("Total",total,"add",tadd,"rem",trem,"turn",ttur,"sort",tsor,"choose",tcho,"rest",total-tadd-trem-ttur-tsor-tcho,"]\n")
        
    }
    
    table <- cbind(model,score)
    table <- table[sort.list(-as.numeric(table[,2])),]
    list(nw=learn(nw,data,prior)$nw,table=table,trylist=trylist)
}

modelstring <- function(x) {
    res <- ""
    g <- function(x) x$name
    for (j in 1:x$n) {
        nd <- x$nodes[[j]]
        res <- paste(res,"[",nd$name,sep="")
        if (length(nd$parents)>0) {
            res <- paste(res,"|",
                               paste(unlist(lapply(x$nodes[nd$parents],g)),
                                     collapse=":"),
                               sep="")
        }
        res <- paste(res,"]",sep="")
    }
        res
}

makenw <- function(tb,template) {
    res <- apply(tb,1,as.network,template)
    class(res) <- "networkfamily"
    nwfsort(res)
}

as.network <- function(nwstring,template) {
    x <- nwstring
    ## x: vector of (modelstring and score)
    ## from 'modelstring' (output from modelstring), create a network
    ## structure (not learned!)
    ## template is a network with the same nodes
    ## Thus, the function inserts the parent-relations that are
    ## described in mstr.
    ## as.network(modelstring(x),x) is the identity
    ## function. Beware though, that the output network needs to be
    ## learned so that the parameters are correct.

    ## first, split into nodes, assuming the correct form
    ## [node1|parent1:parent2][node2][node3|parent1]

    mstr <- x[1]
    score<- x[2]
    
    st <- strsplit(strsplit(mstr,"\\[")[[1]],"\\]")

    ## now, we have a list 2:nw$n+1 with all the nodes
    nw <- template
    for (i in 1:nw$n) {
        cn <- st[[i+1]]
        ## does this node have parents?
        cns <- strsplit(cn,"\\|")[[1]]
        if (length(cns)>1) {
            ## yes, parents are present
            parents <- cns[-1]
            parstr <- strsplit(parents,":")[[1]]
            pidx <- match(parstr,names(nw$nodes))
            pidx <- pidx[!is.na(pidx)]
            nw$nodes[[i]]$parents <- sort(pidx)
        }
        else
            nw$nodes[[i]]$parents <- c()
    }
    nw$score <- as.numeric(score)
    nw
}


addarrow <- function(nw,df,prior,trylist=vector("list",size(nw))) {
    ## Create all networks with one extra arrow
    ## return list of networks (nwl) (Possibly NULL)
    ## trylist: a list of networks wherefrom some learning may be reused
    
    ## used by: autosearch
    ## uses: insert
    ## and network attributes: n
    
    nwl <- list()
    n <- nw$n
    try <- cbind(1:n,rep(1:n,rep(n,n)))
    
    for (i in 1:nrow(try)) {
        newnet <- insert(nw,try[i,1],try[i,2],df,prior,
                         trylist=trylist)
        
        if ( !is.null(newnet$nw) ) { # prevent NULL networks
            nwl[length(nwl)+1] <- list(newnet$nw)
            trylist <- newnet$trylist
        }
        
    }
    class(nwl) <- "networkfamily"
    list(nw=nwl,trylist=trylist)
}



removearrow <- function(nw,df,prior,trylist=vector("list",size(nw))) {
    ## create all networks with one arrow less
    ## return list of networks (possibly NULL)
    ## trylist: a list of networks wherefrom some learning may be reused
    
    ## used by: autosearch
    ## uses: insert, learn
    ## and network attributes: n, nodes$parents
    nwl <- list()
    for (i in 1:nw$n) {
        if (length(nw$nodes[[i]]$parents) > 0) {
            for (j in 1:length(nw$nodes[[i]]$parents)) {
                newnet <- nw
                newnet$nodes[[i]]$parents <- newnet$nodes[[i]]$parents[-j]
                newnet <- learn(newnet,df,prior,i,trylist=trylist)
                trylist <- newnet$trylist
                newnet <- newnet$nw
                nwl[length(nwl)+1] <- list(newnet)
            }
        }
    }
    class(nwl) <- "networkfamily"
    list(nw=nwl,trylist=trylist)
}

turnarrow <- function(nw,df,prior,trylist=vector("list",size(nw))) {
    ## create all networks with one arrow turned
    ## return list of networks (possibly NULL)
    ## trylist: a list of networks wherefrom some learning may be reused
    
    ## used by: autosearch
    ## uses: insert, learn
    ## and network attributes: n, nodes$parents
    
    nwl <- list()
    for (i in 1:nw$n) {
        if (length(nw$nodes[[i]]$parents) > 0) {
            for (j in 1:length(nw$nodes[[i]]$parents)) {
                newnet <- nw
                parent <- nw$nodes[[i]]$parents[j]
                newnet$nodes[[i]]$parents <- newnet$nodes[[i]]$parents[-j]
                newnet <- learn(newnet,df,prior,i,trylist=trylist)
                trylist<- newnet$trylist
                newnet <- newnet$nw
                newnet <- insert(newnet,i,parent,df,prior,trylist=trylist) #parent is learned here
                trylist <- newnet$trylist
                newnet  <- newnet$nw
                if (length(newnet) > 0) { # prevent NULL networks
                    nwl[length(nwl)+1] <- list(newnet) 
                }
            }
        }
    }
    class(nwl) <- "networkfamily"
    list(nw=nwl,trylist=trylist)
}
