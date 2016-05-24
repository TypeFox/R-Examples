## node.R
## Author          : Claus Dethlefsen
## Created On      : Fri Nov 02 21:18:50 2001
## Last Modified By: Claus Dethlefsen
## Last Modified On: Wed Jul 28 09:39:55 2004
## Update Count    : 410
## Status          : OK
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

nodes <- function(nw) nw$nodes
"nodes<-" <- function(nw,value) {nw$nodes<-value;nw}

localprob <- function(nw) lapply(nodes(nw),function(node) node$prob)
"localprob<-" <- function(nw,name,value) {
    names(value) <- names(nodes(nw)[[name]]$prob)
    nodes(nw)[[name]]$prob <- value
    nw
}

localprior <- function(node) node$condprior
localposterior <- function(node) node$condposterior

node <- function(idx,parents,type="discrete",name=paste(idx),
                 levels=2,levelnames=paste(1:levels), position=c(0,0)) {
    ## creator for class 'node'
    
    ## idx:       The unique index of the node
    ## name:      The plotted name
    ## parents:   Vector with indices of parents
    ## type:      "discrete" or "continuous"
    ## levels:    If discrete, the number of levels
    ## levelnames:If discrete, the printed names of the levels
    
    nd         <- list()
    nd$idx     <- idx
    nd$parents <- parents
    nd$type    <- type
    nd$name    <- name
    nd$position<- position
    if (type=="discrete") {
        nd$levels     <- levels
        nd$levelnames <- levelnames
    }
    
    class(nd)  <- "node"
    
    nd
}

print.node <- function(x,filename=NA,condposterior=TRUE,condprior=TRUE,...) {
    
    nd <- x
    str <- paste(nd$idx,nd$name,nd$type,sep="\t")
    str <- paste(str,"(",nd$levels,")",sep="")
    for (i in 1:length(nd$parents)) {
        if (length(nd$parents)>0)
            str <- paste(str,nd$parents[i],sep="\t")
    }
    if (is.na(filename)) cat(str,"\n")
    else cat(str,"\n",file=filename,append=TRUE)
    
    if (condprior)   {
      printline()
        cat("Conditional Prior:",nd$name)
        if (length(nd$parents)>0) {
            cat("| ")
            for (j in 1:length(nd$parents))
                cat(nd$parents[j]," ")
        }
        cat("\n")
        print(nd$condprior)
    }
    
    
    if (condposterior)   {
      printline()
        cat("Conditional Posterior:",nd$name)
        if (length(nd$parents)>0) {
            cat("| ")
            for (j in 1:length(nd$parents))
                cat(nd$parents[j]," ")
        }
        cat("\n")
        print(nd$condposterior)
    }
    
    invisible(nd)
}

plot.node <- function(x,cexscale=10,notext=FALSE,...) {
    
    if (x$type=="discrete") {tt <- 19;col <- "white"} 
    else {tt <- 21;col <- "black"}
    
    points(x$position[1],x$position[2],cex=cexscale,pch=tt,...)
    if (!notext) text(x$position[1],x$position[2],x$name,col=col,...)
}


prob.node <- function(x,df,nw,...) {
    
    data <- df
    node <- x # for compatibility reasons.
    
    ## node: current node
    ## nw: The network - we need the parents of the node
    ## data: for continuous nodes, we need to estimate mu and sigma2
    ## from data. For discrete nodes we need to count the number of
    ## cases for each state.
    ##
    ## Returns: a node with the prob-attribute set to
    ##          for discrete: an array of dimension equal to the levels
    ##          of the discrete parents and value:
    ##          if equalcases=T  1/xx, where xx is
    ##                           the product of the levels.
    
    
    nodelist <- nw$nodes
    
    if (node$type=="discrete") {
        
        vek <- rep(NA,length(node$parents)+1)
        vek[1] <- node$levels
        dnames <- list(node$levelnames)
        if (length(node$parents)>0) {
            for (i in 1:length(node$parents)) { 
                vek[i+1] <- nodelist[[node$parents[i]]]$levels
                dnames <- c(dnames,
                            list(nodelist[[node$parents[i]]]$levelnames))
            }
        }
        node$prob <- array(1/prod(vek),dim=vek)
        dimnames(node$prob) <- dnames
        if (length(node$parents)>0)
            node$prob <- prop.table(node$prob,2:(length(node$parents)+1))
            
    } ## type=="discrete"
    
    if (node$type=="continuous") {
        ## for each product level of discrete parents, calculate
        ## mean and variance from the data.
        
        if (length(node$parents)>0) {
            
            
            parents   <- sort(node$parents)
            if (nw$nd>0)    dparents<- sort(intersect(parents,nw$discrete))
            else dparents <- c()
            if (nw$nc>0)    cparents<- sort(intersect(parents,nw$continuous))
            
            if (length(cparents)>0) {
                
                if (length(dparents)>0) {
                    ## at least one discrete and one continuous parent
                    ## cat("The true mixed case\n")
                    
                    ## find configurations of discrete variables
                    ## for each configuration
                    ##     reduce data
                    ##     do a regression on the cont.parents
                    
                    Dim <- c()
                    dnames <- list()
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
                            label <- paste(label,
                                    nw$nodes[[dparents[j]]]$levelnames[cf[1,j]]
                                           ,sep=":")
                        }
                        lvek <- c(lvek,label)
                    }
                    
                    M <- matrix(NA,TD,2+length(cparents))
                    rownames(M) <- lvek
                    colnames(M) <- c("s2",paste("Intercept",node$name,sep=":"),
                                     names(data)[cparents])
                    
                    for (i in 1:TD) {
                        config <- findex(i,Dim,config=FALSE)
                        obs <- data[,c(dparents,cparents,node$idx)]
                        for (k in 1:ncol(config)) {
                            j <- config[1,k]
                            ## reduce data
                            lev <- nw$nodes[[dparents[k]]]$levelnames[j]
                            obs <- obs[obs[,k]==lev,]
                        }
                        
                        X <- obs[,(length(dparents)+1):(ncol(obs)-1)]
                        y <- obs[,ncol(obs)]
                        lsobj <- lsfit(X,y)
                        
                        beta <- coef(lsobj)
                        s2   <- sum(resid(lsobj)^2)/nrow(data)
                        
                        M[i,] <- c(s2,beta)
                        
                    }
                    
                    node$prob <- M
                }
                else {
                    ## only continuous parents
                    X <- data[,cparents]
                    y <- data[,node$idx]
                    lsobj <- lsfit(X,y)
                    
                    beta <- coef(lsobj)
                    s2   <- sum(resid(lsobj)^2)/nrow(data)
                    
                    node$prob <- c(s2,beta)
                    names(node$prob) <- c("s2",
                                          paste("Intercept",node$name,sep=":")
                                          ,names(data)[cparents])
                }
            }
            else { ## only discrete parents
                
                Dim <- c()
                dnames <- list()
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
                        label <- paste(label,
                                   nw$nodes[[dparents[j]]]$levelnames[cf[1,j]]
                                       ,sep=":")
                    }
                    lvek <- c(lvek,label)
                }
                
                M <- matrix(NA,TD,2)
                rownames(M) <- lvek
                colnames(M) <- c("s2",paste("Intercept",node$name,sep=":"))
                
                for (i in 1:TD) {
                    ## Find configuration of discrete parents
                    ## Find the data that fits
                    ## mean,var of these variables
                    ##     if no data: mean=0, var=100
                    config <- findex(i,Dim,config=FALSE)
                    
                    obs <- data[,c(dparents,node$idx)]
                    for (k in 1:ncol(config)) {
                        j <- config[1,k]
                        ## reduce data
                        lev <- nw$nodes[[dparents[k]]]$levelnames[j]
                        obs <- obs[obs[,k]==lev,]
                    }
                    if (nrow(obs)>1) {
                        n <- nrow(obs)
                        M[i,] <- c(var(obs[,ncol(obs)])*(n-1)/n,
                                   mean(obs[,ncol(obs)]))
                    }
                    else {
                        M[i,] <- c(100,0)
                        if (nrow(obs)==1)
                            M[i,2] <- obs[1,ncol(obs)]
                    } ## else
                } ## for
                node$prob <- M
            } ## else
        } ## if parents
        else { ## no parents
            n <- dim(data)[1]
            node$prob <- c(var(data[,node$idx])*(n-1)/n,mean(data[,node$idx]))
            names(node$prob) <- c("s2",paste("Intercept",node$name,sep=":"))
            
        }
    } ## type=="continuous"
    
    node
} ## function: prob.node

cond.node <- function(node,nw,nw.prior=jointprior(nw)) {
    ## make conditional prior for this node and attach it
    
    thismaster <- localmaster(sort(c(node$idx,node$parents)),
                              nw,nw.prior)
    if (length(node$parents)>0) { ## parents are present
        thiscond <- conditional(node$idx,thismaster,nw)
        
        if (node$type=="continuous") {
            contparents <- intersect(node$parents,nw$continuous)
            if (length(contparents)<1) { ## no cont. parents
                for (k in 1:length(thiscond)) {
                    thiscond[[k]]$tau <- thismaster$nu[k]
                    thiscond[[k]]$mu  <- thismaster$mu[k]
                    thiscond[[k]]$phi <- thismaster$phi[[k]]
                    thiscond[[k]]$rho <- thismaster$rho[k]
                }
            }
        }
    }
    else { ## no parents, so thiscond is just the master
        thiscond <- list(thismaster)
        thiscond[[1]]$tau <- thismaster$nu
    }
    
    ##  node$master  <- thismaster ## only used for debugging
    node$condprior    <- thiscond
    node
}

