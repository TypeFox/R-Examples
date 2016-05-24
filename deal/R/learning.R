## learning.R
## Author          : Claus Dethlefsen
## Created On      : Mon Jan 14 12:24:13 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Mon Jan 12 14:32:52 2004
## Update Count    : 551
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

learn <- function(nw, df, prior=jointprior(nw),
                  nodelist=1:size(nw),trylist=
                  vector("list",size(nw)),
                  timetrace=FALSE
                  ) {
    ## nw: network to be learned (condprior must be present in the nodes)
    ## df: dataframe with observations
    ## nodelist: vector of node-indices of nodes to be learned (default
    ##                          is to learn all nodes) 
    ## trylist: a list of networks wherefrom some learning may be reused
    ##
    ## Returns a network with the following attributes
    ##       score: calculated (or updated) network-score
    ##       for each node in nodelist:
    ##           loglik: the log-likelihood contribution of the node
    ##           cond:   updated posterior parameters
    ##
    ## Uses: cond, learnnode
    ## and network attributes: nodes, score is updated
    ## and node attributes: condprior,condposterior is updated
    ##
    ## Used by: insert,remover,removearrow,turnarrow,
    ##          manualsearch,networkfamily,
    ##          turnrandomarrow,deleterandomarrow (perturb)
    
    
    if (timetrace) {t1 <- proc.time();cat("[Learn.network ")}
    
    old <- df 
    
    for (i in nodelist) {
        node <- nw$nodes[[i]]
        
        ## use trylist
        if (!is.null(trylist[[node$idx]]))  {
            cur <- paste(node$parents,collapse=":")
            curm <- match(cur,trylist[[node$idx]][,1])
            if (!is.na(curm)) {
                nw$nodes[[i]]$loglik <-
                    as.numeric(trylist[[node$idx]][curm,2])
                break
            }
        }
        
        ## learning
        node <- cond.node(node,nw,prior) ## master prior procedure
        
        node$condposterior <- node$condprior ## reset posterior
        node$loglik        <- 0
        node <- learnnode(node,nw,df,timetrace=FALSE)## learn!
        
        ## update trylist
        streng <- paste(node$parents,collapse=":")
        tal    <- node$loglik
        if (is.null(trylist[[i]])) {
            trylist[[i]] <- cbind(streng,tal)
        }
        else
            trylist[[i]] <- rbind(trylist[[i]],cbind(streng,tal))

        ## update network
        nw$nodes[[i]] <- node
    }
    
    ## calculate network score
    nw$score <- 0
    for (i in 1:nw$n) 
        nw$score <- nw$score + nw$nodes[[i]]$loglik
    
    if (timetrace) {
        t2 <- proc.time()
        cat((t2-t1)[1],"]")
    }
    list(nw=nw,trylist=trylist)
}

learnnode <- function(node,nw,df,prior=jointprior(nw),timetrace=FALSE) {
    ## node: node to be learned. condprior must be present
    ## nw:   network
    ## df:   dataframe to learn from
    ##
    ## Returns: node with extra (or updated) attributes:
    ##          loglik: the loglikelihood contribution from this node
    ##          cond:   the posterior parameters
    ##
    ## Uses: udisclik,postc0c,postcc
    ## And network attributes: nc,nd,continuous,discrete,nodes
    ## And node attributes: type,condprior,condposterior
    ## (updated),loglik (updated),parents,levels,idx
    ##
    ## Used by: learn.network
    
    if (timetrace) {t1 <- proc.time();cat("[Learn.node ")}
    
    ## discrete nodes:
    if (node$type=="discrete") {
        
        node$condposterior[[1]]$alpha <- node$condprior[[1]]$alpha+
            as.array(table(df[,sort(c(node$idx,node$parents))]))
        node$loglik <- udisclik(node,nw,df) ## batch update likelihood term
        node <- postdist.node(node,nw)
        if (timetrace) {
            t2 <- proc.time()
            cat((t2-t1)[1],"]")
        }
        return(node)
    }
    
    ## continuous nodes:
    
    ## 0 parents
    if (!length(node$parents)>0) {
        res <- postc0c(node$condposterior[[1]]$mu,
                       node$condposterior[[1]]$tau,
                       node$condposterior[[1]]$rho,
                       node$condposterior[[1]]$phi,
                       df[,node$idx])
        ## Alternatively, use this (pure R)
        ##
        ##        res <- post0(node$condposterior[[1]]$mu,
        ##                       node$condposterior[[1]]$tau,
        ##                       node$condposterior[[1]]$rho,
        ##                       node$condposterior[[1]]$phi,
        ##                       df[,node$idx])
        node$condposterior[[1]]$mu <- res$mu
        node$condposterior[[1]]$tau <- res$tau
        node$condposterior[[1]]$rho <- res$rho
        node$condposterior[[1]]$phi <- res$phi
        node$loglik <- res$loglik
        node <- postdist.node(node,nw)
        return(node)
    }
    parents <- node$parents     
    if (nw$nc>0)    cparents<- sort(intersect(parents,nw$continuous))
    else cparents <- c()
    if (nw$nd>0)    dparents<- sort(intersect(parents,nw$discrete))
    else dparents <- c()
    
    if (length(dparents)>0& (!length(cparents)>0)) {
        ##        cat("Discrete parents, no Cont. parents\n")
        ##        cat("dparents=",dparents,"\n")
        
        mscore <- 0
        Dim <- c()
        for (i in dparents)
            Dim <- c(Dim,nw$nodes[[i]]$levels)
        for (j in 1:prod(Dim)) {
            cf <- findex(j,Dim,config=FALSE)
 
            idx <- 1:nrow(df)
            for (k in 1:length(dparents)) {
                pcf <- nw$nodes[[dparents[k]]]$levelnames[cf[1,k]]
                idx <- idx[df[idx,dparents[k]]==pcf]
            } ## for k

            if (length(idx)>0) {
                mu  <- node$condposterior[[j]]$mu
                tau <- node$condposterior[[j]]$tau
                rho <- node$condposterior[[j]]$rho
                phi <- node$condposterior[[j]]$phi
                y   <- df[idx,node$idx]
                
                res <- postc0c(mu, tau, rho, phi, y)
                ## Alternative (pure R):
                ## res <- post0(mu, tau, rho, phi, y)
                node$condposterior[[j]]$mu <- res$mu
                node$condposterior[[j]]$tau <- res$tau
                node$condposterior[[j]]$rho <- res$rho
                node$condposterior[[j]]$phi <- res$phi
                mscore  <- mscore + res$loglik
            }
        } ## for j
        node$loglik <- mscore
        node <- postdist.node(node,nw)
        return(node)
    }
    
    if (!length(dparents)>0&length(cparents)>0) {
        ##        cat("Continuous parents\n")
        res <- postcc(node$condposterior[[1]]$mu,
                      node$condposterior[[1]]$tau,
                      node$condposterior[[1]]$rho,
                      node$condposterior[[1]]$phi,
                      df[,node$idx],
                      cbind(1,df[,cparents]))
        ## Alternative (pure R):
#        res <- post(node$condposterior[[1]]$mu,
#                      node$condposterior[[1]]$tau,
#                      node$condposterior[[1]]$rho,
#                      node$condposterior[[1]]$phi,
#                      df[,node$idx],
#                      cbind(1,df[,cparents]))
        
        node$condposterior[[1]]$mu <- res$mu
        node$condposterior[[1]]$tau <- res$tau
        node$condposterior[[1]]$rho <- res$rho
        node$condposterior[[1]]$phi <- res$phi
        node$loglik <- res$loglik
        node <- postdist.node(node,nw)
        return(node)
    }
    
    if (length(dparents)>0&length(cparents)>0) {
        ##    cat("Mixed parents\n")

        mscore <- 0
        Dim <- c()
        for (i in dparents)
            Dim <- c(Dim,nw$nodes[[i]]$levels)
        for (j in 1:prod(Dim)) {
            cf <- findex(j,Dim,config=FALSE)
            
            idx <- 1:nrow(df)
            for (k in 1:length(dparents)) {
                pcf <- nw$nodes[[dparents[k]]]$levelnames[cf[1,k]]
                
                idx <- idx[df[idx,dparents[k]]==pcf]
            } ## for k
            if (length(idx)>0) {
                mu <- node$condposterior[[j]]$mu
                tau <- node$condposterior[[j]]$tau
                rho <- node$condposterior[[j]]$rho
                phi <- node$condposterior[[j]]$phi
                y   <- df[idx,node$idx]
                z   <- cbind(1,df[idx,cparents])

                res <- postcc(mu, tau, rho, phi, y, z)
                ## Alternative (pure R):
                ## res <- post(mu, tau, rho, phi, y, z)
                node$condposterior[[j]]$mu <- res$mu
                node$condposterior[[j]]$tau <- res$tau
                node$condposterior[[j]]$rho <- res$rho
                node$condposterior[[j]]$phi <- res$phi
                mscore  <- mscore + res$loglik
            }
            
            
        } ## for j
        node$loglik <- mscore
        node <- postdist.node(node,nw)
        return(node)
        
    }
    
    
}


udisclik <- function(node,nw,df) {
    ## update likelihood term for the discrete nodes

    alpha  <- node$condposterior[[1]]$alpha
    cprior <- node$condprior[[1]]$alpha
    n <- sum(cprior) # img.db size
    N <- sum(alpha)  # n+#obs
    nobs <- N-n
    
    if (length(node$parents)>0) {
        ## we have parents!
        
        idx <- sort(c(node$idx,node$parents))
        cidx <- 1:length(idx)
        pidx <- cidx[-match(node$idx,idx)]
        ## alpha_{+d|i_pa(d)}
        ##      alphaj <- table(cprior,pidx)
        alphaj <- apply(cprior,pidx,sum)
        ## alpha_{+d|i_pa(d)}+n_{+d|i_pa(d)}
        condj <- alphaj + as.array(table(df[,node$parents]))
        
        ##        tres <- prod(gamma(condj)/gamma(alphaj))
        logtres <- -sum( lgamma(condj) - lgamma(alphaj) )
        
        ##      res[[i]] <- tres * prod(gamma(alpha)/gamma(cprior))
        res   <- logtres + sum( lgamma(alpha) - lgamma(cprior) )
        
    }## if parents
    else { ## no parents
        ##      res[[i]] <- prod(gamma(alpha)/gamma(cprior))*gamma(n)/gamma(N)
        res <- sum( lgamma(alpha) - lgamma(cprior)) + lgamma(n)-lgamma(N)
    }
    ##    res[[i]] <- log(res[[i]])
    res
}


