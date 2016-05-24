## jointcont.R
## Author          : Claus Dethlefsen
## Created On      : Wed Mar 06 12:52:57 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Sun Jul 27 15:57:54 2003
## Update Count    : 333
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

jointcont <- function(nw,timetrace=FALSE) {
    ## From the continuous part of nw, the joint distribution is
    ## determined from the local distributions in nodes$prob.
    ##
    ## If eg. x|y,z, y|z, z are given, the joint distribution of x,y,z
    ## is returned
    ##
    
    if (timetrace) {t1 <- proc.time();cat("[jointcont ")}
    
    ## First, determine the discrete nodes and their dimensions
    Dim <- c()
    TD <- 1
    if (nw$nd>0) {
        for (i in nw$discrete) {
            Dim <- c(Dim, nw$nodes[[i]]$levels)
        }
        TD <- prod(Dim)
    }

    ## create labels for the configurations of the discrete variables
    lablist <- c()
    if (nw$nd>0) {
        for (i in 1:TD) {
            cf <- findex( i, Dim, FALSE)
            label <- ""
            for (j in 1:ncol(cf)) {
                label <- paste(label, nw$nodes[[nw$discrete[j]]]$levelnames[cf[1,j]]
                               ,sep=":")
            }
            lablist <- c(lablist,label)
        }
    }
    
    
    ## determine the continuous nodes
    lab <- c()
    for (i in nw$continuous)
        lab <- c(lab,nw$nodes[[i]]$name)
    
    mu     <- matrix(0,TD,nw$nc) 
    sigma2 <- matrix(0,nw$nc,nw$nc)
    sigma2list <- list()
    colnames(mu) <- lab
    rownames(mu) <- lablist
    rownames(sigma2) <- colnames(sigma2) <- lab
    for (i in 1:TD) sigma2list[[i]] <- sigma2
    names(sigma2list) <- lablist
    
    calclist <- c()
    allnodes <- c(nw$continuous)
    
    nidx <- 0
    while ( length( setdiff(allnodes,calclist) )>0 ) {
        ## the main loop. Evaluates nodes sequentially so that the
        ## parents of the current node has already been evaluated
        
        nidx <- nidx%%(nw$nc)+1
        nid  <- nw$continuous[nidx]
        
        if ( length(intersect(nid,calclist))>0) {
            next
        }
        
        node    <- nw$nodes[[nid]] 
        Pn      <- node$prob        ## the local distribution
        parents <- node$parents     ## the parents, 
        if (nw$nc>0)    cparents<- sort(intersect(parents,nw$continuous))
        else cparents <- c()
        if (nw$nd>0)    dparents<- sort(intersect(parents,nw$discrete))
        else dparents <- c()
        
        if ( length( setdiff(cparents,calclist) ) > 0  ) {
            next
        }
        
        
        ## calculate unconditional mu, sigma2 from node|parents
        if (!length(cparents)>0) {
            M <- array(1:TD,dim=Dim)
            if (length(dparents)>0) {
                
                mdim <- c()
                for (i in dparents) 
                    mdim <- c(mdim,nw$nodes[[i]]$levels)
                m <- array(1:TD,dim=mdim) 
                
                ## inflate
                ## first, permute Dim appropriately
                ivek <- c(match(dparents,nw$discrete),
                          match(setdiff(nw$discrete,dparents),nw$discrete))
                jDim <- Dim[ivek]
                bigM <- array(m,jDim)

                ## permute back
                permvek <- match(1:nw$nd,ivek)
                bigM <- aperm(bigM, permvek)
                for (i in 1:length(unique(c(bigM)))) { 
                    theidx <- M[bigM==i]
                    cf <- findex(theidx,Dim,config=FALSE)
                    cfm<- cf[,match(dparents,nw$discrete)]
                    cfm <- matrix(cfm,nrow=length(theidx))
                    theidxm <- findex(cfm,mdim,config=TRUE)
                    paridx  <- match(1:nw$nc,c(nid,cparents))
                    for (k in 1:length(theidx)) {
                        mu[theidx,nidx] <- Pn[theidxm[k],2]
                        sigma2list[[theidx[k]]][nidx,nidx] <- Pn[theidxm[k],1]
                    }
                }
            }
            else { ## no discrete parents
                for (i in 1:TD) {
                    mu[i,nidx] <- Pn[2]
                    sigma2list[[i]][nidx,nidx] <- Pn[1]
                }
            } ## end else (no discrete parents)
            
        }
        else { # we have continuous (and possibly discrete) parents
            
            for (k in 1:TD) {
                if (length(dparents)>0) {
                    mdim <- c()
                    for (i in dparents) 
                        mdim <- c(mdim,nw$nodes[[i]]$levels)
                    
                    Mcf <- findex(k,Dim,config=FALSE)
                    didx <- match(dparents,nw$discrete)
                    dcf <- Mcf[,didx]

                    if (length(dcf)==2) 
                        dcf <- matrix(dcf,ncol=2)
                    
                    kidx <- findex(dcf,mdim,config=TRUE)
                }
                else
                    kidx <- 1
                ## parentidx: index in mu,sigma2list of parents
                ## calcidx:   index in mu,sigma2list of processed nodes
                parentidx <- match(cparents,nw$continuous)
                calcidx <- match(sort(calclist),nw$continuous)
                if (!length(dparents)>0) {        
                    m.ylx <- Pn[2]
                    s2.ylx<- Pn[1]
                    b.ylx <- Pn[3:length(Pn)]
                }
                else {
                    m.ylx <- Pn[kidx,2]
                    s2.ylx<- Pn[kidx,1]
                    b.ylx <- Pn[kidx,3:ncol(Pn)]
                }
                m.x   <- mu[k,parentidx]
                s2.x  <- sigma2list[[k]][parentidx,parentidx]

                pid <- match(parentidx,sort(calclist))
                pid <- pid[!is.na(pid)]
                b.calc <- rep(0,length(calcidx))
                b.calc[pid] <- b.ylx
                s2.calc <- sigma2list[[k]][calcidx,calcidx]

                s.xycalc <- s2.calc %*% b.calc
                
                s.xy  <- s2.x %*% b.ylx
                s2.y  <- s2.ylx + c(s.xy)%*%b.ylx
                
                m.y   <- m.ylx + b.ylx%*%m.x
                
                mu[k,nidx] <- m.y 
                
                sigma2list[[k]][nidx,nidx] <- s2.y
                sigma2list[[k]][calcidx,nidx] <- s.xycalc
                sigma2list[[k]][nidx,calcidx] <- t(s.xycalc)
            }
        }
        
        
        calclist <- c(calclist,nid)
        
    } ## while
    
    if (timetrace) {
        t2 <- proc.time()
        cat((t2-t1)[1],"]")
    }
    
    list(mu=mu,sigma2=sigma2list)
} ## function discjoint

