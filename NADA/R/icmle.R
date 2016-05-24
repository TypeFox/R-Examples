####
# The following is code from the Icens package written by Vandal and Gentleman.
# We only use the EM function and it's related calls. 
####
#
# S-plus/R functions to determine the NPMLE of the distribution for interval-
# censored event time.
# Copyright 1998-2000 Alain Vandal and Robert Gentleman
# University of Auckland
# These functions should work, but their intent is to illustrate the
# concepts involved.  Functions are provided as is, with no guarantee.
# Redistribute freely, without undocumented modification & without charge.
# Queries to vandal@stat.auckland.ac.nz or rgentlem@stat.auckland.ac.nz.

# Returns a list of maximal antichains from a list of real valued intervals
# Arguments:  intvls:  n x 2 matrix;first column contains left endpoints
#				second column contains right endpoints
# Returned value:   list of length m (magnitude of underlying interval order)
#		      - each entry corresponds to one maximal antichain
#		      - each entry contains the row numbers of all intervals
#					  belonging to the maximal antichains
#		      - maximal antichains occur in the list in their natural
#					  linear ordering
# Known bugs: In R, will issue some ignorable warnings if there is
#right-censored data (with an "Inf" right endpoint).
Maclist <- function(intvls, Lopen=TRUE, Ropen=FALSE)
{
    m <- dim(intvls)[1]
    id <- 1:m
    or <- order(intvls[,1])
    maclist <- NULL
    curmac <- id[or[1]]
    minend <- intvls[curmac,2]
    for (i in 2:m) {
    	curintvl <- id[or[i]]
        if( intvls[curintvl,1]>minend ||
           ((Lopen || Ropen) &&  intvls[curintvl,1]==minend ) ) {
                                        # New maximal antichain
            maclist <- c(maclist,list(curmac))
            oldmac <- curmac
            curmac <- NULL
            for (j in 1:length(oldmac))
                if ( intvls[curintvl,1]<intvls[oldmac[j],2] ||
                 (!Lopen && !Ropen &&
                  intvls[curintvl,1]==intvls[oldmac[j],2]) )
                    curmac <- c(curmac,oldmac[j])
            curmac <- c(curmac,curintvl)
            minend <- min(intvls[curmac,2])
        } else {
            curmac <- c(curmac,curintvl)
            minend <- min(minend,intvls[curintvl,2]) }
    }
    c(maclist,list(curmac))
}

# Returns the clique matrix and Petrie pairs of an interval order
# given its list of maximal antichains Arguments: ml: list of maximal
# antichains as returned by Maclist Returned value: object containing
# # - pmat: clique matrix of the underlying interval order, # rows are
# ordered according to the linear ordering of # the maximal antichains
# # - ppairs: Petrie pairs indicate the first and last # maximal
# antichains to which each elements belongs

Macmat <- function(ml)
{
    temp <- NULL
    m <- length(ml)
    for (i in 1:m)
        temp <- c(temp,ml[[i]])
    temp <- sort(unique(temp))
    n <- length(temp)
    ppairs <- matrix(0,2,n)
    retmat <- matrix(0,m,n)
    for (i in 1:m) {
        for (j in ml[[i]]) {
            if (ppairs[1, j]==0)
                ppairs[1, j] <- i
            ppairs[2, j] <- i
        }
        retmat[i, ml[[i]]] <- 1
    }
    dimnames(ppairs) <- list(c("Start","End"),temp)
    dimnames(retmat) <- list(NULL,temp)
    ret <- list(pmat = retmat, ppairs = ppairs)
    class(ret) <- "petrie"
    return(ret)
}

# Produce the mapping of the maximal antichains to their real interval
# representation for an interval order given by real-valued intervals.
# Arguments:	intvls:	see Maclist
#		ml:	list of maximal antichains for the intervals as
#			returned by Maclist
# Returned values:  matrix m x 2 containing the mapping row-wise
#		(1rst row corresponds to 1rst maximal antichains, etc.)
#		     m is the number of maximal antichains,
#		     1rst column contains left endpoints of the mapping
#		     2nd column contains right endpoints of the mapping

MLEintvl <- function(intvls, ml=Maclist(intvls))
{
    if( ncol(intvls) != 2 || any(intvls[,2] < intvls[,1]) )
        stop("only one dimensional intervals can be handled")
    m <- length(ml)
    ret <- matrix(0, m, 2)
    for(i in 1:m) {
        LL <- min(intvls)
        RR <- max(intvls)
        for(j in 1:length(ml[[i]])) {
            LL <- max(LL, intvls[ml[[i]][j], 1])
            RR <- min(RR, intvls[ml[[i]][j], 2])
        }
        ret[i,  ] <- c(LL, RR)
    }
    ret
}

#############
#rescaleP is a function that rescales a prob vector so that elements
# that are negative or less than machine epsilon are set to zero.
###########
rescaleP <- function(pvec, tiny)
{
    pvec<-ifelse(pvec<tiny,0,pvec)
    pvec<-pvec/sum(pvec)
    return(pvec)
}

#The EM algorithm for interval censored data
EM<-function(A, pvec, maxiter = 500, tol = 1e-12)
{
    if( ncol(A)==2 && all(A[,2]>=A[,1]) ) {
        ml <- Maclist(A)
        intmap <- t(MLEintvl(A, ml))
        A <- Macmat(ml)$pmat
    }
    else
        intmap <- NULL
    i<-0
    notdone<-TRUE
    n<-ncol(A)
    Meps<-.Machine$double.eps
    if(missing(pvec))
        pvec <- apply(A, 1, sum)/sum(A)
    pvec<-rescaleP(pvec, Meps)
    while(i<maxiter && notdone) {
        i<-i+1
        dmat<-diag(pvec)
        t1<-dmat%*%A
        t2<-1/(t(A)%*%pvec)
        np<-rescaleP(as.vector(t1%*%t2)/n, Meps)
        if( sum(abs(np-pvec)) < tol )
            notdone<-FALSE
        pvec<-np
    }
    #if( notdone ) warning("EM may have failed to converge")
    ret <- list(pf=pvec, numiter=i,
                converge=!notdone, intmap=intmap)
    class(ret) <- "icsurv"
    return(ret)
}

