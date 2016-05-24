# Do error checking and preprocssing of the variance list that was
#  given as an argument to coxme.
coxme.varcheck <- function(ncluster, varlist, n, gvars, groups, sparse,
                           rescale, pdcheck) {
    if (is.null(varlist)) {
	varlist <- vector('list', ncluster)
	names(varlist) <- gvars
	}
    else {
	if (is.matrix(varlist)) {
	    if (ncluster >1) 
		    stop("Matrix given in 'varlist' goes with which term?")
	    varlist <- list(varlist)
	    names(varlist) <- gvars
	    }
	else if (is.function(varlist)) {
	    if (ncluster >1) 
		    stop("Matrix given in 'varlist' goes with which term?")
	    varlist <- list(varlist)
	    names(varlist) <- gvars
	    }
	else if (is.list(varlist)) {
	    if (ncluster==1 && !is.list(varlist[[1]])) {
		# the list need not be named elements
		varlist <- list(varlist)
		}
            names(varlist) <- gvars
	    }
	else stop("Illegal varlist, must be a list or matrix")
	}

    vindex <- match(names(varlist), gvars)  #for each variance matrix, the term
    if (length(vindex) && any(is.na(vindex)))
        stop(paste("varlist element", (names(varlist))[is.na(vindex)],
                   "not found in the list of random effects"))
    # put the varlist in the same order as gvars
    vindex <- match(gvars, names(varlist), nomatch=0) +1
    varlist <- (c(list(NULL), varlist))[vindex]

    kindex <- matrix(0, nrow=n, ncol=ncluster)  # final group indicators
    ntheta <- rep(1, ncluster)  #total number of variances to estimate/term
    # Go through the grouping variables one term at at time
    for (i in 1:ncluster) {
        if (ncluster==1) tgrp <- groups
        else             tgrp <- groups[[i]]
        tgrp <- tgrp[, drop=T]     #throw away any unused levels
        gnames <- levels(tgrp)
        nfrail <- length(gnames)
        if (is.null(varlist[[i]])) {
            # No varlist was provided for this term    
            # At this point, I assume that the random term was 1| something
            #   Eventually, we will figure out the pdMat stuff of lme, and use
            #   it as the 'prototype' for construction
            # Usually, we want to construct a sparse bdsmatrix form of the
            #   identity matrix as our variance term.  If the number of random
            #   effects is small then we can use a dense matrix -- there is
            #   no need to trade accuracy for speed in the Cox fits.
            # Otherwise, don't use sparseness for any term that accounts
            #  for more than sparse[2]% of the population.  If there were 50
            #  groups, I'd expect each to account for about 2% of the total. 
            #  This is again a "saftey" net for the Newton-Raphson.
            indx <- as.numeric(tgrp)
            if (nfrail <  sparse[1]) {  #non-sparse
                kmat <- bdsmatrix(blocksize=nfrail, 
                                  blocks = rep(0, (nfrail*(nfrail+1))/2),
                                  dimnames=list(gnames, gnames))
                diag(kmat) <- rep(1, nfrail)
                }
            else {
                temp <- table(tgrp)/length(tgrp)
                large <- (temp > sparse[2])
                if (any(large)) {
                    nlarge  <- sum(large)    #these terms are not sparse
                    nsparse <- sum(!large)
                    temp <- c(gnames[!large], gnames[large])
                    indx <- match(tgrp, temp)  
                    gnames <- temp
		    kmat <- bdsmatrix(blocksize=rep(1, nsparse), 
                                     blocks=rep(1., nsparse),
                                     rmat=rbind(matrix(0., nsparse, nlarge),
                                                diag(nlarge)),
                                     dimnames=list(temp,temp))
                    }
                else kmat <- bdsI(gnames)  # all sparse
                }
            # Add the component to the end of the variance list
            varlist[[i]] <- list(kmat)
            kindex[,i] <- indx
            }
      
        else {
            # The caller supplied a variance matrix or list for this
            #  term. 
            tlist <- varlist[[i]]
            tlist <- bdsmatrix.reconcile(tlist, levels(tgrp))
            if (is.list(tlist)) {  #complex variance structure
		# It will be a list, not a list of lists
                ntheta[i] <- length(tlist)
                if (rescale) {
                    for (j in 1:ntheta[i]) {
			kmat <- tlist[[j]]
                        # scale the kinship matrix to have a diagonal of 1's.  
                        #  This may already have been done by the user
                        # With inbreeding, the diagonal might not be constant.
                        # Until we figure out the right thing to do in 
                        #  that case, we complain
			if (pdcheck) {
			    temp <- gchol(kmat)
			    if (any(diag(temp) <0)) 
				 stop("A variance matrix is not positive defininte")
			    }
                        temp <- diag(kmat)
#                        if (any(temp != temp[1])) 
#                        warning("Diagonal of variance matrix is not constant")
			if (max(temp)==0)
			    stop("Diagonal of a variance matrix is 0!")
                        if (max(temp) !=1) {
                            kmat <- kmat/max(temp)
                            tlist[[j]] <- kmat
                            }
                        }
		    }
		else if (pdcheck) {
		    for (kmat in tlist) {
			temp <- gchol(kmat)
			if (any(diag(temp) < 0))	
			    stop("A variance matrix is not non-negative definite")
			}
                    }
                kindex[,i] <-match(tgrp, (dimnames(tlist[[1]]))[[1]])
                varlist[[i]] <- tlist                   
                }
            else {  # tlist is just a matrix
		if (pdcheck) tempg <- gchol(tlist)
		if (rescale) {
		    temp <- diag(tlist)
#		    if (any(temp != temp[1])) 
#			warning("Diagonal of kmat is not constant")
		    if (temp[1] !=1)  tlist <- tlist/temp[1]
		    if (pdcheck &&any(diag(tempg) <0))
		      stop("A variance matrix is not non-negative definite")
		    }
		else if (pdcheck && any(diag(tempg) <0))
		      stop("A variance matrix is not positive definite")
                kindex[,i] <- match(tgrp, (dimnames(tlist))[[1]])
                varlist[[i]] <- list(tlist)
                }
            }
        }
    list(varlist=varlist, kindex=kindex, ntheta=ntheta)
    }
