#
#
# The two variance functions for the attributable risk function.
#
attribrisk_jack <- function(attribrisk, match.id, X, Y, weights, 
                          offset, xbase, xgroup, casewt) {
    nxval <- max(xgroup)
    ar.xfold <- double(nxval)

    # This function is one of the few in R called with data sets that have
    #  actual case weights.  In that case we mimic leaving one out by
    #  reducing the weights by 1.
    if (casewt)  {
        for(i in 1:nrow(X)){
            temp <- weights
            temp[i] <- temp[i] -1
            if (is.null(match.id))
                ar.xfold[i] <- attribrisk.fit(X, Y, w=temp, offset,
                                              match = NULL, 
                                              xbase= xbase, fit=FALSE)
            else 
                ar.xfold[i] <- attribrisk.fit(X, Y, w=temp, offset,
                                              match= match.id,
                                              xbase= xbase, fit=FALSE)
        }
        avar <- sum(weights*(ar.xfold - attribrisk)^2)
    } 
    else {
        for(i in 1:nxval){
            keep <- which(xgroup != i)
            if (is.null(match.id))
                ar.xfold[i] <- attribrisk.fit(X[keep,,drop=F], Y[keep], 
                                              weights[keep], offset[keep], 
                                              match = NULL, 
                                              xbase= xbase[keep,,drop=F], 
                                              fit=FALSE)
            else 
                ar.xfold[i] <- attribrisk.fit(X[keep,,drop=F], Y[keep], 
                                          weights[keep], offset[keep], 
                                          match= match.id[keep], 
                                          xbase = xbase[keep,,drop=F],
                                          fit=F)
        }
        avar <- sum((ar.xfold - attribrisk)^2)
    }
    avar
}


# Variance function based on bootstrap

attribrisk_boot <- function(match.id, X, Y, weights, 
                          offset, xbase, casewt, xgroup, conf, nboot, bci) {
    
    if (casewt) {
        # We expand out the data to full size
        indx <- rep(1:nrow(X), weights)
        X <- X[indx,,drop=FALSE]
        Y <- Y[indx]
        weights <- weights[indx]
        offset <- offset[indx]
        xbase <- xbase[indx,, drop=FALSE]
        weights <- rep(1, nrow(X))
        xgroup <- seq(along=weights)
    }


    if (is.null(match.id)) {
        # logistic regression case
        stat <- function(d, i, myargs){
            attribrisk.fit(myargs$X[i,,drop=F], myargs$Y[i], myargs$weights[i],
                                myargs$offset[i], 
                                match =NULL, 
                                xbase= myargs$xbase[i,,drop=F], fit=FALSE)
         }
    }
    else {
        n <- nrow(X)
        mlist <- split(1:n, match.id)
        nper <- table(match.id)  
        #
        # The trick here is that we need to create a new data set
        #  where i indexes the matched sets (usually pairs but not always)
        #  rather than the subjects
        #  and treats each such set as a new one.
        # Say that match.id were (1,1, 2, 3,3,3, 4,4), n=8,  and that the stat
        #  function below were called with i= 1, 3, 3,2, 4
        # The row selected will be {(1,2), (4,5,6), (4,5,6), 2, (7,8)}
        #  the strata for the call should be 1,1, 2,2,2, 3,3,3, 4, 5,5
        #  It needs to look like new independent samples of the sets
        #
        stat <- function(d, i, myargs){
            j <- unlist(mlist[i])
            k <- rep(1:length(i), nper[i])
            attribrisk.fit(myargs$X[j,,drop=F], myargs$Y[j], myargs$weights[j],
                                myargs$offset[j], 
                                match = k, 
                                xbase = myargs$xbase[j,,drop=F], fit=FALSE)
                           
        }
    }
    fit.ar <- boot(data=unique(xgroup), stat, R=nboot, stype='i', 
                   myargs = list(attribrisk=attribrisk, match.id=match.id,
                                 X=X, Y=Y, weights=weights, 
                                 offset=offset, xbase=xbase))
    b.ci <- do.call(boot.ci, c(list(boot.out= fit.ar, conf=conf), bci))

    list(var=var(fit.ar$t), boot=fit.ar, boot.ci=b.ci)
}
