# gibbs.msbvar.R Functions for Gibbs sampling an MSBVAR model based on
# a SZ prior.

# Patrick T. Brandt
# 20081113 : Initial version
# 20100325 : Added the computation for the log marginal data densities
# 20100615 : Cleaned up permutation / flipping code.  This now has its
#            own function that is called inside the Gibbs sampler.
#            Also reorganized the code in this file so that is easier
#            to read, document and follow.
# 20100617 : Wrote separate functions for the Beta, Sigma, and e block
#            updates.  Cleans up the code for the gibbs.msbvar into
#            managable chunks.
# 20120113 : Replaced filtering-sampler steps for the state space with
#            compiled Fortran code.
# 20140609 : Added gc() calls for some improved memory usage


####################################################################
# Utility functions for managing matrices in the Gibbs sampler
####################################################################
#
# vech function for efficiently sorting and subsetting the unique
# elements of symmetric matrices (like covariance Sigma)
#
vech <- function (x)
{
    x <- as.matrix(x)
    if (dim(x)[1] != dim(x)[2]) {
        stop("Non-square matrix passed to vech().\n")
    }
    output <- x[lower.tri(x, diag = TRUE)]
    dim(output) <- NULL
    return(output)
}

####################################################################
# xpnd function for reconstructing symmetric matrices from vech
####################################################################

xpnd <- function (x, nrow = NULL)
{
    dim(x) <- NULL
    if (is.null(nrow))
        nrow <- (-1 + sqrt(1 + 8 * length(x)))/2
    output <- matrix(0, nrow, nrow)
    output[lower.tri(output, diag = TRUE)] <- x
    hold <- output
    hold[upper.tri(hold, diag = TRUE)] <- 0
    output <- output + t(hold)
    return(output)
}

####################################################################
# tcholpivot function for computing Choleski's of known PD matrices
# with correction for underflow via pivoting.  Necessary to handle
# some numerically unstable possible draws that the Gibbs sampler can
# take with low probability.  Returns the draws after undoing the
# pivoting and transposing for correct usage in scaling draws from the
# relevant MVN densities.
####################################################################
tcholpivot <- function(x)
{
    tmp <- chol(x, pivot=TRUE)
    return(t(tmp))
}


####################################################################
# Sampler blocks for the MSBVAR Gibbs sampler
####################################################################
# Metropolis-Hastings sampler for the transition matrix Q
Q.drawMH <- function(Q, SStrans, prior, h)
{
    Q.old <- Q

    # compute the density ordinate based on the counts.
    alpha <- SStrans + prior - matrix(1, h, h)

    Q.new <- t(apply(alpha, 1, rdirichlet, n=1))

    eta.new <- eta.ergodic(Q.new,h)
    eta.old <- eta.ergodic(Q.old,h)

    # acceptance rate
    # assume initial state is state 1
    A <- eta.new[1] / eta.old[1]

    if (runif(1) < A) {
        # need to add acceptance ratio calculations here
        # e.g., accsum <- accsum + 1
        return(Q.new)
    } else {
        return(Q.old)
    }
}

# Calculate invariant probability distribution, eta,
# when assumed distribution for initial state is ergodic
# Fruhwirth-Schnatter 2006, page 306
eta.ergodic <- function(Q, h)
{
    A <- rbind(diag(h)-Q, 1)
    b <- solve(t(A) %*% A) %*% t(A)
    eta <- b[,h+1]  # h by 1 column vector
    return(eta)
}


# Gibbs sampler for the transition matrix Q
Q.drawGibbs <- function(Q, SStrans, prior, h)
{
    # compute the density ordinate based on the counts.
    alpha <- SStrans + prior - matrix(1, h, h)
    # Sample
    Q.new <- t(apply(alpha, 1, rdirichlet, n=1))
    return(Q.new)
}


####################################################################
# State-space sampling for the regimes -- this is the FFBS
# implementation.  Actual work is done by calling compiled Fortran
# code in the package
####################################################################
## oldSS.ffbs <- function(e, bigt, m, p, h, sig2, Q)
## {
##     TT <- bigt-p
##     SStmp <- .Fortran("FFBS",
##                       bigt=as.integer(bigt),
##                       m = as.integer(m), p = as.integer(p),
##                       h = as.integer(h), e = e, sig2 = sig2,
##                       Q = Q, f = double(1),
##                       filtprSt = matrix(0,as.integer(TT),as.integer(h)),
##                       SS = matrix(as.integer(0),
##                                   as.integer(TT),as.integer(h)),
##                       transmat = matrix(as.integer(0),
##                                         as.integer(h),as.integer(h))
##                       )

##     # Output
##     ss <- list(SS=SStmp$SS, transitions=SStmp$transmat)
##     return(ss)
## }

SS.ffbs <- function(e, bigt, m, p, h, sig2, Q)
{

    # Hamilton's forward filter

    ff <- .Fortran('ForwardFilter',
                   nvar=as.integer(m), e=e,
                   bigK=as.integer(h), bigT=as.integer(bigt),
                   nbeta=as.integer(1+m*p),
                   sig2, Q,
                   llh=double(1),
                   pfilt=matrix(0,bigt+1,h))

    # Backwards multi-move sampler

    # Draw random Unif(0,1) to feed into backwards sampler.
    # Doing this on the R side (rather than using RAND() in Fortran)
    # allows to control the seed.
    rvu <- runif(bigt+1)

    bs <- .Fortran('BackwardSampler',
                   bigK=as.integer(h),bigT=as.integer(bigt),
                   ff$pfilt,Q,rvu,
                   backsamp.bigS=integer(bigt+1),
                   transmat=matrix(as.integer(0),h,h))

    # Now, convert the draws in backsamp.bigS (1,2,3,...)
    # to matrix of zeros and ones.
    # Includes initial state, we remove in the below ([-1])
    matSS <- matrix(as.integer(matrix(seq(1,h), bigt, h,
                                      byrow=TRUE)==bs$backsamp.bigS[-1]),
                    bigt, h)

    # Output:
    # 1) SS = state-space
    # 2) transitions = transition matrix
    # 3) llf = log-likelihood computed by the filter
    ss <- list(SS=matSS, transitions=bs$transmat, llf=ff$llh)
    return(ss)
}




####################################################################
# Sample the error covariances for each state
####################################################################

Sigma.draw <- function(m, h, ss, hreg, Sigmai)
{
    # Draw the variances for the MSVAR from an inverse Wishart --
    # draws h matrices which are returned as an m^2 x h matrix.  These
    # then need to be unwound into the matrices we need.

    # Unscaled inverse Wishart draws
    df <- hreg$df
    tmp <- array(unlist(lapply(sapply(1:h,
                                      function(i) {rwishart(N=1,
                                                            df=df[i]+m+1,
                                                            Sigma=diag(m))},
                                      simplify=FALSE), solve)), c(m,m,h))

    wishscale <- hreg$Sigmak

    # Scale the draws for each regime
    for (i in 1:h)
    {
#        dc <- tcholpivot(hreg$Sigmak[,,i])
        dc <- t(chol(hreg$Sigmak[,,i]))
        Sigmai[,,i] <- dc%*%(tmp[,,i]*(df[i]+m+1))%*%t(dc)
        wishscale[,,i] <- wishscale[,,i]*(df[i]+m+1)
    }


    return(list(Sigmai=Sigmai, wishscale=wishscale))
}

####################################################################
# Sample the regressors for each state
####################################################################

Beta.draw <- function(m, p, h, Sigmai, hreg, init.model, Betai)
{
    mmp1 <- m*(m*p+1)
    Beta.cpprec <- Beta.cpv <- array(NA, c(mmp1,mmp1,h))

    for(i in 1:h)
    {
        # Set up some the moments we need
        # 1) (X'X)^{-1} with the prior included -- see regression step
        # 2) Use the fact that the inverses and Kroneckers can be
        # switched around.
        # 3) Compute the posterior precision and variance on the fly
        # for later use

        # Robust version that uses stable choleskys to do the
        # computation
        Beta.cpprec[,,i] <- kronecker(Sigmai[,,i], hreg$moment[[i]]$Sxx
                                 + init.model$H0)
        Beta.cpv[,,i] <- chol2inv(chol(Beta.cpprec[,,i]))

        bcoefs.se <- t(chol(Beta.cpv[,,i]))

        # Sample the regressors
        Betai[,,i] <- hreg$Bk[,,i] +
            matrix(bcoefs.se%*%matrix(rnorm(m^2*p+m), ncol=1), ncol=m)
    }

    # Returns
    # 1) Betai = the MCMC draw
    # 2) Beta.cpm = conditional posterior mean
    # 3) Beta.cpprec = conditional posterior precision
    # 4) Beta.cpv = conditional posterior variance

    return(list(Betai=Betai, Beta.cpm=hreg$Bk, Beta.cpprec=Beta.cpprec,
                Beta.cpv=Beta.cpv))
}

####################################################################
# Update the residuals for each state
# Start with m+2 residual, since the others are part of the dummy
# observations for the prior.
####################################################################
residual.update <- function(m, h, init.model, Betai, e)
{
    for(i in 1:h)
    {
        e[,,i] <- init.model$Y[(m+2):nrow(init.model$Y),] -
            init.model$X[(m+2):nrow(init.model$X),]%*%Betai[,,i]
    }
    return(e)
}

######################################################################
# Function to do the random permuation / state labeling steps for the
# models.  This function handles both the random permutation steps and
# the identified models based on either tbe Beta.idx or the Sigma.idx
# input args.
######################################################################
PermuteFlip <- function(x, h, permute, Beta.idx, Sigma.idx)
{
    # Random permutation step
    if(permute==TRUE){
        rp <- sort(runif(h), index.return=TRUE)$ix

        Sigmaitmp <- x$Sigmai; Betaitmp <- x$Betai
        etmp <- x$e; sstmp <- x$ss; df <- x$df

        for(i in 1:h)
        {
            Sigmaitmp[,,i] <- x$Sigmai[,,rp[i]]
            Betaitmp[,,i] <- x$Betai[,,rp[i]]
            etmp[,,i] <- x$e[,,rp[i]]
            sstmp$SS[,i] <- x$ss$SS[,rp[i]]
        }

        Q <- x$Q[rp,rp]
        sstmp$transitions <- x$ss$transitions[rp, rp]
        df <- df[rp]
        return(list(Betai=Betaitmp, Sigmai=Sigmaitmp,
                    ss=sstmp, e=etmp, Q=Q, df=df))
    }

    # Sort regimes based on values of regressors in Beta using the
    # Beta.idx indices

    if(permute==FALSE & is.null(Beta.idx)==FALSE){
        # Find the sort indices for the objects
        ix <- sort(x$Betai[Beta.idx[1], Beta.idx[2],], index.return=TRUE)$ix

        Sigmaitmp <- x$Sigmai
        Betaitmp <- x$Betai
        etmp <- x$e
        sstmp <- x$ss
        df <- x$df

        for(i in 1:h)
        {
            Sigmaitmp[,,i] <- x$Sigmai[,,ix[i]]
            Betaitmp[,,i] <- x$Betai[,,ix[i]]
            etmp[,,i] <- x$e[,,ix[i]]
            sstmp$SS[,i] <- x$ss$SS[,ix[i]]
        }

        Q <- x$Q[ix,ix]
        sstmp$transitions <- x$ss$transitions[ix, ix]
        df <- df[ix]
        return(list(Betai=Betaitmp, Sigmai=Sigmaitmp, ss=sstmp,
                    e=etmp, Q=Q, df=df))
    }

    # Sort regimes based on values of the Sigma matrix element
    # selected.

    if(permute==FALSE & is.null(Sigma.idx)==FALSE){
        # Find the sort indices for the objects
        ix <- sort(x$Sigmai[Sigma.idx, Sigma.idx,], index.return=TRUE)$ix

        Sigmaitmp <- x$Sigmai
        Betaitmp <- x$Betai
        etmp <- x$e
        sstmp <- x$ss
        df <- x$df

        for(i in 1:h)
        {
            Sigmaitmp[,,i] <- x$Sigmai[,,ix[i]]
            Betaitmp[,,i] <- x$Betai[,,ix[i]]
            etmp[,,i] <- x$e[,,ix[i]]
            sstmp$SS[,i] <- x$ss$SS[,ix[i]]
        }

        Q <- x$Q[ix,ix]
        sstmp$transitions <- x$ss$transitions[ix, ix]
        df <- df[ix]
        return(list(Betai=Betaitmp, Sigmai=Sigmaitmp, ss=sstmp,
                    e=etmp, Q=Q, df=df))
    }
}

######################################################################
# marginal.moments() Computes the marginal moments for the conditional
# posterior of Beta and Sigma for the given prior.  These computations
# are necessary for later evaluation of the marginal likelihood.  So
# this function is called conditionally after the regression draws are
# made in order to get the conditional moments.
######################################################################

## marginal.moments <- function(h, Betai, Sigmai, hreg, init.model)
## {
##     # Constants we need up front
##     tmp <- dim(Betai)
##     mp1 <- tmp[1]
##     m <- tmp[2]
##     h <- tmp[3]

##     prior.b0m <- as.vector(rbind(diag(m), matrix(0, m*(p-1)+1, m)))

##     # Storage

##     for(i in 1:h)
##     {
##         Sinv <- chol2inv(chol(Sigmai[,,i]))
##         SiginvXX <- kronecker(Sinv, hreg$moment[[i]]$Sxx)

##         # Compute the conditional posterior precision, B|S^{-1}
##         beta.prec <- kronecker(Sinv, init.model$H0) + SiginvXX

##         # Compute conditional posterior variance, B|S
##         beta.cpv <- chol2inv(chol(beta.prec))

##         # Compute the conditional posterior mean
##         bhat <- kronecker(Sinv, diag(mp1)) %*% as.vector(hreg$moment[[i]]$Sxy)
##         beta.cpm <- beta.cpv %*% kronecker(Sinv, init.model$H0) + bhat

##         # Flatten things out for later storage

##     }
##     return()
## }

####################################################################
# Full MCMC sampler for MSBVAR models with SZ prior
####################################################################

gibbs.msbvar <- function(x, N1=1000, N2=1000,
                         permute=TRUE,
                         Beta.idx=NULL, Sigma.idx=NULL,
                         Q.method="MH")
{
    # Do some sanity / error checking on the inputs so people know
    # what they are getting back

    if(permute==TRUE & is.null(Beta.idx)==FALSE){
        cat("You have set permute=TRUE which means you requested a random permutation sample.\n Beta.idx argument will be ignored")
        }

    if(permute==TRUE & is.null(Sigma.idx)==FALSE){
        cat("You have set permute=TRUE which means you requested a random permutation sample.\n Sigma.idx argument will be ignored")
        }
    if(permute==FALSE & (is.null(Beta.idx)==TRUE &
        is.null(Sigma.idx)==TRUE)){
        cat("You have set permute=FALSE, but failed to provide an identification to label the regimes.\n Set Beta.idx or Sigma.idx != NULL")
    }

    # Set the sampler method for Q based on inputs
    if(Q.method=="Gibbs") Qsampler <- Q.drawGibbs
    if(Q.method=="MH") Qsampler <- Q.drawMH

    # Initializations
    init.model <- x$init.model
    hreg <- x$hreg
    Q <- x$Q
    fp <- x$fp
    m <- x$m
    p <- x$p
    h <- x$h
    alpha.prior <- x$alpha.prior
    TT <- nrow(fp)

    # Initialize the Gibbs sampler parameters
    e <- hreg$e
    Sigmai <- hreg$Sigma
    Betai <- hreg$Bk

    # Burnin loop
    for(j in 1:N1)
    {
        # Smooth / draw the state-space, with checks for the
        # degeneracy of the state-space on the R side.
##         oldtran <- matrix(0,h,h)

##         while(sum(diag(oldtran)==0)>0)
##         {
##             ss <- SS.ffbs(e, (TT+p), m, p, h, Sigmai, Q) # may need p=0
##             oldtran <- ss$transitions
## #            print(oldtran)
##         }

        ss <- SS.ffbs(e, TT, m, p, h, Sigmai, Q)

#        print(ss$transitions)
        # Draw Q
        Q <- Qsampler(Q, ss$transitions, prior=alpha.prior, h)

        # Update the regression step
        hreg <- hregime.reg2(h, m, p, ss$SS, init.model)

        # Draw the variance matrices for each regime
        Sout <- Sigma.draw(m, h, ss, hreg, Sigmai)
        Sigmai <- Sout$Sigmai

        # Sample the regression coefficients
        Bout <- Beta.draw(m, p, h, Sigmai, hreg, init.model, Betai)
        Betai <- Bout$Betai

        # Update the residuals
        e <- residual.update(m, h, init.model, Betai, e)

        # Now do the random permutation of the labels with an MH step
        # for the permutation or sort the regimes based on the
        # Beta.idx or Sigma.idx conditions.
        pf.obj <- PermuteFlip(x=list(Betai=Betai, Sigmai=Sigmai,
                                     ss=ss, e=e, Q=Q, df=hreg$df),
                              h, permute, Beta.idx, Sigma.idx)

        # Replace sampler objects with permuted / flipped ones
        Betai <- pf.obj$Betai
        Sigmai <- pf.obj$Sigmai
        ss <- pf.obj$ss
        e <- pf.obj$e
        Q <- pf.obj$Q

        # Print out some iteration information
        if(j%%1000==0) cat("Burn-in iteration : ", j, "\n")

    }

    # End of burnin loop!

    gc(); gc()

    # Declare the storage for the return objects
    ss.storage <- vector("list", N2)
    transition.storage <- array(NA, c(h,h,N2))
    Beta.storage <- matrix(NA, N2, (m^2*p+m)*h)
    Sigma.storage <- matrix(NA, N2, m*(m+1)*0.5*h)
    Q.storage <- matrix(NA, N2, h^2)
    llf <- matrix(NA, N2, 1)

    # Storage for the conditional posterior moments when permute==TRUE
    if(permute==TRUE)
    {
        # Conditional posterior moments
        Beta.cpm <- matrix(NA, N2, (m^2*p+m)*h)
        Sigma.wishscale <- matrix(NA, N2, 0.5*m*(m+1)*h)
        tmpdim <- m*((m*p)+1)
        Beta.cpprec <- Beta.cpv <- matrix(NA, N2, 0.5*tmpdim*(tmpdim+1)*h)
        # Sigma / Wishart df moments
        df.storage <- matrix(NA, N2, h)
        # Q, conditional posterior
        Q.cp <- matrix(NA, N2, h^2)
    }

    # Main loop after the burnin -- these are the sweeps that are kept
    # for the final posterior sample.
    # Note, we need to do the storage AFTER the permutation steps so
    # we get the labeling correct.
    for(j in 1:N2)
    {
        # Draw the state-space
        # Smooth / draw the state-space, with checks for the
        # degeneracy of the state-space on the R side.
        ## oldtran <- matrix(0,h,h)

        ## while(sum(diag(oldtran)==0)>0)
        ## {
        ##     ss <- SS.ffbs(e, (TT+p), m, p, h, Sigmai, Q) # may need p=0
        ##     oldtran <- ss$transitions
        ## }

        ss <- SS.ffbs(e, TT, m, p, h, Sigmai, Q)

        # Draw Q
        Q <- Qsampler(Q, ss$transitions, prior=alpha.prior, h)

        # Update the regression step
        hreg <- hregime.reg2(h, m, p, ss$SS, init.model)

        # Draw the variance matrices for each regime
        Sout <- Sigma.draw(m, h, ss, hreg, Sigmai)
        Sigmai <- Sout$Sigmai

        # Sample the regression coefficients
        Bout <- Beta.draw(m, p, h, Sigmai, hreg, init.model, Betai)
        Betai <- Bout$Betai

        # Update the residuals
        e <- residual.update(m, h, init.model, Betai, e)

        # Now do the random permutation of the labels with an MH step
        # for the permutation or sort the regimes based on the
        # Beta.idx or Sigma.idx conditions.
        pf.obj <- PermuteFlip(x=list(Betai=Betai, Sigmai=Sigmai,
                              ss=ss, e=e, Q=Q, df=hreg$df),
                              h, permute, Beta.idx, Sigma.idx)

        # Replace sampler objects with permuted / flipped ones
        Betai <- pf.obj$Betai
        Sigmai <- pf.obj$Sigmai
        ss <- pf.obj$ss
        e <- pf.obj$e
        Q <- pf.obj$Q
        df <- pf.obj$df

        # Store things -- after flip / identification
        Sigma.storage[j,] <- as.vector(apply(Sigmai, 3, vech))
        Beta.storage[j,] <- as.vector(Betai)
        ss.storage[[j]] <- as.bit.integer(as.integer(ss$SS[,1:(h-1)]))
        transition.storage[,,j] <- ss$transitions
        Q.storage[j,] <- as.vector(Q)

        # Store the llf
        llf[j,1] <- ss$llf

        # Store conditional posterior moments
        if(permute==TRUE)
        {
            Beta.cpm[j,] <- as.vector(Bout$Beta.cpm)
            Beta.cpprec[j,] <- as.vector(apply(Bout$Beta.cpprec, 3, vech))
            Beta.cpv [j,] <- as.vector(apply(Bout$Beta.cpv, 3, vech))
            Sigma.wishscale[j,] <- as.vector(apply(Sout$wishscale, 3, vech))
            df.storage[j,] <- df
            Q.cp[j,] <- as.vector(ss$transitions + alpha.prior)
        }

        # Print out some iteration information
        if(j%%1000==0) cat("Final iteration : ", j, "\n")
    }

    gc(); gc()

    # Now make an output object and provide classing
    class(ss.storage) <- c("SS")

    if(permute==FALSE)
    {
        output <- list(Beta.sample=mcmc(Beta.storage),
                       Sigma.sample=mcmc(Sigma.storage),
                       Q.sample=mcmc(Q.storage),
                       transition.sample=transition.storage,
                       ss.sample=ss.storage,
                       llf=llf,
                       init.model=init.model,
                       alpha.prior=alpha.prior,
                       h=h,
                       p=p,
                       m=m)
    } else {
        output <- list(Beta.sample=mcmc(Beta.storage),
                       Sigma.sample=mcmc(Sigma.storage),
                       Q.sample=mcmc(Q.storage),
                       transition.sample=transition.storage,
                       ss.sample=ss.storage,
                       llf=llf,
                       init.model=init.model,
                       alpha.prior=alpha.prior,
                       h=h,
                       p=p,
                       m=m,
                       Beta.cpm=Beta.cpm,
                       Beta.cpprec=Beta.cpprec,
                       Beta.cpv=Beta.cpv,
                       Sigma.wishscale=Sigma.wishscale,
                       df=df.storage,
                       Q.cp=Q.cp)
    }

    # Class the equation name information
    class(output) <- c("MSBVAR")
    attr(output, "eqnames") <- attr(init.model, "eqnames")

    # Keep track of the regime identification properties from
    # estimation.
    attr(output, "permute") <- permute
    attr(output, "Beta.idx") <- Beta.idx
    attr(output, "Sigma.idx") <- Sigma.idx
    attr(output, "Qsampler") <- Q.method

    # ts attributes for inputs to later plotting and summary functions
    # for the states.
    tmp <- tsp(init.model$y)
    attr(output, "start") <- tmp[1]
    attr(output, "end") <- tmp[2]
    attr(output, "freq") <- tmp[3]

    return(output)
}

################################################################
# Deprecated code
################################################################

# R code for legacy purposes
#SS.ffbs(e, TT+p, m, p, h, Sigmai, Q)
#SS.ffbs <- function(e, bigt, m, p, h, sig2, Q)
#{
#  TT <- bigt-p
#
#  # Forward-filtering
#  fHam <- .Fortran("HamiltonFilter",
#                   bigt=as.integer(bigt),
#                   m = as.integer(m), p = as.integer(p), h = as.integer(h),
#                   e = e,
#                   sig2 = sig2,
#                   Q = Q,
#                   f = double(1),
#                   filtprSt = matrix(0,as.integer(TT),as.integer(h))
#                   )
#  fpH <- fHam$filtprSt
#
#  # Backwards-sampling
#  SS <- matrix(0, nrow=TT, ncol=h)
#  SS[TT,] <- bingen(fpH[TT,], Q, 1, h)
#  for (t in (TT-1):1)
#  {
#    SS[t,] <- bingen(fpH[t,], Q, which(SS[t+1,]==1), h)
#  }
#
#  # Construct STT (TTx1 vector of regimes)
#  STT <- rep(0,TT)
#  for (ist in 1:TT) {
#    STT[ist] <- which(SS[ist,]==1)
#  }
#
#  # Construct transition matrix
#  transmat <- matrix(0,h,h)
#
#  for (ist in 1:(TT-1)) {
#    transmat[STT[ist+1], STT[ist]] <- transmat[STT[ist+1], STT[ist]] + 1
#  }
#
#  ss <- list(SS=SS,
#             transitions=transmat)
#
#  return(ss)
#}
#
#bingen <- function(prob, Q, st1, h)
#{
#  i <- 1
#  while(i<h)
#  {
#    pr0 <- prob[i]*Q[st1,i]/sum(prob[i:h]*Q[st1,i:h])
#    if(runif(1)<=pr0)
#    { return(diag(h)[i,]) } else { st1 <- i <- i+1 }
#  }
#  return(diag(h)[h,])
#}
