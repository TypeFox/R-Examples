# posterior.fit.MSBVAR
# Bridge sampler for MSBVAR marginal likelihoods.
# 2012-05-21 : Adapted from WRD code

# Initial functions for density evaluations

# Evaluate log density of multivariate normal
logdmvnorm <- function(x, mean, siginv, n, logdet){
    ymu <- x - mean
    distval <- crossprod(ymu,siginv) %*% ymu
    logretval <- -.5*(n*log(2*pi) + logdet + distval)
    return(logretval)
}

# Evaluate log density of the inverse Wishart
logdinvwish <- function(Y, alpha, S)
{
    d <- nrow(S)
    lgammapart <- 0
    for (j in 1:d) {
      lgammapart <- lgammapart + lgamma((2*alpha + 1 - j)/2)
    }
    denom <- lgammapart + (d * (d - 1)/4)*log(pi)

    logdetS <- determinant(S)$modulus[1]
    invY <- chol2inv(chol(Y))
    logdetinvY <- determinant(invY)$modulus[1]
    trSinvY <- sum(diag(S %*% invY))
    num <- alpha*logdetS + (alpha+(d+1)/2)*logdetinvY + -trSinvY
    return(num-denom)
}

# Log density of the MVN regressors under the importance density.
qSlogdmvnorm <- function(nbeta,nvar,bigK,S,betadraw,
                         S.beta.cpm,S.beta.cpinv,S.beta.logdet) {

    qSlogdmvnorm <- .Fortran('qSlogdmvnorm',
                             nbeta=as.integer(nbeta),
                             nvar=as.integer(nvar),
                             bigK=as.integer(bigK),
                             S=as.integer(S),
                             betadraw,
                             S.beta.cpm,
                             S.beta.cpinv,
                             S.beta.logdet,
                             qSmvn=matrix(0,S,1)) #double(S))

    return(qSlogdmvnorm = qSlogdmvnorm$qSmvn)
}

# Log density of iW for variances under the importance density
qSlogdinvwish <- function(nvar,bigK,S,sigdraw,
                          S.bigSig.cp.wishdf,S.bigSig.cp.wishscale) {

    qSlogdinvwish <- .Fortran('qSlogdinvwish',
                              nvar=as.integer(nvar),
                              bigK=as.integer(bigK),
                              S=as.integer(S),
                              sigdraw,
                              S.bigSig.cp.wishdf,
                              S.bigSig.cp.wishscale,
                              qSliw=double(S))

    return(qSlogdinvwish = qSlogdinvwish$qSliw)
}

# Forward filter with internal regression steps for regimes.

ForwardFilterReg <- function(nvar, bigY, bigX, bigK, bigT,
                             nbeta, betadraw, sigdraw, xidraw)
{

  ff <- .Fortran('ForwardFilterReg',
                 nvar=as.integer(nvar),
                 bigY,bigX,
                 bigK=as.integer(bigK), bigT=as.integer(bigT),
                 nbeta=as.integer(nbeta),
                 betadraw, sigdraw, xidraw,
                 llh=double(1),
                 pfilt=matrix(0,bigT+1,bigK))

  return(list(llh=ff$llh, pfilt=ff$pfilt))
}

# Random Wishart generator
"rwish" <-
  function(v, S) {
    if (!is.matrix(S))
      S <- matrix(S)
    if (nrow(S) != ncol(S)) {
      stop(message="S not square in rwish().\n")
    }
    if (v < nrow(S)) {
      stop(message="v is less than the dimension of S in rwish().\n")
    }
    p <- nrow(S)
    CC <- chol(S)
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, v:(v-p+1)))
    if(p > 1) {
      pseq <- 1:(p-1)
      Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
    }
    return(crossprod(Z %*% CC))
  }


# riwish generates a draw from the inverse Wishart distribution
# (using the above Wishart generator)

"riwish" <-
  function(v, S) {
    return(solve(rwish(v,solve(S))))
  }


###########################################################
# posterior.fit.MSBVAR()
#
# x = permuted gibbs.msbvar() object
###########################################################

posterior.fit.MSBVAR <- function(x, maxiterbs=500)
{
    #########################################################
    # Logical checks for user inputs
    #########################################################
    if(attr(x, "permute")==FALSE) {
        stop("Must use random permutation with gibbs.msbvar() before creating inputs to this function.  Set permute=TRUE in you earlier call to gibbs.msbvar()\n")
    }

    # Set up constants / preset data that we need for the
    # computations.

    im <- x$init.model
    h <- x$h; p <- x$p; m <- x$m
    N <- nrow(x$llf)
    Q.prior <- x$alpha.prior
    wdof <- x$df
    mp1 <- m*p+1

    # Regression-step data
    bigY <- im$Y
    bigX <- im$X
    bigT <- im$pfit$capT

    # Evaluate log prior density first

    # Static log prior steps -- these are the conjugate prior moments
    # with SZ prior.
    prior.B0v <- solve(im$H0)
    prior.b0m <- as.vector(rbind(diag(m), matrix(0, m*(p-1)+1, m)))
    prior.B0prec <- kronecker(diag(m), im$H0)

    # Storage
    log.prior <- matrix(NA, N, 1)

    for(i in 1:N)
    {
        #########################################################
        # Evaluate the log prior for the posterior sample
        #########################################################

        # Dirichlet for the transition matrix, SFS, p. 432, A5
        lpQ <- sum(lgamma(rowSums(Q.prior)) +
                   rowSums((Q.prior - 1) * log(matrix(x$Q.sample[i,], h, h))
                           - lgamma(Q.prior)))

        lpB <- lpS <- 0
        Sigtmp <- matrix(x$Sigma.sample[i,], ncol=h)
        Btmp <- matrix(x$Beta.sample[i,], ncol=h)

        for(j in 1:h)
        {
            # Beta / regression step
            logdetval <- log(det(kronecker(xpnd(Sigtmp[,j]), prior.B0v)))
            lpB <- lpB + logdmvnorm(Btmp[,j], prior.b0m,
                                    prior.B0prec, mp1, logdetval)


            lpS <- lpS + logdinvwish(xpnd(Sigtmp[,j]), m+1, im$S0)
        }

        # Now add up the elements and store them

        log.prior[i,1] <- lpQ + lpB + lpS

        #########################################################
        # End evaluate log prior for the posterior sample
        #########################################################

    } # Ends loop over the MCMC sample for prior.


    ###################################################################
    # Construct importance density q for the marglik
    # SFS 2006, page 145, eq. (5.36)
    #
    # For K=1 and sig2 fixed, the importance density simplifies to
    # sums of normals since joint posterior is normal.
    #
    # The importance density q() is constructed from S
    # randomly selected MCMC draws of the allocation
    # vector boldS.
    #
    # First, we determine what S, the number of importance samples
    # should be.  Since q() needs to only be a rough approximation to
    # the complete-data posterior density, typically S is much smaller
    # than the total number of MCMC draws. i.e., S << N2 = N
    #
    # The following follows page 51 of SFS's Matlab companion manual.
    #
    # SFS recommends choosing S by the following formula:
    #    S = min(M0*K!, Smax, M)
    #
    # - M0 is the expected number of times that the construction of
    #   the importance density q will be based on a particular mode
    #   of the mixture posterior.
    # - K is the number of states/regimes.
    # - Smax is an upper limit for S.
    # - M is the number of MCMC draws (num.mcmc)
    ####################################################################

    Smax <- 2000

    if (h <= 3) {
        M0 <- 100
    } else if (h==4) {
        M0 <- 25
    } else {
        M0 <- 5
    }

    # Compute K! (K factorial) using n! = Gamma(n+1)
    # use lgamma() function to avoid potential underflow
    S <- round(M0*exp(lgamma(h+1)))

    # Finalize value of S
    S <- min(S, Smax, N)

    #########################################################
    # Now, we extract the posterior moments corresponding
    # to S components
    #
    # Obtain a random permutation of the integers from
    # 1 to N.

    idx <- sample(N, N)

    # Extract S randomly chosen posterior moments
    # indexchoose is S by 1 vector of randomly selected
    # integers from 1-N, without replacement and
    # sorted from smallest to largest
    # e.g., indexchoose = 14,15,...,4867,4933, etc.

    indexchoose <- sort(idx[(N-S+1):N])

    ################################################################
    # Obtain a subset of the conditional posterior moments from the
    # MCMC draws. Use randomly permuted output (not identified
    # output).  Note that this is checked at the top of the function.
    ################################################################
    # Subsets of the posterior from the permuted MCMC output
    ################################################################
    # Transition matrix
    S.xi.post  <- array(x$Q.cp[indexchoose,], c(S, h, h))

    # Transition matrix draws
    xi.draws <- array(x$Q.sample, c(N, h, h))

    # Conditional posterior mean of regressors, pre-draw
    S.beta.cpm <- array(x$Beta.cpm[indexchoose,], c(S,mp1*m,h))


    # Conditional posterior variances of parameters.  Note that we are
    # putting these into arrays for faster access and indexing.
    # Before we do this we need to rebuild the matrices since vech()
    # was used to simplify and reduce memory storage

    wishdf <- x$df[indexchoose,]  # S.bigSig.cp.wishdf

    tmp.beta.cpv <- x$Beta.cpv[indexchoose,]

    tmp.beta.cpprec <- x$Beta.cpprec[indexchoose,]

    tmp.wishscale <- x$Sigma.wishscale[indexchoose,]

    S.bigSig.cp.wishscale <- array(NA, c(S, m, m, h))
    for(i in 1:S)
    {
        tmp <- apply(matrix(tmp.wishscale[i,], ncol=h), 2, xpnd)
        tmp1 <- array(tmp, c(m,m,h))

        for(j in 1:h)
        {
            S.bigSig.cp.wishscale[i,,,j] <- tmp1[,,j]
        }
    }

    # Translations from WRD to PTB objects
    #
    # Note, new objects eliminate any unnecessary array dimensions.
    #
    # S.xi.post = same
    # S.beta.cpm = same
    # S.beta.cpv = c(S,nbeta*nvar,nbeta*nvar,bigK)
    # S.bigSig.cp.wishdf = c(S,1,bigK)
    # S.bigSig.cp.wishscale = c(S,nvar,nvar,bigK)

    # Now rebuild the relevant matrices and compute necessary
    # functions of all of these matrices.  Some of this has been moved
    # here from lower sections of the earlier code so we can do more
    # in the same set of loops.

    ###############################################################
    # Define interim storage and computations
    # Pre-compute information matrix (inverse of variance -- same as
    # the precision!), Cholesky decomposition and log determinant for
    # speed
    #
    # Conditional posterior variance and precision
    # Note: S.beta.cpinv = S.beta.cpprec
    #

    S.beta.cpv <- S.beta.cpprec <- S.beta.vchol <- array(NA, c(S, m*mp1, m*mp1, h))
    S.beta.logdet <- array(NA, c(S, 1, h))

    ################################################################


    ################################################################
    # Loop to do all of the computations outlined above...
    ################################################################

    for(iters in 1:S)
    {
        tmp1 <- matrix(tmp.beta.cpv[iters,], ncol=h)
        tmp2 <- matrix(tmp.beta.cpprec[iters,],ncol=h)
            ####################################################
            # Populate the moments for the conditional posterior
            # variance and precision.
            ####################################################
        S.beta.cpv[iters,,,] <- array(apply(tmp1, 2, xpnd),
                                      c(m*mp1, m*mp1, h))
        S.beta.cpprec[iters,,,] <- array(apply(tmp2, 2, xpnd),
                                         c(m*mp1, m*mp1, h))
        S.beta.vchol[iters,,,] <- array(apply(S.beta.cpv[iters,,,], 3,
                                              chol, pivot=TRUE),
                                        c(m*mp1, m*mp1, h))
        S.beta.logdet[iters,1,] <- matrix(unlist(apply(S.beta.cpv[iters,,,], 3,
                                               determinant,logarithm=TRUE)),
                                        ncol=h)[1,]
    }


    ################################################################
    # Initialize the marginal likelihood storage
    ################################################################

    # IS = Importance Sampling
    q.IS.log.prior <- rep(0,N)
    q.IS.log.lik   <- rep(0,N)
    q.IS.log.q     <- rep(0,N)

    # RI = Reciprocal Importance Sampling
    # n.b., evaluated at MCMC draws (SFS 2006 pg 148)
    q.RI.log.prior <- log.prior
    q.RI.log.lik   <- rep(0,N)
    q.RI.log.q     <- rep(0,N)

    # n.b., cannot just use log likelihood values from Gibbs
    # since pmax might influence results

    ################################################################
    # Choose the components for sampling from the mixture importance
    # density sample N random integers between 1 and S, with
    # replacement.
    #
    # Here qs is (N by 1) col vector of random integers between 1 and
    # S
    ################################################################
    qs <- sample(S, N, replace=TRUE)

    # Storage for the draws of regressors
    q.beta.draw <- array(0, c(m*mp1, 1, h))

    cat('\n Begin Construction of Importance Density... \n')

    iterq <- 1

    # Begin loop for importance density construction
    for (iterq in 1:N) {

    # qS.RI : q(thetamod^(m)_K) in equation (5.43) page 148, reciprocal IS
    # qs.IS : q(thetamod^(l)_K) in equation (5.39) page 146, importance sampling
        qS.IS <- rep(0,S)
        qS.RI <- rep(0,S)

        if (h > 1) {
    # Draw xi from q for Importance Sampling
    # Gibbs draw for xi from Dirichlet(S.xi.post) using Gamma
    # distribution and then normalizing, SFS 2006 pg 433
            tmpxi <- rgamma(h*h, S.xi.post[qs[iterq],,], 1)
            q.xi.draw <- matrix(tmpxi, h, h)
            q.xi.draw <- q.xi.draw / rowSums(q.xi.draw)


    # Evaluate the above draw using transition matrix part
    # of importance density q for Importance Sampling
    # see page 432, A.5

            tmpeta <- array(NA, c(S, h, h))
            for (iterxi1 in 1:h) {
                for (iterxi2 in 1:h) {
                    tmpeta[,iterxi2,iterxi1] <-
                        q.xi.draw[iterxi2,iterxi1]
                }}

    # looks like there are some rounding difference between this and Matlab
            qS.IS <- qS.IS + rowSums(lgamma(apply(S.xi.post, 2, rowSums)) +
                                     apply((S.xi.post - 1) *
                                           log(tmpeta) -
                                           lgamma(S.xi.post), 2, rowSums))


    # Evaluate MCMC xi draws using using transition matrix part
    # of importance density q for Importance Sampling
    # see page 432, A.5
            tmpeta <- array(NA, c(S, h, h))
            for (iterxi1 in 1:h) {
                for (iterxi2 in 1:h) { tmpeta[,iterxi2,iterxi1] <-
                                           xi.draws[iterq,iterxi2,iterxi1]
                                   }}

            qS.RI <- qS.RI + rowSums(lgamma(apply(S.xi.post, 2, rowSums)) +
                                     apply((S.xi.post - 1) *
                                           log(tmpeta) - lgamma(S.xi.post),
                                           2, rowSums))

        } else {
            q.xi.draw <- 1
        }


  # Draw beta from q for Importance Sampling
  # simulate Beta from its multivariate normal conditional
        for (iterk in 1:h) {
            q.beta.draw[,,iterk] <- S.beta.cpm[qs[iterq],,iterk] +
                crossprod(S.beta.vchol[qs[iterq],,,iterk], rnorm(m*mp1))
        }

        qS.IS <- qS.IS + qSlogdmvnorm(mp1, m, h, S, q.beta.draw, S.beta.cpm,
                                      S.beta.cpprec, S.beta.logdet)

        beta.draws <- matrix(x$Beta.sample[iterq,], ncol=h)

        qS.RI <- qS.RI + qSlogdmvnorm(mp1, m, h, S, beta.draws,
                                      S.beta.cpm, S.beta.cpprec,
                                      S.beta.logdet)


  # For h=1, sig2 fixed, the qs.IS and qs.RI will all have the same values
  # in each entry, but when h>1, these will be different.

  # Draw bigSig for Importance Sampling
        q.bigSig.draw <- array(NA, c(m,m,h))
        stmp <- matrix(tmp.wishscale[qs[iterq],], ncol=h)
        stmp <- array(apply(stmp, 2, xpnd), c(m,m,h))
        dfs <- wishdf[qs[iterq],]
        for (iterk in 1:h) {
            q.bigSig.draw[,,iterk] <- riwish(dfs[iterk]+m+1,
                                             stmp[,,iterk])
        }

        qS.IS <- qS.IS + qSlogdinvwish(m, h, S, q.bigSig.draw,
                                       array(wishdf, c(S,1,h)),
                                       S.bigSig.cp.wishscale)

        Sigma.draws <- array(apply(matrix(x$Sigma.sample[iterq,], ncol=h),
                                   2, xpnd), c(m,m,h))

        qS.RI <- qS.RI + qSlogdinvwish(m, h, S, Sigma.draws,
                                       array(wishdf, c(S,1,h)), S.bigSig.cp.wishscale)

        # Average over qS, i.e., eq (5.36) on page 145
        max.qS.IS <- max(qS.IS)
        max.qS.RI <- max(qS.RI)
        q.IS.log.q[iterq] <- max.qS.IS + log(mean(exp(qS.IS - max.qS.IS)))
        q.RI.log.q[iterq] <- max.qS.RI + log(mean(exp(qS.RI - max.qS.RI)))

        # Importance Sampling

        # Evaluate likelihood for qsample at iteration iterq

        ff <- ForwardFilterReg(m, bigY, bigX, h, bigT, mp1,
                               q.beta.draw, q.bigSig.draw, q.xi.draw)

        # Store the likelihood
        q.IS.log.lik[iterq] <- ff$llh

        # Evaluate q draw(s) at MVN|iW prior density, summing over
        # regimes k

        q.log.beta.prior <- 0
        q.log.bigSig.prior <- 0
        for (iterk in 1:h) {
            logdetval <- log(det(kronecker(q.bigSig.draw[,,iterk],
                                           prior.B0v)))

            q.log.beta.prior <- q.log.beta.prior +
                logdmvnorm(q.beta.draw[,,iterk], prior.b0m,
                           prior.B0prec,
                           mp1, logdetval)

            q.log.bigSig.prior <- q.log.bigSig.prior +
                logdinvwish(q.bigSig.draw[,,iterk],
                            m+1,
                            im$S0)
        }

  # If h > 1, then evaluate prior on transition matrix xi
  # see page 432, A.5
        q.log.xi.prior <- 0
        if (h > 1) {
            q.log.xi.prior <- sum(lgamma(rowSums(x$alpha.prior)) +
                                  rowSums((x$alpha.prior - 1) * log(q.xi.draw) -
                                          lgamma(x$alpha.prior)))
        }

  # Sum up the above
        q.IS.log.prior[iterq] <- q.log.beta.prior + q.log.bigSig.prior + q.log.xi.prior


  # Reciprocal IS
  # only have to evaluate likelihood, since can get prior from Gibbs
  # this likelihood is necessary if p does not equal pmax

  # Evaluate likelihood for qsample
        ff <- ForwardFilterReg(m, bigY, bigX, h, bigT, mp1,
                               beta.draws, Sigma.draws,
                               xi.draws[iterq,,])

        q.RI.log.lik[iterq] <- ff$llh

  # write output to screen
        if (iterq %% 1000 == 0) cat('   Importance Density: ', iterq, '\n',sep='')

    } # End importance sampling loop

    cat(' End Construction of Importance Density \n')

# Define values that go into the ratio calculations for IS and RI

# IS = Importance Sampling, all (N by 1)
    priorIS  <- q.IS.log.prior
    loglikIS <- q.IS.log.lik
    qIS      <- q.IS.log.q

# RI = Reciprocal Importance Sampling, all (N by 1)
    priorRI  <- q.RI.log.prior
    loglikRI <- q.RI.log.lik
    qRI      <- q.RI.log.q

# Construct IS and RI estimators (eqs 5.39 and 5.45, resp)
    ratio.IS <- q.IS.log.lik + q.IS.log.prior - q.IS.log.q
    maxratio.IS <- max(ratio.IS)
    marglik.IS <- maxratio.IS + log(mean(exp(ratio.IS-maxratio.IS)))

    ratio.RI <- q.RI.log.q - q.RI.log.lik - q.RI.log.prior
    maxratio.RI <- max(ratio.RI)
    marglik.RI <- -1*(maxratio.RI + log(mean(exp(ratio.RI-maxratio.RI))))

# Iterative Bridge Sampling, starting from IS and RI
# Step (c) of Algorithm 5.5 on page 152

    marglik.iter.BS <- matrix(0, maxiterbs, 2)
    marglik.iter.BS[1,1] <- marglik.IS
    marglik.iter.BS[1,2] <- marglik.RI

    iterbs <- 2

    # Bridge sampling loop
    cat('\n Begin Bridge Sampling \n')
    for (iterbs in 2:maxiterbs) {

  # Using IS for starting value
  # Nonnormalized mixture posterior
        logpostIS <- loglikIS + priorIS - marglik.iter.BS[iterbs-1,1]
        logpostRI <- loglikRI + priorRI - marglik.iter.BS[iterbs-1,1]
  # Ratio
        maxqRI <- apply(cbind(qIS, logpostIS), 1, max)  # N by 1
        rIS <- logpostIS - maxqRI - log(exp(log(N) + logpostIS -
                                            maxqRI) + exp(log(N) + qIS
                                                          - maxqRI))


        maxqRI <- apply(cbind(qRI, logpostRI), 1, max)  # N by 1
  # log(exp(... below is way off
        rRI <- qRI - maxqRI - log(exp(log(N) + logpostRI - maxqRI) +
                                  exp(log(N) + qRI - maxqRI))

  # max(rRI) and second log(mean()) is way off in line below
        marglik.iter.BS[iterbs,1] <- marglik.iter.BS[iterbs-1,1] +
            max(rIS) + log(mean(exp(rIS-max(rIS)))) - max(rRI) -
                log(mean(exp(rRI-max(rRI))))

  # Using RI for starting value

  # Nonnormalized mixture posterior
        logpostIS <- loglikIS + priorIS - marglik.iter.BS[iterbs-1,2]
        logpostRI <- loglikRI + priorRI - marglik.iter.BS[iterbs-1,2]
  # Ratio
        maxqRI <- apply(cbind(qIS, logpostIS), 1, max)  # N by 1

        rIS <- logpostIS - maxqRI - log(exp(log(N) + logpostIS -
                                            maxqRI) + exp(log(N) + qIS - maxqRI) )

        maxqRI <- apply(cbind(qRI, logpostRI), 1, max)  # N by 1
        rRI <- qRI - maxqRI - log(exp(log(N) + logpostRI - maxqRI) +
                                  exp(log(N) + qRI - maxqRI))

        marglik.iter.BS[iterbs,2] <- marglik.iter.BS[iterbs-1,2] +
            max(rIS) + log(mean(exp(rIS-max(rIS)))) - max(rRI) -
                log(mean(exp(rRI-max(rRI))))

        if (iterbs %% 100 == 0) cat('   Bridge Sampling: ', iterbs, '\n', sep='')
    }
    cat(' End Bridge Sampling \n')

# Choose to use the estimate using RI for final value
    marglik.BS <- marglik.iter.BS[nrow(marglik.iter.BS),2]
    marglik.iter.BS[nrow(marglik.iter.BS),1]


# Bridge sampling using importance sampling
    psIS <- loglikIS + priorIS
    mIS  <- max(qIS)
    fISbs <- psIS - mIS - log(exp(qIS - mIS) + exp(psIS - marglik.BS - mIS))

# Bridge sampling using reciprocal IS
    psRI <- loglikRI + priorRI
    mRI  <- max(qRI)
    fRIbs <- psRI - mRI - log(exp(qRI - mRI) + exp(psRI - marglik.BS - mRI))

    r <- 50

    h <- matrix(0, N, 2)
    h[,1] <- fISbs; h[,1] <- exp(h[,1] - max(h[,1]))
    h[,2] <- fRIbs; h[,2] <- exp(h[,2] - max(h[,2]))

    tmp <- ryan.autocov(h,r,2,N)
    hm <- tmp$hm; varhhat <- tmp$varhhat

    seblock <- t(((diag(varhhat)^.5)/hm))
    psi <- c(hm[1], -hm[2])
    sechib <- ((t(1/psi) %*% varhhat %*% (1/psi))^.5)

    marglik.BS.se <- sechib

# Standard error for Importance Sampling
    fIS <- loglikIS + priorIS - qIS

    h <- matrix(0, N, 1)
    h[,1] <- fIS; h[,1] <- exp(h[,1] - max(h[,1]))

    tmp <- ryan.autocov(h,r,1,N)
    hm <- tmp$hm; varhhat <- tmp$varhhat

    seblock <- t(((diag(varhhat)^.5)/hm))
    sechib <- ((t(1/hm) %*% varhhat %*% (1/hm))^.5)
    marglik.IS.se <- sechib

# Standard error for Reciprocal IS
    fRI <- qRI - loglikRI - priorRI

    h <- matrix(0, N, 1)
    h[,1] <- fRI; h[,1] <- exp(h[,1] - max(h[,1]))

    tmp <- ryan.autocov(h,r,1,N)
    hm <- tmp$hm; varhhat <- tmp$varhhat

    seblock <- t(((diag(varhhat)^.5)/hm))
    sechib <- ((t(1/hm) %*% varhhat %*% (1/hm))^.5)
    marglik.RI.se <- sechib

    output <- list(marglik.IS=marglik.IS,
                   marglik.IS.se=marglik.IS.se,
                   marglik.RI=marglik.RI,
                   marglik.RI.se=marglik.RI.se,
                   marglik.BS=marglik.BS,
                   marglik.BS.se=marglik.BS.se)

    class(output) <- 'posterior.fit.MSBVAR'
    return(output)
}




############################################
# Compute standard errors
############################################

# Create function for autocovariance
ryan.autocov <- function(h,r,k,n) {

    hm <- colMeans(h)
    hhat <- matrix(0, k, n)
    for (iterk in 1:k) { hhat[iterk,] <- hm[iterk] }

    omega <- array(0, c(k,k,r+1))
    omegat <- omega
    svec <- array(0, c(k,k,r+1))

    h <- t(h)
    for (s in 1:(r+1)) {
        tmpom1 <- matrix((h[,s:n]-hhat[,1:(n-s+1)]), nrow=k, ncol=(n-s+1))
        tmpom2 <- matrix(t(h[,1:(n-s+1)]-hhat[,1:(n-s+1)]), nrow=(n-s+1), ncol=k)
        omegas <- matrix((tmpom1 %*% tmpom2)/n, k, k)
        omega[,,s] <- omegas
        omegat[,,s] <- t(omegas)
        svec[,,s] <- (1-(s-1)/(r+1))/n
    }
    svec[,,1] = svec[,,1]*0.5

    sd <- (diag(matrix(omega[,,1], k, k)))^.5

    tmp <- svec*(omega+omegat)
    tmpsum <- matrix(0,k,k)
    for (itertmp in 1:(r+1)) { tmpsum <- tmpsum + tmp[,,itertmp] }
    varhhat <- tmpsum

    return(list(hm=hm, varhhat=varhhat))

} # end autocov function
