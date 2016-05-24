#' Construct prior from calibration model estimates
#'
#' Builds prior from appropriately structured output of the calibration model
#' from Blocker & Airoldi (2011). Handles all formatting so result can be fed
#' directly to \code{\link{bayesianDynamicFilter}}.
#'
#' @param xHat matrix (n x k) of estimates for OD flows from calibration model,
#'      one time point per row
#' @param varHat matrix (n x k) of estimated variances for OD flows from
#'      calibration, one time point per row
#' @param phiHat numeric vector (length n) of estimates for phi from calibration
#'      model
#' @param Y matrix (n x l) of observed link loads, one time point per row
#' @param A routing matrix (l x k) for network; must be of full row rank
#' @param phiPriorDf numeric prior convolution parameter for independent
#'      inverse-gamma priors on phi_t
#' @param rho numeric fixed autoregressive parameter for dynamics on lambda; see
#'      reference for details
#' @param backward logical to activate construction of reversed prior (for
#'      smoothing applications)
#' @param lambdaMin numeric value at which to floor estimated OD flows for prior
#'      construction
#' @param ipfp.maxit integer maximum number of iterations for IPFP 
#' @param ipfp.tol numeric tolerance for convergence of IPFP iterations
#' @return list containing priors for lambda and phi, consisting of:
#'      \itemize{
#'          \item mu, a matrix (n x k) containing the prior means for the
#'              log-change in each lambda at each time
#'          \item sigma, a matrix (n x k) containing the prior standard
#'              deviations for the log-change in each lambda at each time
#'          \item a list phi, containing the numeric prior \code{df} and a
#'              numeric vector \code{scale} of length n
#'      }
#' @keywords models multivariate ts
#' @references A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
#' time series on a graph. Proceedings of the Twenty-Seventh Conference Annual
#' Conference on Uncertainty in Artificial Intelligence (UAI-11) 51-60, 2011.
#' @export
#' @family bayesianDynamicModel
buildPrior <- function(xHat, varHat, phiHat, Y, A, rho=0.9,
                       phiPriorDf=ncol(A)/2, backward=FALSE, lambdaMin=1,
                       ipfp.maxit=1e6, ipfp.tol=1e-6) {
    # Extract dimensions of problem
    n <- nrow(Y)
    k <- ncol(A)
    l <- nrow(A)

    # Floor estimates at lambdaMin for stability on log scale
    lambdaPrior <- apply(xHat, 2, pmax, lambdaMin)

    # Apply IPFP to prior to ensure stability (issues arise when large negative
    # estimates are floored)
    for (i in 1:n) {
        tmp <- ipfp(Y[i,], A, lambdaPrior[i,], maxit=ipfp.maxit, tol=ipfp.tol)
        lambdaPrior[i,] <- tmp
    }

    # Floor again at lambdaMin
    lambdaPrior <- pmax(lambdaPrior, lambdaMin)

    # Smooth lambdaPrior with running median
    lambdaPrior <- apply(lambdaPrior, 2, runmed, 5)

    # Reverse if requested
    if (backward) {
        lambdaPrior <- lambdaPrior[n:1,]
        varHat <- varHat[n:1,]
        phiHat <- phiHat[n:1]
    }

    # Setup priors
    prior <- list()

    # Compute prior mean for changes in log-intensities
    prior$mu <- matrix(0, nrow(lambdaPrior), k)
    prior$mu[1,] <- log(lambdaPrior[1,])
    prior$mu[-1,] <- log(tail(lambdaPrior, -1)) - rho*log(head(lambdaPrior, -1))

    # Compute prior standard deviations for changes in log-intensities
    prior$sigma <- {c(1,rep(sqrt(1-rho^2),nrow(lambdaPrior)-1))*
                    sqrt(log(1+varHat/lambdaPrior^2))}

    # Construct prior for phi
    prior$phi <- list()
    prior$phi$df <- phiPriorDf
    prior$phi$scale <- prior$phi$df*log(1+phiHat)

    # Output prior
    return(prior)
}


#' Move step of sample-resample-move algorithm for multilevel state-space model
#' 
#' Function to execute single MCMC-based move step for
#' \code{\link{bayesianDynamicFilter}}. This can use two types of stopping
#' rules: number of iterations or number of accepted moves for the X particles.
#' The former is used by default, but the latter adapts better to low acceptance
#' rates (sometimes with substantial computational cost). Most updates in this
#' algorithm are Metropolis-Hastings with customized proposals.
#'
#' @param y numeric vector (length l) of observed link loads
#' @param X matrix (m x k) of particles for OD flows, one particle per row, in
#'      pivoted order
#' @param tme integer time index currently used in estimation
#' @param lambda matrix (m x k) of particles for OD intensities, one particle
#'      per row, in pivoted order
#' @param phi numeric vector (length m) of particles for phi
#' @param lambdatm1 lambda matrix (m x k) of particles for OD intensities from
#'      previous time, one particle per row, in pivoted order
#' @param phitm1 numeric vector (length m) of particles for phi from previous
#'      time
#' @param prior list containing priors for hyperparameters; see
#'      \code{\link{bayesianDynamicFilter}} for details
#' @param A routing matrix (l x k) for network
#' @param A1_inv inverse of full-rank portion of routing matrix (l x l)
#' @param A2 remainder of routing matrix (l x k-l)
#' @param rho numeric fixed autoregressive parameter for dynamics on lambda; see
#'      reference for details
#' @param tau numeric fixed power parameter for variance structure on truncated
#'      normal noise; see reference for details
#' @param m integer number of particles
#' @param l integer number of observed link loads
#' @param k integer number of OD flows to infer
#' @param ndraws integer number of draws to perform (can be overriden by
#'      minAccepts)
#' @param minAccepts integer minimum number of acceptances before results are
#'      returned; activates alternative stopping rule if >= 1
#' @param verbose logical activates verbose diagnostic output
#' @return list containing updated values of X, lambda, and phi
#' @keywords models multivariate ts
#' @references A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
#' time series on a graph. Proceedings of the Twenty-Seventh Conference Annual
#' Conference on Uncertainty in Artificial Intelligence (UAI-11) 51-60, 2011.
#' @export
#' @family bayesianDynamicModel
move_step <- function(y, X, tme, lambda, phi,
                      lambdatm1, phitm1, prior,
                      A, A1_inv, A2,
                      rho, tau,
                      m=ncol(X), l=nrow(A1_inv), k=length(lambda),
                      ndraws=10, minAccepts=0,
                      verbose=FALSE) {

    # Setup X1 and X2 matrices (want each particle in a column)
    X1 <- t(X[,seq(1,l)])
    X2 <- t(X[,seq(l+1,k)])
    X <- t(X)

    X1prop <- X1
    X2prop <- X2
    Xprop <- X

    # Setup lambda matrices
    lambda <- t(lambda)
    lambdatm1 <- t(lambdatm1)

    # Initialize information for acceptance rates
    xAccepts <- matrix(0,m,k-l)
    lambdaAccepts <- matrix(0,m,k)
    phiAccepts <- rep(0,m)

    # Setup matrices for mean, std dev, and truncation adjustment for truncated
    # normal
    meanMat <- exp(lambda)
    sdMat <- sqrt( t( t(meanMat)^tau * (exp(phi^2) - 1) ) )
    truncMat <- matrix(pnorm(0, meanMat, sdMat, lower.tail=FALSE, log.p=TRUE),
                       k, m)

    # Loop for MCMC iterations
    iter <- 0
    while (iter < ndraws || (min(xAccepts) < minAccepts)) {
        # Draw lambda, phi | X via MH steps

        # MH step on lambda
        varPropMat <- 1 / (1/prior$sigma[tme,]^2 +
                           t(matrix(1/phi^2, m, k)))
        meanPropMat <- lambda

        lambdaProp <- matrix(rnorm(k*m, meanPropMat, sqrt(varPropMat)), k, m)

        varPropMatRev <- varPropMat
        meanPropMatRev <- lambdaProp

        # Calculate LLR and LIR
        meanMatProp <- exp(lambdaProp)
        sdMatProp <- sqrt( t( t(meanMatProp)^tau * (exp(phi^2) - 1) ) )
        truncMatProp <- matrix(pnorm(0, meanMatProp, sdMatProp,
                                     lower.tail=FALSE, log.p=TRUE),
                               k, m)

        llr <- matrix(dnorm(lambdaProp, prior$mu[tme,] + rho*lambdatm1,
                            prior$sigma[tme,], log=TRUE) +
                      dnorm(X, meanMatProp, sdMatProp, log=TRUE) - truncMatProp,
                      k, m)
        llr <- llr - matrix(dnorm(lambda, prior$mu[tme,] + rho*lambdatm1,
                                  prior$sigma[tme,], log=TRUE) +
                            dnorm(X, meanMat, sdMat, log=TRUE) - truncMat,
                            k, m)

        lir <- matrix(dnorm(lambdaProp, meanPropMat, sqrt(varPropMat),
                            log=TRUE),
                      k, m)
        lir <- lir - matrix(dnorm(lambda, meanPropMatRev, sqrt(varPropMatRev),
                                  log=TRUE),
                            k, m)

        # Accept with correct probability
        logAcceptProb <- (llr - lir)
        logAcceptProb[is.na(logAcceptProb)] <- -Inf

        acceptMat <- ( log(runif(k*m)) < logAcceptProb )

        # Handle acceptances and store necessary information
        lambda[acceptMat] <- lambdaProp[acceptMat]
        lambdaAccepts <- lambdaAccepts + t(acceptMat)

        # Update means and variances for X
        meanMat <- exp(lambda)
        sdMat <- sqrt( t( t(meanMat)^tau * (exp(phi^2) - 1) ) )
        truncMat <- matrix(pnorm(0, meanMat, sdMat,
                                 lower.tail=FALSE, log.p=TRUE),
                           k, m)

        ###

        # MH step on phi
        alphaProp <- (k/2 + prior$phi$df)
        betaProp <- colSums((log(X)-lambda)^2)/2 + prior$phi$scale[tme]
        phiProp <- 1 / rgamma(m, alphaProp, betaProp)
        phiProp <- sqrt(phiProp)

        # Calculate LLR & LIR
        meanMatProp <- exp(lambda)
        sdMatProp <- sqrt( t( t(meanMatProp)^tau * (exp(phiProp^2) - 1) ) )
        truncMatProp <- matrix(pnorm(0, meanMatProp, sdMatProp,
                                     lower.tail=FALSE, log.p=TRUE),
                               k, m)

        llr <- (dgamma(1/phiProp^2, prior$phi$df, prior$phi$scale[tme],
                       log=TRUE) - 2*log(phiProp^2) +
                colSums(matrix(dnorm(X, meanMatProp, sdMatProp, log=TRUE) -
                               truncMatProp, k, m)))
        llr <- llr - (dgamma(1/phi^2, prior$phi$df, prior$phi$scale[tme],
                             log=TRUE) - 2*log(phi^2) +
                      colSums(matrix( dnorm(X, meanMat, sdMat, log=TRUE)
                             - truncMat, k, m)))

        lir <- (dgamma(1/phiProp^2, alphaProp, betaProp, log=TRUE) -
                2*log(phiProp^2))
        lir <- lir - ( dgamma(1/phi^2, alphaProp, betaProp, log=TRUE) -
                      2*log(phi^2) )

        # Accept with correct probability
        logAcceptProb <- (llr - lir)
        logAcceptProb[is.na(logAcceptProb)] <- -Inf

        acceptVec <- ( log(runif(m)) < logAcceptProb )

        # Handle acceptances and store necessary information
        phi[acceptVec] <- phiProp[acceptVec]
        phiAccepts <- phiAccepts + acceptVec

        # Update means and variances for X
        meanMat <- exp(lambda)
        sdMat <- sqrt( t( t(meanMat)^tau * (exp(phi^2) - 1) ) )
        truncMat <- matrix(pnorm(0, meanMat, sdMat, lower.tail=FALSE,
                                 log.p=TRUE), k, m)

        ###

        # Draw X | lambda, phi

        # # Using coordinate direction algorithm of Smith (1984)
        # # Draw (uniform) random coordinate, then propose uniformly along
        # Using random direction algorithm of Smith (1984)
        # Draw (uniform) random direction, then propose uniformly along
        # feasible line in that direction
        # Symmetric proposal -> random walk MH

        for (j in 1:(k-l)) {
            # Draw random directions (uniform on unit sphere)
            # by normalizing MV normal RVs
            dMat <- matrix( rnorm((k-l)*m), k-l, m )
            dMat <- t( t(dMat) / sqrt(colSums(dMat*dMat)) ) 

            # Calculate relevant matrices for feasible region calculations
            remainderMat <- X1
            adjMat <- A1_inv %*% A2 %*% dMat
            adjMat.nonzero <- which(abs(adjMat)>.Machine$double.eps)
            interceptMat <- matrix(0, nrow=nrow(remainderMat),
                                   ncol=ncol(remainderMat))
            interceptMat[adjMat.nonzero] <- remainderMat[adjMat.nonzero] /
                                            adjMat[adjMat.nonzero]

            # Find vector of maxima for \delta X2 st X1 >= 0
            maxVec1 <- sapply(1:m, function(i)
                              if (max(interceptMat[,i])>0) 
                                  return(min(interceptMat[,i][adjMat[,i]>0]))
                              else
                                  return(0)
                              )
            # Find vector of maxima for \delta X2 st X2 >= 0
            maxVec2 <- sapply(1:m, function(i)
                              if (max(-X2[,i]/dMat[,i])>0) 
                                  return(min( (-X2[,i]/dMat[,i])[-dMat[,i]>0] ))
                              else
                                  return(0)
                              )
            # Merge maxima
            maxVec <- pmin(maxVec1, maxVec2)

            # Find vector of minima for \delta X2 st X1 >= 0
            minVec1 <- sapply(1:m, function(i)
                              if (min(interceptMat[,i])<0)
                                  return(max(interceptMat[,i][adjMat[,i]<0]))
                              else
                                  return(0)
                              )
            # Find vector of minima for \delta X2 st X2 >= 0
            minVec2 <- sapply(1:m, function(i)
                              if (min(-X2[,i]/dMat[,i])<0)
                                  return(max( (-X2[,i]/dMat[,i])[-dMat[,i]<0] ))
                              else
                                  return(0)
                              )
            # Merge minima
            minVec <- pmax(minVec1, minVec2)

            proposal <- runif(m, minVec, maxVec)
            X2prop <- X2 + t( t(dMat) * proposal )
            X1prop <- X1 - t( t(adjMat) * proposal )
            Xprop[1:l,] <- X1prop
            Xprop[(l+1):k,] <- X2prop

            llr <- colSums( matrix(dnorm(Xprop, meanMat, sdMat, log=TRUE) -
                                   dnorm(X, meanMat, sdMat, log=TRUE), k, m) )
            lir <- 0

            logAcceptProb <- rep(-Inf, m)
            validVec <- !is.na(llr-lir)
            logAcceptProb[validVec] <- (llr - lir)[validVec]

            acceptVec <- (log(runif(m)) < logAcceptProb)

            X[,acceptVec] <- Xprop[,acceptVec]

            xAccepts[,j] <- xAccepts[,j] + acceptVec
        }

        iter <- iter + 1
    }

    # If verbose, print acceptance rate summaries
    if (verbose) {
        print(iter)
        # print(colMeans(xAccepts))
        print(summary(rowSums(xAccepts)))
        print(colMeans(lambdaAccepts))
        print(mean(phiAccepts))
    }

    # Return updated results
    return(list( X=t(X), lambda=t(lambda), phi=phi ))
}

#' Function for inference with multilevel state-space model
#'
#' Particle filtering with sample-resample-move algorithm for multilevel
#' state-space model of Blocker & Airoldi (2011). This has log-normal
#' autoregressive dynamics on OD intensities, log-normal emission distributions,
#' and truncated normal observation densities.  This can return full (all
#' particles) output, but it is typically better to aggregate results as you go
#' to reduce memory consumption. It can also run forward or backward filtering
#' for smoothing. These results are combined via a separate function for
#' smoothing; however, this procedure typically performs poorly due to
#' differences between the distributions of particles from forward and reverse
#' filtering.
#'
#' @param Y matrix (n x l) of observed link loads over time, one observation per
#'      row
#' @param A routing matrix (l x k) for network; must be of full row rank
#' @param prior list containing priors for lambda and phi; must have
#'      \itemize{
#'          \item mu, a matrix (n x k) containing the prior means for the
#'              log-change in each lambda at each time
#'          \item sigma, a matrix (n x k) containing the prior standard
#'              deviations for the log-change in each lambda at each time
#'          \item a list phi, containing the numeric prior \code{df} and a
#'              numeric vector \code{scale} of length n
#'      }
#' @param lambda0 numeric vector (length k) of time 0 prior means for OD flows
#' @param sigma0 numeric vector (length k) of time 0 prior standard deviations
#'      for OD flows
#' @param phi0 numeric starting value for phi at time 0
#' @param rho numeric fixed autoregressive parameter for dynamics on lambda; see
#'      reference for details
#' @param tau numeric fixed power parameter for variance structure on truncated
#'      normal noise; see reference for details
#' @param m integer number of particles to use
#' @param verbose logical activates verbose diagnostic output
#' @param Xdraws integer number of draws to perform for \code{xsample} RDA
#       proposals; must be larger than \code{m}
#' @param Xburnin integer number of burnin draws to discard for \code{xsample}
#'      proposals RDA in addition to baseline number of draws
#' @param Movedraws integer number of iterations to run for each move step
#' @param nThresh numeric effective number of independent particles below which
#'      redraw will be performed
#' @param aggregate logical to activate aggregation of MCMC results; highly
#       recommended for large networks
#' @param backward logical to activate reverse filtering (for smoothing
#       applications)
#' @param tStart integer time index to begin iterations from
#' @return list containing:
#'      \itemize{
#'          \item xList
#'          \item lambdaList
#'          \item phiList
#'          \item y
#'          \item rho
#'          \item prior
#'          \item n
#'          \item l
#'          \item k
#'          \item A
#'          \item A_qr
#'          \item A1
#'          \item A1_inv
#'          \item A2
#'          \item nEff
#'          \item tStart
#'          \item backward
#'          \item aggregate
#'      }
#' @references A.W. Blocker and E.M. Airoldi. Deconvolution of mixing time
#' series on a graph. Proceedings of the Twenty-Seventh Conference Annual
#' Conference on Uncertainty in Artificial Intelligence (UAI-11) 51-60, 2011.
#' @export
#' @family bayesianDynamicModel
bayesianDynamicFilter <- function (Y, A, prior,
                                   lambda0, sigma0, phi0, 
                                   rho=0.1, tau=2, m=1e3,
                                   verbose=FALSE,
                                   Xdraws=5*m, Xburnin=m, Movedraws=10,
                                   nThresh=10,
                                   aggregate=FALSE,
                                   backward=FALSE,
                                   tStart=1) {
    # Calculate dimensions
    k <- ncol(A)
    l <- nrow(A)
    n <- nrow(Y)

    # Decompose A matrix into full-rank and singular parts; retain pivot info
    A_qr <- qr(A)
    pivot <- A_qr$pivot
    A_pivot <- A[,pivot]
    A1 <- A_pivot[,seq(A_qr$rank)]
    A1_inv <- solve(A1)
    A2 <- A_pivot[,seq(A_qr$rank+1,k)]

    # If verbose, print pivot information
    if (verbose) {
        cat('Pivot information:\n', file=stderr())
        print(pivot)
        cat(sprintf('rank(A) = %d\n', A_qr$rank), file=stderr())
    }

    # Reconfigure prior based on pivot (only using components pertaining to
    # X2, as X1 is determined based on data)
    prior$mu <- prior$mu[,pivot]
    prior$sigma <- prior$sigma[,pivot]

    # Reconfigure initial values based on pivot
    lambda0 <- lambda0[pivot]
    sigma0 <- sigma0[pivot]

    # If running backwards filter, reverse relevant matrices
    # Assumes prior has been input reversed; only reversing y
    if (backward) {
        Y <- Y[n:1,]
    }

    # Setup data structure for estimated OD flows
    xList <- list()
    lambdaList <- list()
    phiList <- list()
    nEffVec <- rep(0, n)

    # Setup data structures for particle filter
    lambda <- matrix(NA, m, k)
    w <- numeric(m)

    # Initialize particles from prior at time 0
    lambda <- matrix( rnorm(k*m, lambda0, sigma0), m, k, byrow=TRUE)
    phi <- rep(phi0, m)

    # Iterate over times
    for (tme in tStart:n) {
        # If verbose, print time
        if (verbose) {
            cat(sprintf('t = %d\n', tme), file=stderr())
            startTime <- proc.time()
        }

        # Retain last set of particles for move step
        lambdatm1 <- lambda
        phitm1 <- phi

        # Draw phi from prior
        phi_prop <- 1/rgamma(m, prior$phi$df, prior$phi$scale[tme])
        phi_prop <- sqrt(phi_prop)

        # Draw set of lambdas
        lambda_prop <- (rho*lambdatm1 +
                        matrix(rnorm(k*m, prior$mu[tme,], prior$sigma[tme,]),
                               m, k, byrow=TRUE) )

        # Draw X's from truncated normal on feasible region via xsample, 
        # then thin
        propMean <- exp( rho*colMeans(lambdatm1) + prior$mu[tme,] )
        propSD <- sqrt( (exp(prior$phi$scale[tme]/prior$phi$df)-1) *
                       propMean^tau )

        # Check for known OD flows
        odRanges <- xranges(E=A_pivot, F=Y[tme,], ispos=TRUE)
        activeOD <- (odRanges[,2]-odRanges[,1]>0)
        #
        A.active <- A_pivot[,activeOD]
        At.active.qr <- qr(t(A.active))
        if (At.active.qr$rank < nrow(A.active)) {
            knownOD <- !activeOD

            activeLink <- At.active.qr$pivot[1:At.active.qr$rank]
            A.active <- A.active[activeLink,,drop=FALSE]
            
            x.adj <- rep(0, ncol(A_pivot))
            x.adj[knownOD] <- odRanges[knownOD, 2]
            Y.active <- Y[tme,] - A_pivot %*% x.adj
            Y.active <- Y.active[activeLink]
        } else {
            activeLink <- 1:nrow(A_pivot)
            Y.active <- Y[tme,activeLink]
        }
        #
        nActive <- sum(activeOD)

        if (verbose) {
            cat(sprintf("nActive = %d\n", nActive), file=stderr())
            if (verbose > 1)
                print(activeLink)
        }

        # Start with IPF beginning from propMean
        x0Active <- ipfp(Y.active, A.active, propMean[activeOD], tol=1e-6,
                         maxit=Xdraws*1e2, verbose=FALSE)

        # Now, lsei for refinement
        x0Active <- lsei(A=diag(nActive), B=x0Active, E=A.active,
                         F=Y.active, G=diag(nActive), H=rep(0,nActive),
                         tol=sqrt(.Machine$double.eps), type=1)$X

        # Switch to RDA algorithm
        # Correct distribution (truncated normal) and faster
        X_prop_active <- xsample(A=diag(nActive), B=propMean[activeOD],
                                 E=A.active, F=Y.active,
                                 G=diag(nActive), H=rep(0,nActive),
                                 sdB=propSD[activeOD], iter=Xdraws,
                                 outputlength=m, burninlength=Xburnin,
                                 x0=x0Active, type='rda')

        if (verbose) {
            cat(sprintf("Accepted ratio = %g\n", X_prop_active$acceptedratio),
                file=stderr())
        }

        if (nActive == k) {
            X_prop <- X_prop_active$X
        } else {
            X_prop <- matrix(0, m, k)
            X_prop[,activeOD] <- X_prop_active$X
            #
            X_prop <- t(X_prop)
            X_prop[!activeOD,] <- odRanges[!activeOD,2]
            X_prop <- t(X_prop)
        }

        # Calculate weights
        w <- rep(0,m)

        # Drawing lambda from prior; drawing X from truncated normal on
        # feasible region
        meanMat <- exp(lambda_prop)
        sdMat <- sqrt( meanMat^tau * (exp(phi_prop^2) - 1) )
        truncMat <- matrix(pnorm(0, meanMat, sdMat, lower.tail=FALSE,
                                 log.p=TRUE), m, k)

        w <- w + rowSums(matrix(dnorm(X_prop, meanMat, sdMat, log=TRUE) -
                                truncMat, m, k))
        w <- w - colSums(matrix(dnorm( t(X_prop), propMean, propSD, log=TRUE),
                                k, m))

        w <- exp(w - max(w))
        # Handle numerical errors
        w[is.na(w)] <- 0

        # Normalize importance weights
        w <- w / sum(w)

        # Check effective number of particles
        nEff <- 1/sum(w*w)

        # If verbose, print importance weight diagnostics
        if (verbose) {
            cat(sprintf('Effective number of particles: %g\n', nEff),
                file=stderr())
            if (nEff < nThresh)
                cat('Redrawing lambda from prior\n', file=stderr())
        }

        # If too few particles, redraw lambda from prior
        if (nEff < nThresh && tme > 1) {
            meanPrior <- colSums( rho^(tme:1 - 1) * prior$mu[1:tme,] )
            sdPrior <- 1/(1-rho^2) * prior$sigma[tme,]

            lambda_prop <- matrix(rnorm(k*m, colSums(rho^(tme:1 - 1) *
                                                     prior$mu[1:tme,] ),
                                        1/(1-rho^2) * prior$sigma[tme,]), m, k,
                                  byrow=TRUE )

            # Draw X's from truncated normal on feasible region via xsample, 
            # then thin
            propMean <- exp( meanPrior )
            propSD <- sqrt( (exp(prior$phi$scale[tme]/prior$phi$df)-1) *
                           propMean^tau )

            # Start with IPF beginning from propMean
            x0Active <- ipfp(Y.active, A.active, propMean[activeOD], tol=1e-6,
                             maxit=Xdraws*1e2, verbose=FALSE)

            # Now, lsei for refinement
            x0Active <- lsei(A=diag(nActive), B=x0Active, E=A.active,
                             F=Y.active, G=diag(nActive), H=rep(0,nActive),
                             tol=sqrt(.Machine$double.eps), type=1)$X

            # Switch to RDA algorithm
            # Correct distribution (truncated normal) and faster
            X_prop_active <- xsample(A=diag(nActive), B=propMean[activeOD],
                                     E=A.active, F=Y.active,
                                     G=diag(nActive), H=rep(0,nActive),
                                     sdB=propSD[activeOD], iter=Xdraws,
                                     outputlength=m, burninlength=Xburnin,
                                     x0=x0Active, type='rda')

            if (verbose) {
                cat(sprintf("Accepted ratio = %g\n",
                            X_prop_active$acceptedratio),
                    file=stderr())
            }

            if (nActive == k) {
                X_prop <- X_prop_active$X
            } else {
                X_prop <- matrix(0, m, k)
                X_prop[,activeOD] <- X_prop_active$X
                #
                X_prop <- t(X_prop)
                X_prop[!activeOD,] <- odRanges[!activeOD,2]
                X_prop <- t(X_prop)
            }

            # Drawing lambda from prior; drawing X from truncated normal on
            # feasible region
            meanMat <- exp(lambda_prop)
            sdMat <- sqrt( meanMat^tau * (exp(phi_prop^2) - 1) )
            truncMat <- matrix(pnorm(0, meanMat, sdMat, lower.tail=FALSE,
                                     log.p=TRUE), m, k)

            w <- rep(0, m)
            w <- w + rowSums(matrix(dnorm(X_prop, meanMat, sdMat, log=TRUE) -
                                    truncMat, m, k))
            w <- w - colSums(matrix(dnorm(t(X_prop), propMean, propSD,
                                          log=TRUE), k, m))

            w <- exp(w - max(w))

            # Handle numerical errors
            w[is.na(w)] <- 0

            # Normalize importance weights
            w <- w / sum(w)

            # Check effective number of particles
            nEff <- 1/sum(w*w)

            if (verbose) {
                cat(sprintf('Effective number of particles after redraw: %g\n',
                            nEff))
            }
        }

        # Resample particles according to weights
        ind <- sample(m, m, replace=TRUE, prob=w)

        lambda <- lambda_prop[ind,]
        phi <- phi_prop[ind]
        X <- X_prop[ind,]

        if (verbose > 1) {
            # Check constraints
            yHat <- A_pivot %*% t(X)
            err.l1 <- apply(yHat-Y[tme,], 2, function(x)
                            sum(abs(x)))
            cat("L_1 deviation from constraints\n", file=stderr())
            print(summary(err.l1))
        }

        # Move step
        move.out <- move_step(y=Y[tme,], X=X, tme, lambda=lambda, phi=phi,
                              lambdatm1=lambdatm1, phitm1=phitm1, prior=prior,
                              A=A_pivot, A1_inv=A1_inv, A2=A2, rho=rho, tau=tau,
                              m=m, l=l, k=k, ndraws=Movedraws, verbose=verbose)

        # Get revised results
        lambda <- move.out$lambda
        X <- move.out$X
        phi <- move.out$phi

        # Store results from this iteration
        nEffVec[tme] <- nEff

        if (aggregate) {
            xList[[tme]] <- agg(X)
            lambdaList[[tme]] <- agg(lambda)
            phiList[[tme]] <- agg(phi)
        } else {
            xList[[tme]] <- X
            lambdaList[[tme]] <- lambda
            phiList[[tme]] <- phi
        }

        # If verbose, print parameter summaries
        if (verbose) {
            endTime <- proc.time()
            print(mean(phi))
            print( t(data.frame(lambdaHat=exp(colMeans(lambda)),
                                xMean=colMeans(X), xMed=apply(X,2,median),
                                xSD=apply(X,2,sd))) )
            cat(sprintf("Runtime for iteration %d:\n", tme), file=stderr())
            print(endTime-startTime)
        }
    }

    # Un-reverse results, if needed
    if (backward) {
        Y <- Y[n:1,]
        xList <- xList[n:1]
        lambdaList <- lambdaList[n:1]
        phiList <- phiList[n:1]
    }

    # Return list containing results
    unpivot <- numeric(k)
    unpivot[pivot] <- seq(k)

    xList <- lapply(xList, function(mat) mat[,unpivot])
    lambdaList <- lapply(lambdaList, function(mat) mat[,unpivot])

    retval <- list(xList=xList, lambdaList=lambdaList, phiList=phiList,
                   y=Y, rho=rho, prior=prior, n=n, l=l, k=k,
                   A=A, A_qr=A_qr, A1=A1, A1_inv=A1_inv, A2=A2,
                   nEff=nEffVec,
                   tStart=tStart, backward=backward, aggregate=aggregate)

    return(retval)
}
