
#' Objective function of Zhang et al. 2003
#' 
#' Requires bounded optimization to maintain positive OD flows, and only those 
#' flows that are not deterministically zero should be included in the 
#' estimation.
#' 
#' @param xt length-k numeric vector of point-to-point flows
#' @param yt length-m numeric vector of observed aggregate flows
#' @param A m x k routing matrix, yt = A xt
#' @param srcDstInd list of source and destination flow indices corresponding to
#'   each point-to-point flow, as produced by \code{\link{getSrcDstIndices}}
#' @param lambda regularization parameter for mutual information prior. Note 
#'   that this is scaled by the squared total traffic in the objective function 
#'   before scaling the mututal information prior.
#' @return numeric value of objective function to minimize in tomogravity 
#'   estimation
obj.tomogravity <- function(xt, yt, A, srcDstInd, lambda) {
    # Compute fitted values
    yHat <- A %*% xt
    pos <- which(xt > sqrt(.Machine$double.eps))
    N <- sum(yt) / 2
    g <- yHat[srcDstInd$src] * yHat[srcDstInd$dst] / sum(xt)

    # Compute mutual information
    mi  <- sum(xt[pos]/sum(xt) * log(xt[pos]/g[pos]))

    # Compute overall objective function
    return(sum((yt-yHat)^2) + lambda * N^2 * mi)
}

#' Analytic gradient of objective function of Zhang et al. 2003
#' 
#' Requires bounded optimization to maintain positive OD flows, and only those 
#' flows that are not deterministically zero should be included in the 
#' estimation.
#' 
#' @param xt length-k numeric vector of point-to-point flows
#' @param yt length-m numeric vector of observed aggregate flows
#' @param A m x k routing matrix, yt = A xt
#' @param srcDstInd list of source and destination flow indices corresponding to
#'   each point-to-point flow, as produced by \code{\link{getSrcDstIndices}}
#' @param lambda regularization parameter for mutual information prior. Note 
#'   that this is scaled by the squared total traffic in the objective function 
#'   before scaling the mututal information prior.
#' @return numeric vector of length k containing gradient of objective function 
#'   with respect to xt
dobj.dxt.tomogravity <- function(xt, yt, A, srcDstInd, lambda) {
    # Compute fitted values
    yHat <- A %*% xt
    pos <- which(xt > sqrt(.Machine$double.eps))
    N <- sum(xt)
    g <- yHat[srcDstInd$src] * yHat[srcDstInd$dst] / N
    
    # Compute derivative of quadratic part
    dd <- -2 * t(A) %*% (yt - yHat)

    # Compute derivative of mutual information part
    dmi <- rep(0, length(xt))
    dmi[pos] <- log(xt[pos]/g[pos])/N + 1/N

    # Combine parts of derivative
    dd <- dd + lambda * N^2 * dmi
    
    return(dd)
}

#' Tomogravity estimation for a single time point using L-BFGS-B
#' 
#' @param yt length-m numeric vector of observed aggregate flows at time t
#' @param A m x k routing matrix
#' @param srcDstInd list of source and destination flow indices corresponding to
#'   each point-to-point flow, as produced by \code{\link{getSrcDstIndices}}
#' @param lambda regularization parameter for mutual information prior. Note 
#'   that this is scaled by the squared total traffic in the objective function 
#'   before scaling the mututal information prior.
#' @param N total traffic for normalization. Unused if normalized is FALSE.
#' @param normalize If TRUE, xt and yt are scaled by N. Typically used in 
#'   conjunction with calcN to normalize traffic to proportions, easing the 
#'   tuning of lambda.
#' @param lower Component-wise lower bound for xt in L-BFGS-B optimization.
#' @param control List of control information for optim.
#' @return A list as returned by optim, with element \code{par} containing the
#'   estimated point-to-point flows and elementer \code{gr} containing the
#'   analytic gradient evaluated at the estimate.
#' @keywords models multivariate ts
#' @export
#' @family tomogravity
#' @examples
#' data(cmu)
#' srcDstInd <- getSrcDstIndices(cmu$A.full)
#' estimate <- tomogravity.fit(yt=cmu$Y.full[1, ], A=cmu$A.full,
#'      srcDstInd=srcDstInd, lambda=0.01)
tomogravity.fit <- function(yt, A, srcDstInd, lambda, N=1, normalize=FALSE,
                            lower=0, control=list()) {
    # Store dimensions
    k <- ncol(A)
    l <- nrow(A)
    
    # Normalize if requested
    if (normalize)
        norm <- N
    else
        norm <- 1

    ## Find positive flows
    #pos <- which(getActive(yt, A))
    #nActive <- length(pos)
    pos <- seq(ncol(A))
    nActive <- ncol(A)
    
    # Initialize via IPFP
    xt.init <- ipfp(yt, A[, pos], rep(1, nActive))

    # Run optimization
    result <- optim(par=xt.init/norm, fn=obj.tomogravity, gr=NULL,
                    yt=yt/norm, A=A[, pos], srcDstInd=srcDstInd,
                    lambda=lambda,
                    method="L-BFGS-B",
                    lower=rep(lower, nActive),
                    control=control)
    
    # Setup return value
    init <- rep(0, k)
    init[pos] <- xt.init
    #
    xt.hat <- init
    xt.hat[pos] <- result$par

    # Unnormalize if needed
    xt.hat <- xt.hat * norm

    retval <- c(result, list(xt.init=init))
    retval$par <- xt.hat
    retval$gr <- dobj.dxt.tomogravity(result$par, yt=yt/norm, A=A[, pos],
                                      srcDstInd=srcDstInd, lambda=lambda)
    return(retval)
}

#' Run tomogravity estimation on complete time series of aggregate flows
#' 
#' The aggregate flows Y and their corresponding routing matrix A must include 
#' all aggregate source and destination flows.
#' 
#' @param Y n x m matrix contain one vector of observed aggregate flows per row.
#'   This should include all observed aggegrate flows with none removed due to 
#'   redundancy.
#' @param A m x k routing matrix. This need not be of full row rank and must 
#'   include all source and destination flows.
#' @param lambda Regularization parameter for mutual information prior. Note 
#'   that this is scaled by the squared total traffic in the objective function 
#'   before scaling the mututal information prior.
#' @param lower Component-wise lower bound for xt in L-BFGS-B optimization.
#' @param normalize If TRUE, xt and yt are scaled by N. Typically used in 
#'   conjunction with calcN to normalize traffic to proportions, easing the 
#'   tuning of lambda.
#' @param .progress name of the progress bar to use, see 
#'   \code{\link{create_progress_bar}} in plyr documentation
#' @param control List of control information for optim.
#' @return A list containing three elements: \itemize{ \item resultList, a list 
#'   containing the output from running \code{\link{tomogravity.fit}} on each 
#'   timepoint \item changeFromInit, a vector of length n containing the 
#'   relative L_1 change between the initial (IPFP) point-to-point flow 
#'   estimates and the final tomogravity estimates \item Xhat, a n x k matrix 
#'   containing a vector of estimated point-to-point flows (for each time point)
#'   per row }
#' @keywords models multivariate ts
#' @export
#' @family tomogravity
#' @examples
#' data(cmu)
#' estimate <- tomogravity(Y=cmu$Y.full[1, , drop=FALSE], A=cmu$A.full,
#'                         lambda=0.01, .progress='text')
tomogravity <- function(Y, A, lambda, lower=0, normalize=FALSE,
                        .progress="none", control=list()) {
    # Decompose routing matrix A
    Adecomp <- decomposeA(A)

    # Calculate total traffic for each time
    Nvec <- apply(Y, 1, calcN, A1=Adecomp$A1)

    # Calculate independence flows for each time
    srcDstInd <- getSrcDstIndices(A)

    # Run tomogravity inference for each time
    resultList <- llply(seq(nrow(Y)), function(tme)
                        tomogravity.fit(yt=Y[tme, ], A=A,
                                         srcDstInd=srcDstInd, lambda=lambda,
                                         normalize=normalize,
                                         lower=lower,
                                         control=control),
                        .progress=.progress)

    # Compute change from initalization as diagnostic
    changeFromInit <- sapply(resultList, function(result)
                             sum(abs(result$par - result$xt.init)) /
                             sum(result$par))

    # Extract results into nice matrix
    Xhat <- t(sapply(resultList, function(result) result$par))
    
    # Return results
    retval <- list(resultList=resultList,
                   changeFromInit=changeFromInit,
                   Xhat=Xhat)
    return(retval)
}

