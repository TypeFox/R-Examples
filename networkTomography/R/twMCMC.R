
#' Function to run MCMC sampling for model of Tebaldi & West (1998)
#'
#' Runs MCMC sampling for the gamma-Poisson model presented in Tebaldi & West
#' (1998). The algorithm used is a modification of that presented in the
#' original paper. It uses a joint proposal for (x_k, lambda_k) to greatly
#' accelerate convergence.
#'
#' @param Y numeric vector of observed link loads at a single time (length k)
#' @param A routing matrix of dimension (k x n); needs to be full row rank
#' @param prior parameters for conjugate gamma prior (convolution and rate)
#' @param ndraws integer number of draws for sampler to produce (excluding
#'  burn-in)
#' @param burnin integer number of additional draws to discard as burnin
#' @param verbose integer level of verbosity; levels > 1 have no effect
#'  currently
#' @return list consisting of matrix of draws for X \code{XDraws},
#'  matrix of draws for X \code{lambdaDraws}, and vector of acceptances per OD
#'  flow \code{accepts}
#' @references C. Tebaldi and M. West. Bayesian inference on network traffic
#' using link count data. Journal of the American Statistical Association,
#' 93(442):557-573, 1998.
#' @keywords models multivariate
#' @export
#' @examples
#' data(bell.labs)
#' # Quick, simple run to test the function
#' prior <- list(a=rep(1, ncol(bell.labs$A)), b=rep(0, ncol(bell.labs$A)))
#' mcmcOut <- twMCMC(Y=bell.labs$Y[1,], A=bell.labs$A, prior=prior,
#'                   ndraws=1000, burnin=100,
#'                   verbose=0)
#' print(summary(mcmcOut$XDraws))
#' print(mcmcOut$accepts)
twMCMC <- function(Y, A, prior, ndraws=1.2e5, burnin=2e4, verbose=0) {
    # Format verification
    Y <- as.numeric(Y)

    # Initialize X and lambda
    n <- ncol(A)
    k <- nrow(A)
    X <- matrix(0, ndraws+burnin, n)

    # Decompose A matrix into full-rank and singular parts; retain pivot info
    A_qr <- qr(A)
    pivot <- A_qr$pivot
    A_pivot <- A[,pivot]
    A1 <- A_pivot[,seq(A_qr$rank)]
    A1_inv <- solve(A1)
    A2 <- A_pivot[,seq(A_qr$rank+1,n)]

    # Pivot prior
    prior$a <- prior$a[pivot]
    prior$b <- prior$b[pivot]

    # Initialize X via simplex method (gives integer solutions with linear
    # integer constraints for transportation problems)
    # New and improved! Now gives actual interior integer solutions!
    bounds <- xranges(E=A_pivot, F=Y, ispos=TRUE)
    #
    obj <- rnorm(n)
    const_mat <- rbind(A_pivot, diag(n))
    const_dir <- c(rep('==',k), rep('>=',n))
    int_vec <- seq(n)

    const_rhs <- c(Y, round(bounds[,'max']/2))
    #
    tries <- 0
    status <- 1
    while ( status > 0 ) {
        init <- Rglpk_solve_LP(obj, const_mat, const_dir, const_rhs)
        status <- init$status
        const_rhs[(k+1):(n+k)] <- round(const_rhs[(k+1):(n+k)]/2)
        tries <- tries + 1
        #
        if (status == 0) {
            x.init <- init$solution
            err <- A_pivot %*% x.init - Y
            if (sum(abs(err)) > 0) {
                # In event of FUBAR, restart
                obj <- rnorm(n)
                const_rhs[(k+1):(n+k)] <- round(bounds[,'max']/2)
            }
        }
    }

    X[1,] <- x.init
    X1 <- X[1,seq(k)]
    X2 <- X[1,seq(k+1,n)]

    # If verbose, print initial solution
    if (verbose) print(X[1,])

    # Initialize lambda with draw from conditional posterior
    lambda <- matrix(0, ndraws+burnin, n)
    lambda[1,] <- rgamma(n, prior$a+X[1,], prior$b+1)

    # Setup adjList
    adjList <- lapply(1:(n-k), function(j) A1_inv %*% A2[,j])
    posList <- lapply(adjList, function(x) which(x>0))
    negList <- lapply(adjList, function(x) which(x<0))

    # Loop for MCMC iterations
    accepts <- rep(0, n-k)
    for (iter in seq(2,ndraws+burnin)) {
        # Draw X2_j, lambda_j+k
        for (j in seq(n-k)) {
            # MH step on X2_j

            # Propose uniformly along feasible region
            adjVec <- adjList[[j]]
            posVec <- posList[[j]]
            negVec <- negList[[j]]
            remainderVec <- X1 + adjVec * X2[j]
            limitVec <- remainderVec/adjVec
            maxX2 <- max(min(limitVec[posVec]), 0)
            minX2 <- max(max(limitVec[negVec]), 0)
            proposal <- sample(floor(maxX2)-ceiling(minX2)+1, 1)
            proposal <- proposal - 1 + ceiling(minX2)

            # Propose lambda_j+k | X2_j
            lambdaProp <- rgamma(1, proposal + prior$a[j+k] + 1,
                                 prior$b[j+k] + 1)

            # If feasible, MH step
            llr <- (log(lambdaProp)*proposal - lgamma(proposal+1) -
                    log(lambda[iter-1,j+k])*X2[j] + lgamma(X2[j]+1) )
            # Check feasibility (issues with integer constraints)
            llr <- llr + log(min(remainderVec - adjVec*proposal)>0)
            lir <- 0

            # Accept/reject part
            if (!is.na(llr-lir)) {
                if (log(runif(1)) < llr - lir) {
                    X2[j] <- proposal
                    accepts[j] <- accepts[j] + 1
                    X1 <- remainderVec - adjVec * X2[j]
                    lambda[iter,j+k] <- lambdaProp
                }
            }
        }

        X[iter,] <- c(X1,X2)

        # Draw lambda_1, ..., lambda_k
        lambda[iter,1:k] <- rgamma(k, X[iter,1:k] + prior$a[1:k] + 1,
                                   prior$b[1:k] + 1)

        # If verbose, print updates every 1e3 iterations
        if (verbose>0 & iter %% 1e3 == 0) {
            cat(sprintf('Iter %d\n', iter))
            if (verbose > 1) {
                print(lambda[iter,])
                print(X[iter,])
            }
        }
    }

    inv <- numeric(n)
    inv[pivot] <- seq(n)

    lambda <- lambda[,inv]
    X <- X[,inv]
    accepts <- accepts[inv]

    out <- tail(data.frame(X, lambda),ndraws)
    colnames(out) <- c(paste('X', seq(n), sep=''),
                       paste('lambda', seq(n), sep='') )

    return(list(lambdaDraws=mcmc(out[,seq(n+1,2*n)]), XDraws=mcmc(out[,seq(n)]),
                accepts=accepts))
}

