#' Kendall's tau measure of association.
#'
#' This function provides a Monte Carlo approximation to Kendall's tau
#' measure of association.
#'
#' @param x a vector.
#' @param y a vector.
#' @param nmc an integer number of Monte Carlo simulations.
#'
#' @return A Monte Carlo approximation to Kendall's tau measure of
#'     association.
#'
#' @author Peter Hoff.
#'
#' @examples
#' mu <- rexp(30)
#' tensr:::kendalltau(rpois(30, mu), rpois(30, mu))
#'
kendalltau <- function(x, y, nmc = 1e+05) {
    i <- sample(length(y), nmc, replace = TRUE)
    j <- sample(length(y), nmc, replace = TRUE)
    s1 <- cbind(x[i], y[i])
    s2 <- cbind(x[j], y[j])
    cp <- sum( (s1[, 1] - s2[, 1]) * (s1[, 2] - s2[, 2]) > 0, na.rm = TRUE)
    dp <- sum( (s1[, 1] - s2[, 1]) * (s1[, 2] - s2[, 2]) < 0, na.rm = TRUE)
    return( (cp - dp) / sum(cp + dp))
}





#' Commutation matrix.
#'
#' Construct the communtation matrix.
#'
#' This function constructs the commutation matrix \code{K}, which maps
#' \code{c(A)} to \code{c(t(A))} for an \eqn{m} by \eqn{n} matrix
#' \code{A}.
#'
#' @references Magnus, J. R., & Neudecker,
#'     H. (1979). \href{http://www.janmagnus.nl/papers/JRM005.pdf}{The
#'     commutation matrix: some properties and
#'     applications}. \emph{The Annals of Statistics}, 381-394.
#'
#' Tracy, D. S., & Dwyer,
#' P. S. (1969). \href{http://www.che.iitm.ac.in/~naras/ch544/maxmin_matrixderivatives.pdf}{Multivariate
#' maxima and minima with matrix derivatives}. \emph{Journal of the
#' American Statistical Association}, 64(328), 1576-1594.
#'
#' @param m a natural number.
#' @param n another natural number.
#'
#' @return \code{K} The \code{m * n} by \code{m * n} commutation
#'     matrix.
#'
#' @author Peter Hoff.
#' @export
#' @examples
#' m <- 5 ; n <- 4
#' A <- matrix(stats::rnorm(m * n), m, n)
#' Kom(5, 4) %*% c(A) - c(t(A))
Kom <- function(m, n) {
    K <- matrix(0, m * n, m * n)
    for (i in 1:m) {
        for (j in 1:n) {
            ei <- rep(0, m)
            ei[i] <- 1
            ej <- rep(0, n)
            ej[j] <- 1
            K <- K + kronecker(ei %*% t(ej), ej %*% t(ei))
        }
    }
    return(K)
}



#' Normal scores.
#'
#' This function applies a quantile-quantile transformation to the
#' data, resulting in a distribution that is approximately normal but
#' has the same ranks as the original data.
#'
#' @param y A vector.
#' @param ties.method The option \code{ties.method} in the \code{rank}
#'     function.
#'
#' @return A vector of the same length as \code{y}.
#'
#' @author Peter Hoff.
#'
#' @examples
#' y <- rexp(100)
#' z <- tensr:::zscores(y)
#' par(mfrow = c(1, 3))
#' hist(y)
#' hist(z)
#' plot(y,z)
zscores <- function(y, ties.method = "average") {
    z <- stats::qnorm(rank(y, na.last = "keep", ties.method = ties.method) / (sum(!is.na(y)) + 1))
    names(z) <- names(y)
    m <- dim(y)
    if (length(m) == 2) {
        z <- matrix(z, nrow = m[1], ncol = m[2])
        dimnames(z) <- dimnames(y)
    }
    if (length(m) >= 3) {
        z <- array(z, dim = m)
        dimnames(z) <- dimnames(y)
    }
    return(z)
}





#' Tucker sum.
#'
#' Computes the Tucker sum of an array and a list of matrices.
#'
#' @export
#'
#' @param X A real array.
#' @param A A list of real matrices.
#'
tsum <- function(X, A) {
    ### needs work
    m <- sapply(A, function(x) {
        dim(x)[1]
    })
    K <- length(m)
    XA <- array(0, dim = c(m, dim(X)[-(1:length(m))]))

    for (k in 1:K) {
        XA <- sweep(XA, c(k, K + 1), A[[k]] %*% apply(X, c(k, K + 1), sum), "+")
    }
    return(XA)
}

#' Wishart simulation.
#'
#' Simulate a Wishart-distributed random matrix.
#'
#' This function simulates a Wishart random matrix using Bartletts
#' decomposition, as described in Everson and Morris (2000).
#'
#' @param S0 a positive definite matrix.
#' @param nu a positive scalar.
#' @author Peter Hoff.
#' @keywords multivariate simulation
#' @examples
#' # simulate several matrices and compute the mean.
#' SS <- matrix(0, 5, 5)
#' for(s in 1:1000) { SS <- SS + tensr:::rwish(diag(5), 3) }
#' SS / s
rwish <- function(S0, nu = dim(as.matrix(S0))[1] + 1) {
    S0 <- as.matrix(S0)
    S0h <- eigen(S0, symmetric = TRUE)
    S0h <- S0h$vec %*% diag(sqrt(S0h$val), nrow = length(S0h$val)) %*% t(S0h$vec)

    p <- dim(S0)[1]
    T <- matrix(0, p, p)
    T[lower.tri(T)] <- stats::rnorm(p * (p - 1) / 2)
    diag(T) <- sqrt(stats::rgamma(p, (nu - (1:p) + 1) / 2, 1 / 2))
    return(S0h %*% T %*% t(T) %*% S0h)
}


#' Multivariate normal simulation.
#'
#' Simulate a multivariate normal random matrix.
#'
#' This function simulates multivariate normal random vectors.
#'
#' @param n number of mvnormal vectors to simulate.
#' @param mu mean vector.
#' @param Sigma covariance matrix.
#' @param Sigma.chol Cholesky decomposition of \code{Sigma}.
#' @author Peter Hoff.
#' @keywords multivariate simulation
#' @examples
#' # Simulate several matrices and compute the mean.
#' Y <- tensr:::rmvnorm(100, c(1, 2, 3), matrix(c(3, 0, 1, 0, 1, -1, 1, -1, 2), 3, 3))
#' colMeans(Y)
#' cov(Y)
rmvnorm <- function(n, mu, Sigma, Sigma.chol = chol(Sigma)) {
    E <- matrix(stats::rnorm(n * length(mu)), n, length(mu))
    X <- t(t(E %*% Sigma.chol) + c(mu))
    if (n == 1) {
        X <- c(X)
    }
    return(X)
}


#' Array indices.
#'
#' Generates indices corresponding to subarrays.
#'
#' This function generates a matrix corresponding to all combinations
#' of a list of indices, to be used in subsetting arrays.
#'
#' @param saidx either a vector of the dimensions of a potential
#'     array, or a list of the indices in the subarray.
#' @author Peter Hoff.
#' @export
#' @examples
#' # all indices of an array
#' arrIndices(c(4, 3, 2))
#' # indices of a subarray
#' arrIndices(list(c(1, 3), c(4, 5), c(2, 3, 6)))
arrIndices <- function(saidx) {
    if (is.list(saidx)) {
        m <- sapply(saidx, length)
        midx <- which(array(0, m) == 0, arr.ind = TRUE)
        for (k in 1:length(saidx)) {
            midx[, k] <- saidx[[k]][midx[, k]]
        }
    }
    if (!is.list(saidx)) {
        midx <- which(array(0, saidx) == 0, arr.ind = TRUE)
    }
    return(midx)
}

#' Array normal conditional distributions.
#'
#' Conditional mean and variance of a subarray.
#'
#' This function calculates the conditional mean and variance in the array
#' normal model. Let \eqn{Y} be array normal and let \eqn{Y_a} be a subarray of
#' \eqn{Y}. Then this function will calculate the conditional means and
#' variances of \eqn{Y_a}, conditional on every other element in \eqn{Y}.
#'
#' @param Y A real valued array.
#' @param M Mean of \code{Y}.
#' @param S List of mode-specific covariance matrices of \code{Y}.
#' @param saidx List of indices for indexing sub-array for which the conditional
#'   mean and variance should be computed. For example, \code{said_x = list(1:2,
#'   1:2, 1:2)} will compute the conditional means and variances for the \eqn{2}
#'   by \eqn{2} by \eqn{2} sub-array Y[1:2, 1:2, 1:2]. This is conditional on
#'   every other element in \code{Y}.
#' @author Peter Hoff.
#' @keywords multivariate
#'
#' @references Hoff, P. D. (2011).
#'   \href{http://arxiv.org/abs/1008.2169}{Separable covariance arrays via the
#'   Tucker product, with applications to multivariate relational data}.
#'   \emph{Bayesian Analysis}, 6(2), 179-196.
#' @export
#' @examples
#' p <- c(4, 4, 4)
#' Y <- array(stats::rnorm(prod(p)), dim = p)
#' saidx <- list(1:2, 1:2, 1:2)
#' true_cov <- tensr::start_ident(p)
#' true_mean <- array(0, dim = p)
#' cond_params <- anorm_cd(Y = Y, M = true_mean, S = true_cov, saidx = saidx)
#'
#' ## Since data are independent standard normals, conditional mean is 0 and
#' ##    conditional covariance matrices are identities.
#' cond_params$Mab
#' cond_params$Sab
anorm_cd <- function(Y, M, S, saidx) {
    Yab <- Y
    Sab <- S
    Mab <- M

    for (k in 1:length(dim(Y))) {
        if (length(saidx[[k]]) < dim(Y)[k]) {
            mab <- dim(Yab)
            idx <- arrIndices(mab)

            a <- saidx[[k]]
            b <- (1:mab[k])[-a]
            ma <- mab
            ma[k] <- length(a)
            mb <- mab
            mb[k] <- length(b)

            Ma <- array(Mab[is.element(idx[, k], a)], ma)
            Mb <- array(Mab[is.element(idx[, k], b)], mb)
            Yb <- array(Yab[is.element(idx[, k], b)], mb)

            iSab <- lapply(Sab, solve)
            Mab <- Ma + amprod(Yb - Mb, Sab[[k]][a, b] %*% solve(Sab[[k]][b, b, drop = FALSE]), k)
            Sab[[k]] <- solve(iSab[[k]][a, a])
            Yab <- array(Yab[is.element(idx[, k], a)], ma)
        }
    }
    return(list(Mab = Mab, Sab = Sab))
}



#' Standard normal array.
#'
#' Generate an array of iid standard normal variables.
#'
#' This functions generates an array of dimension \code{dim} filled
#' with iid standard normal variables.
#'
#' @param dim a vector of positive integers.
#' @author Peter Hoff.
#' @keywords simulation multivariate
#' @examples
#' tensr:::rsan(c(5,4,3))
rsan <- function(dim) {
    return(array(stats::rnorm(prod(dim)), dim))
}

#' Top K elements of a vector.
#'
#' Identify top K elements of a vector.
#'
#' This function returns the indices corresponding to the top elements
#' of a vector.
#'
#' @param x The vector.
#' @param K The number of indices to return.
#' @param ignoreties If \code{FALSE}, will return a vector of the
#'     indices whose elements are greater than or equal to the Kth
#'     largest element, resulting in a vector possibly of length
#'     greater than \code{K} in the case of ties.
#' @author Peter Hoff.
#' @examples
#' x <- c(3, 6, 2, 4, 1)
#' tensr:::topK(x, 3)
topK <- function(x, K = 1, ignoreties = TRUE) {
    iK <- which(x >= sort(x, decreasing = TRUE)[K])
    iK <- (iK[order(-x[iK])])
    if (ignoreties) {
        iK <- iK[1:K]
    }
    return(iK)
}
