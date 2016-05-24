#' Trace of a matrix.
#'
#' Returns the sum of the diagonal elements of a matrix.
#'
#' This returns the trace of a matrix, which is just the sum of its
#' diagonal elements.
#'
#' @param X A matrix whose diagonal elements will be added together.
#'
#' @return The sum of the diagonal elements of X.
#'
#' @export
#'
#' @author Peter Hoff.
#'
#'
#' @examples
#' X <- matrix(1:4, nrow = 2, ncol = 2)
#' X
#' tr(X)
tr <- function(X) {
    return(sum(diag(X)))
}

#' Unfold a matrix.
#'
#' \code{mat} returns a matrix version of a provided tensor.
#'
#' Applies the matrix unfolding operator (also called 'matricization' or 'matrix
#' flattening' operator) on a provided tensor. There are multiple ways one could
#' do this. This function performs the matrix unfolding described in
#' \href{http://epubs.siam.org/doi/abs/10.1137/07070111X}{ Kolda and Bader
#' (2009)}.
#'
#' @param A An array to be unfolded.
#' @param k The mode, or dimension, along which the unfolding is to be applied.
#'
#' @return  A matrix  whose rows  index the  \eqn{k}th mode  and whose columns
#'   index every other mode.  The ordering of the columns is in lexicographical
#'   order of the indices of the array \eqn{A}.
#'
#' @author Peter Hoff.
#'
#' @export
#'
#' @references Kolda, T. G., & Bader, B. W. (2009).
#'   \href{http://epubs.siam.org/doi/abs/10.1137/07070111X}{Tensor
#'   decompositions and applications}. \emph{SIAM review}, 51(3), 455-500.
#'
#'
#' @examples
#' A <- array(1:8, dim = c(2,2,2))
#' mat(A, 1)
#' mat(A, 2)
#' mat(A, 3)
mat <- function(A, k) {
    Ak <- t(apply(A, k, "c"))
    if(nrow(Ak) != dim(A)[k]){
      Ak <- t(Ak)
    }
    return(Ak)
}

#' \eqn{k}-mode product.
#'
#' \code{amprod} returns the \eqn{k}-mode product of an array with a
#' matrix.
#'
#' The \eqn{k}-mode product of a tensor \eqn{A} with a matrix \eqn{M}
#' results in a tensor whose \eqn{k}-mode unfolding is \eqn{M} times
#' the \eqn{k}-mode unfolding of \eqn{A}. That is
#' \code{mat(amprod(A,M,k)) = M \%*\% mat(A,k)}.  More details of the
#' \eqn{k}-mode product can be found in
#' \href{http://epubs.siam.org/doi/abs/10.1137/07070111X}{ Kolda and
#' Bader (2009)}.
#'
#' @param A A real valued array.
#' @param M A real matrix.
#' @param k An integer. The mode along which \code{M} is to be
#'     multiplied to \code{A}.
#'
#' @return An array whose \eqn{k}-mode unfolding is \code{M \%*\%
#'     mat(A,k)}.
#'
#' @seealso \code{\link{atrans}} for applying multiple \eqn{k}-mode
#'     products.
#'
#' @author Peter Hoff.
#'
#'
#' @references Kolda, T. G., & Bader, B. W. (2009).
#'   \href{http://epubs.siam.org/doi/abs/10.1137/07070111X}{Tensor
#'   decompositions and applications}. \emph{SIAM review}, 51(3), 455-500.
#'
#' @export
#'
#' @examples
#' A <- array(1:8, dim = c(2,2,2))
#' M <- matrix(1:4, nrow = 2, ncol = 2)
#' Y <- amprod(A, M, 2)
#' Y
#' identical(M %*% mat(A,2), mat(Y,2))
amprod <- function(A, M, k) {
    K <- length(dim(A))
    AM <- M %*% mat(A, k)
    AMA <- array(AM, dim = c(dim(M)[1], dim(A)[-k]))
    return(aperm(AMA, match(1:K, c(k, (1:K)[-k]))))
}

#' Tucker product.
#'
#' Performs the Tucker product between an array and a list of matrices.
#'
#' The Tucker product between a list of matrices \code{B} and an array \code{A}
#' is formally equivalent to performing the \eqn{k}-mode product between
#' \code{A} and each list element in \code{B}. For example, if the dimension of
#' \code{A} is three, then \code{atrans(A,B) =
#' amprod(amprod(amprod(A,B[[1]],1),B[[2]],2),B[[3]],3)}.  The ordering of this
#' \eqn{k}-mode product does not matter. See
#' \href{http://epubs.siam.org/doi/abs/10.1137/07070111X}{ Kolda and Bader
#' (2009)} for details.
#'
#' @param A An array of dimension \eqn{K}.
#' @param B A list of matrices of length \eqn{K}. It must be that
#'   \code{ncol(B[[k]]) == dim(A)[k]}.
#'
#' @seealso \code{\link{amprod}} for multiplying one matrix along one mode of an
#'   array.
#'
#' @author Peter Hoff.
#'
#' @references Kolda, T. G., & Bader, B. W. (2009).
#'   \href{http://epubs.siam.org/doi/abs/10.1137/07070111X}{Tensor
#'   decompositions and applications}. \emph{SIAM review}, 51(3), 455-500.
#'
#' @export
#'
#' @examples
#' A <- array(1:8, dim = c(2,2,2))
#' B <- list()
#' B[[1]] <-matrix(1:4, nrow = 2)
#' B[[2]] <- matrix(1:6, nrow = 3)
#' B[[3]] <- matrix(1:2, nrow = 1)
#' atrans(A,B)
atrans <- function(A, B) {
    X <- A
    for (k in 1:length(B)) {
        X <- amprod(X, B[[k]], k)
    }
    return(X)
}


#' The symmetric square root of a positive definite matrix.
#'
#' Returns the unique symmetric positive definite square root matrix
#' of a provided symmetric positive definite matrix.
#'
#' @param M A symmetric positive definite matrix.
#'
#' @return The unique symmetric positive definite matrix \eqn{X} such
#'     that \eqn{XX = M}.
#'
#' @author Peter Hoff.
#'
#' @export
#'
#' @examples
#' Y <- matrix(stats::rnorm(4), nrow = 2)
#' M <- Y %*% t(Y)
#' X <- mhalf(M)
#' X
#' identical(M, X %*% X)
mhalf <- function(M) {
    ## symmetric square root of a pos def matrix
    tmp <- eigen(M)
    return(tmp$vec %*% sqrt(diag(tmp$val, nrow = nrow(M))) %*% t(tmp$vec))
}

#' Frobenius norm of an array.
#'
#' Calculates the Frobenius norm of an array.
#'
#' The Frobenius norm of an array is the square root of the sum of its
#' squared elements. This function works for vector and matrix
#' arguments as well.
#'
#' @param X An array, a matrix, or a vector.
#'
#' @return The square root of the sum of the squared elements of
#'     \code{X}.
#'
#' @export
#'
#' @author David Gerard.
#'
#' @examples
#' X <- c(1:8)
#' Y <- matrix(1:8, nrow = 2)
#' Z <- array(1:8, dim = c(2, 2, 2))
#' fnorm(X)
#' fnorm(Y)
#' fnorm(Z)
fnorm <- function(X) {
    return(sqrt(sum(X ^ 2)))
}

#' Log-likelihood of array normal model.
#'
#' \code{ldan} calculates the log-likelihood of the array normal
#' model, minus a constant.
#'
#' @export
#'
#' @author David Gerard.
#'
#' @param E An array. This is the data.
#' @param Sig A list of symmetric positive definite matrices. These
#'     are the component covariance matrices.
ldan <- function(E, Sig) {
    m <- dim(E)
    ld <- 0
    for (k in 1:length(dim(E))) {
        ld <- ld - 0.5 * prod(m[-k]) * log(det(Sig[[k]]))
    }
    Z <- atrans(E, lapply(lapply(Sig, mhalf), solve))
    return(ld - 0.5 * sum(Z ^ 2))
}

#' QR Decomposition.
#'
#' QR decomposition, constraining the R matrix to have non-negative diagonal
#' entries.
#'
#' This function is almost a wrapper for \code{qr()}, \code{qr.R()}, and
#' \code{qr.Q()}, except it constrains the diagonal elements of \code{R} to be
#' non-negative. If \code{X} is full rank with fewer columns than rows, then
#' this is sufficient to gaurantee uniqueness of the QR decomposition
#' (Proposition 5.2 of
#' \href{https://books.google.com/books?id=WyvvAAAAMAAJ}{Eaton (1983)}).
#'
#' @param X A matrix of dimension \eqn{n} by \eqn{p} where \eqn{n \ge p}
#'
#' @author David Gerard.
#'
#' @return \code{Q} An \eqn{n} by \eqn{p} matrix with orthonormal columns.
#'
#' \code{R} A \eqn{p} by \eqn{p} upper-triangular matrix with non-negative
#' diagonal elements.
#'
#' @seealso \code{\link[base]{qr}}, \code{\link[base]{qr.Q}}, and
#'   \code{\link[base]{qr.R}} for the base methods on the obtaining the QR
#'   decomposition. \code{\link{lq}} for the related LQ decomposition.
qr2 <- function(X) {
    # just re-jiggers qr decomposition to make diagonals all positive
    qr_X <- qr(X)
    R <- qr.R(qr_X)
    Q <- qr.Q(qr_X)

    # re-jigger to make diagonals all positive
    sign_vec <- sign(diag(R))

    Q <- t(t(Q) * sign_vec)  # multiply the columns of Q
    R <- R * sign_vec  # multiply the rows of R
    return(list(R = R, Q = Q))
}

#' Truncates small numbers to 0.
#'
#' Given an array, matrix, or vector, \code{trim} will truncate all
#' elements smaller than \code{epsilon} (in absolute value) to zero.
#'
#' All elements in \code{X} that are smaller than \code{epsilon} (in
#' absolute value) will be set to zero then returned.
#'
#' @param X An array, a matrix, or a vector.
#' @param epsilon A numeric.
#'
#' @author David Gerard.
#'
#' @export
#'
#' @examples
#' X <- c(0, 1, 10^-7, -1, -10^-7)
#' X
#' trim(X)
trim <- function(X, epsilon = 10 ^ -6) {
    ## places 0's whereever within epsilon of 0
    X[abs(X) < epsilon] <- 0
    return(X)
}

#' Element-wise matrix products between two lists.
#'
#' Given two lists of matrices with conformable dimensions,
#' \code{listprod} returns a list whose elements are the matrix
#' products of the elements of these two lists.
#'
#' @param A A list of matrices.
#' @param B A second list of matrices.
#'
#' @author David Gerard.
#'
#' @export
#'
#' @return A list \code{C} such that \code{C[[i]] = A[[i]] \%*\% B[[i]]}.
listprod <- function(A, B) {
    ## A and B are lists of matrices with appropriate dimensions
    C <- list()
    for (index in 1:length(A)) {
        C[[index]] <- A[[index]] %*% B[[index]]
    }
    return(C)
}

#' LQ decomposition.
#'
#' Computes the LQ decomposition of a matrix.
#'
#' If \eqn{X} is an \eqn{n} by \eqn{p} matrix with \eqn{n \le p}, then
#' \code{lq} computes the LQ decomposition of \eqn{X}. That is, \eqn{X
#' = LQ'} where \eqn{Q} is \eqn{p} by \eqn{n} with orthonormal columns
#' and \eqn{L} is \eqn{n} by \eqn{n} lower triangular with positive
#' diaognal entries.
#'
#' @param X A \eqn{n} by \eqn{p} matrix of rank \eqn{n}.
#'
#' @return \code{L} An \eqn{n} by \eqn{n} lower triangular matrix with
#'     positive diagonal entries.
#'
#' \code{Q} An \eqn{n} by \eqn{p} matrix with orthonormal columns.
#'
#' The returned values satisfy \code{X = L \%*\% t(Q)}, up to
#' numerical precision.
#'
#' @seealso \code{\link{qr2}} for the related QR decomposition.
#'
#' @export
#'
#' @author David Gerard.
#'
#' @examples
#' X <- matrix(stats::rnorm(12), nrow = 3)
#' lq_X <- lq(X)
#' L <- lq_X$L
#' Q <- lq_X$Q
#' L
#' Q
#' trim(t(Q) %*% Q)
#' trim(X - L%*%t(Q))
lq <- function(X) {
    # computes LQ decomposition of a matrix.
    qr_tx <- qr(t(X))
    r <- qr.R(qr_tx)
    q <- qr.Q(qr_tx)

    # re-jigger to make diagonals all positive
    sign_vec <- sign(diag(r))

    q <- t(t(q) * sign_vec)  # multiply the columns of Q
    r <- r * sign_vec  # multiply the rows of R

    # Note: X = l %*% t(q)
    return(list(L = t(r), Q = q))
}

#' Collapse multiple modes into one mode.
#'
#' Given an array \code{X} and a vector of integers \code{m},
#' \code{collapse_mode} returns an array of lower order where the
#' first mode indexes the modes indicated in \code{m}.
#'
#' Transforms an array into another array where the provided modes are
#' collapsed into one mode. The indexing along this new mode is in
#' lexicographical order of the indices of the collapsed modes. The
#' collapsed mode is the first mode unless \code{length(m) == 1}, then
#' \code{collapse_mode} simply returns \code{X}.
#'
#' @param X An array whose modes we are collapsing.
#' @param m A vector of integers giving the modes to collapse.
#'
#' @return If \eqn{X} is of order \eqn{K} and \code{length(m) = q},
#'     then returns an array \eqn{Y} of order \eqn{K - q + 1}, where
#'     the modes indicated in \code{m} are combined to be the first
#'     mode in \eqn{Y}.
#'
#' @export
#'
#' @author David Gerard.
#'
#' @examples
#' X <- array(rep(c(1, 2), 8), dim = c(2, 2, 2, 2))
#' X
#' #mode 1 is now mode 2, modes 2, 3, and 4 are combined to be mode 1.
#' collapse_mode(X, c(2, 3, 4))
#' collapse_mode(X, c(2, 4)) ## another example.
#' collapse_mode(X, 4) #returns X
collapse_mode <- function(X, m) {
    p <- dim(X)
    n <- length(p)
    if (any(n < m) | length(m) < 2) {
        return(X)
    }

    return(apply(X, (1:n)[-m], c))
}

#' The left polar decomposition.
#'
#' \code{polar} calculates the left polar decomposition of a matrix.
#'
#' \code{polar} Takes a matrix \eqn{X}, of dimensions \eqn{n} by
#' \eqn{p}, and returns two matrices \eqn{P} and \eqn{Z} such that
#' \eqn{X = PZ}. \eqn{P} is a symmetric positive definite matrix of
#' dimension \eqn{n} by \eqn{n} and \eqn{Z} is an \eqn{n} by \eqn{p}
#' matrix with orthonormal rows.
#'
#' @param X A matrix.
#'
#' @return \code{P} A \eqn{n} by \eqn{n} symmetric positive definite
#'     matrix.
#'
#' \code{Z} A \eqn{n} by \eqn{p} matrix with orthonormal rows.
#'
#' Note that \code{X == P \%*\% Z}, up to numerical precision.
#'
#' @export
#'
#' @author David Gerard.
#'
#' @examples
#' X <- matrix(1:6, nrow = 2)
#' polar_x <- polar(X)
#' P <- polar_x$P
#' Z <- polar_x$Z
#' P
#' Z
#' trim(Z %*% t(Z))
#' trim(X - P %*% Z)
polar <- function(X) {
    ## calculates left polar decomposition
    x_svd <- svd(X)
    P <- x_svd$u %*% diag(x_svd$d) %*% t(x_svd$u)
    Z <- x_svd$u %*% t(x_svd$v)
    return(list(P = P, Z = Z))
}
