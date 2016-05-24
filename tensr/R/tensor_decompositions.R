#' Calculate the (truncated) higher-order SVD (HOSVD).
#'
#' Calculates the left singular vectors of each matrix unfolding of an array,
#' then calculates the core array. The resulting output is a Tucker
#' decomposition.
#'
#' If \code{r} is equal to the rank of \code{Y}, then \code{Y} is equal to
#' \code{atrans(S, U)}, up to numerical accuracy.
#'
#' More details on the HOSVD can be found in
#' \href{http://epubs.siam.org/doi/abs/10.1137/S0895479896305696}{ De Lathauwer
#' et. al. (2000)}.
#'
#' @param Y An array of numerics.
#' @param r A vector of integers. The rank of the truncated HOSVD.
#'
#' @return \code{U} A list of matrices with orthonormal columns. Each matrix
#'   contains the mode-specific singular vectors of its mode.
#'
#'   \code{S} An all-orthogonal array. This is the core array from the HOSVD.
#'
#' @export
#'
#' @author Peter Hoff.
#'
#' @references De Lathauwer, L., De Moor, B., & Vandewalle, J. (2000).
#'   \href{http://epubs.siam.org/doi/abs/10.1137/S0895479896305696}{A
#'   multilinear singular value decomposition}. \emph{SIAM journal on Matrix
#'   Analysis and Applications}, 21(4), 1253-1278.
#'
#' @keywords decompositions
#'
#' @examples
#' #Generate random data.
#' p <- c(2, 3, 4)
#' X <- array(stats::rnorm(prod(p)), dim = p)
#'
#' #Calculate HOSVD.
#' hosvd_x <- hosvd(X)
#' S <- hosvd_x$S
#' U <- hosvd_x$U
#'
#' #Recover X.
#' trim(X - atrans(S, U))
#'
#' #S is all-orthogonal.
#' trim(mat(S, 1) %*% t(mat(S, 1)))
#' trim(mat(S, 2) %*% t(mat(S, 2)))
#' trim(mat(S, 3) %*% t(mat(S, 3)))
hosvd <- function(Y, r = NULL) {
    ## higher order svd
    m <- dim(Y)
    K <- length(m)
    if (is.null(r)) {
        r <- pmin(m, prod(m) / m)
    }

    ## get standard tsvd
    U <- list()
    for (k in 1:K) {
        U[[k]] <- svd(mat(Y, k))$u[, 1:r[k], drop = FALSE]
    }
    S <- atrans(Y, lapply(U, t))


    ## now put into identifiable form
    if (all(r > 1)) {
        G <- sign(S)
        a <- list()
        for (k in 1:K) {
            a[[k]] <- rep(1, r[k])
            if (k == 1) {
                a[[k]][1] <- G[1]
            }
            idx <- matrix(1, nrow = r[k] - 1, ncol = K)
            idx[, k] <- 2:r[k]
            a[[k]][2:r[k]] <- a[[k]][1] * G[idx] / G[1]
        }
        aprod <- a[[1]]
        for (k in 2:K) {
            aprod <- outer(aprod, a[[k]])
        }
        U <- mapply(function(X, b) {
            X %*% diag(b)
        }, U, a, SIMPLIFY = FALSE)
        S <- S * aprod
    }

    list(U = U, S = S)
}

#' Calculate the incredible higher-order LQ decomposition (HOLQ).
#'
#' This function calculates a generalization of the LQ decomposition to tensors.
#' This decomposition has a close connection to likelihood inference in
#' Kronecker structured covariande models.
#'
#' Given an array \code{X}, the default version of this function will calculate
#' (1) \code{L} a list of lower triangular matricies with positive diagonal
#' elements and unit determinant, \code{Z} an array of the same dimensions as
#' \code{X} that has special orthogonal properties, and (3) \code{sig} a numeric
#' such that \code{X} is the same as \code{sig * atrans(Z,L)} up to numeric
#' precision.
#'
#' This output (1) can be considered a generalization of the LQ decomposition to
#' tensors, (2) solves an optimization problem which the matrix LQ decomposition
#' solves, and (3) has a special connection to likelihood inference in the array
#' normal model.
#'
#' There are options to constrain the matrices in \code{L} to either be
#' diagonal, lower triangular with unit diagonal, or the identity matrix. Each
#' of these correspond to submodels in Kronecker structured covariance models.
#' The core array corresponding to each of these options has different
#' properities (see \href{http://arxiv.org/abs/1410.1094}{Gerard and Hoff
#' (2014)}). These more constrained tensor decompositions are called HOLQ
#' juniors.
#'
#' The MLE of the \eqn{i}th component covariance matrix under \emph{any}
#' elliptically contoured Kronecker structured covariance model is given by
#' \code{L[[i]] \%*\% t(L[[i]])}. The MLE for the total variation pamarameter
#' will be different depending on the distribution of the array, but for the
#' array normal it is \code{sig ^ 2 / prod(p)} (the "variance" form for the
#' total variation parameter).
#'
#' The likelihood ratio test statistic depends only on \code{sig} and can be
#' implemented in \code{\link{lrt_stat}}.
#'
#' The algorithm used to fit the HOLQ iteratively repeats the LQ decomposition
#' along each mode.
#'
#' For more details on the incredible HOLQ, see
#' \href{http://arxiv.org/abs/1410.1094}{Gerard and Hoff (2014)}.
#'
#' @param X An array of numerics.
#' @param tol A numeric. The maximum difference in frobenius norm between two
#'   successive core arrays before stopping. Or maximum difference of the ratio
#'   of sigs from 1 before stopping (which depends on the value of
#'   \code{use_sig}).
#' @param itermax An integer. The maximum number of iterations of the LQ
#'   decomposition to do before stopping.
#' @param mode_rep A vector of integers. The optional mode(s) to be considered
#'   identity matrices.
#' @param mode_diag A vector of integers. The optional mode(s) to be considered
#'   as independent but heteroscedastic.
#' @param mode_ldu A vector of integers. The optional modes(s) to be considered
#'   as having unit diagonal.
#' @param print_diff A logical. Should the updates be printed to the screen each
#'   iteration?
#' @param start_vals Determines how to start the optimization algorithm. Either
#'   'identity' (default), or 'residuals', which results in using the cholesky
#'   square roots of the sample covariance matrices along each mode scaled to
#'   have unit determinant. You can also use your own start values.
#' @param use_sig A logical. What stopping criterion should we use? Frobenius
#'   norm of difference of cores (FALSE) or absolute value of difference of
#'   ratio of \code{sig} from 1 (TRUE).
#'
#' @return \code{Z} The core array with scaled all-orthonormality property.
#'
#'   \code{A} A list of matrices.
#'
#'   \code{sig} A numeric. The total variation parameter. This is the "standard
#'   deviation" form.
#'
#' @seealso \code{\link{array_bic_aic}} for using the output of \code{holq} to
#'   calculate AIC and BIC,
#'
#'   \code{\link{get_isvd}} for using the output of \code{holq} to calculate a
#'   tensor generalization of the singular value decomposition.
#'
#'   \code{\link{lq}} for the matrix LQ decomposition.
#'
#'   \code{\link{lrt_stat}} for using the output of \code{holq} to calculate
#'   likelihood ratio test statistics.
#'
#'   \code{\link{mle_from_holq}} for using the output of \code{holq} to
#'   calculate the maximum likelihood estimates of the component covariance
#'   matrices under the array normal model.
#'
#' @export
#'
#' @references Gerard, D. C., & Hoff, P. D. (2014).
#'   \href{http://arxiv.org/abs/1410.1094}{A higher-order LQ decomposition for
#'   separable covariance models}. \emph{arXiv preprint arXiv:1410.1094.}
#'
#' @author David Gerard.
#'
#' @keywords decompositions likelihood
#'
#' @examples
#' #Genrate random data.
#' p <- c(2, 3, 4)
#' X <- array(stats::rnorm(prod(p)), dim = p)
#'
#' #Calculate HOLQ with unit diagonal on 2nd mode,
#' #  and diagonal along 3rd mode.
#' holq_x <- holq(X, mode_ldu = 2, mode_diag = 3)
#' Z <- holq_x$Z
#' A <- holq_x$A
#' sig <- holq_x$sig
#'
#' #Reconstruct X
#' trim(X - sig * atrans(Z, A), 10^-5)
#'
#' #Properties of core
#' #First mode has orthonormal rows.
#' trim(mat(Z, 1) %*% t(mat(Z, 1)), 10^-5)
#'
#' #Second mode has orthogonal rows.
#' trim(mat(Z, 2) %*% t(mat(Z, 2)), 10^-5)
#'
#' #Third mode has unit diagonal (up to scale).
#' diag(mat(Z, 3) %*% t(mat(Z, 3)))
holq <- function(X, tol = 10 ^ -9, itermax = 1000, mode_rep = NULL, mode_diag = NULL, mode_ldu = NULL, print_diff = TRUE, start_vals = "identity",
    use_sig = TRUE) {

    p <- dim(X)  # the dimension of X
    n <- length(p)  # the number of modes (order)
    if (is.null(mode_rep)) {
        mode_rep <- n + 1
    }
    if (is.null(mode_diag)) {
        mode_diag <- n + 1
    }
    if (is.null(mode_ldu)) {
        mode_ldu <- n + 1
    }
    if (start_vals[1] == "identity") {
        ## start values are the identity
        A <- list()
        for (index in 1:n) {
            A[[index]] <- diag(p[index])
        }

        X_fnorm <- fnorm(X)

        Z_new <- X / X_fnorm
        sig <- X_fnorm
    } else if (start_vals[1] == "residuals") {
        ## start values are the empirical residuals
        resid_X <- start_resids(X, mode_rep)
        A <- resid_X$Sig
        Z_new <- atrans(X, resid_X$Sig_inv)
        Z_fnorm <- fnorm(Z_new)
        Z_new <- Z_new / Z_fnorm
        sig <- Z_fnorm
    } else {
        ## if supplied own starting values
        stopifnot(p == sapply(start_vals, nrow))
        A_det1 <- list()
        A_inv <- list()
        for (k in 1:n) {
            A_det1[[k]] <- start_vals[[k]] / det(start_vals[[k]]) ^ (1 / p[k])
            A_inv[[k]] <- backsolve(A_det1[[k]], diag(p[k]), upper.tri = FALSE)
        }
        Z_new <- atrans(X, A_inv)
        Z_fnorm <- fnorm(Z_new)
        Z_new <- Z_new / Z_fnorm
        sig <- Z_fnorm
        A <- A_det1
    }

    crit_vec <- c(tol + 1, tol + 1)

    iter_index <- 0
    while (crit_vec[use_sig + 1] > tol & iter_index <= itermax) {
        Z_old <- Z_new
        sig_old <- sig
        sig_old <- sig

        for (mat_index in (1:n)[-mode_rep]) {

            if (any(mat_index == mode_diag)) {
                ## assume unstructured along a mode multiply in sig first so that don't have to worry about it later
                Z_mat_temp <- sig * diag(A[[mat_index]]) * mat(Z_new, mat_index)
                L_k <- diag(sqrt(p[mat_index] * rowSums(Z_mat_temp ^ 2) / prod(p)))
                Z_mat_temp <- diag(1 / diag(L_k)) %*% Z_mat_temp
                Z_arr <- array(Z_mat_temp, dim = c(p[mat_index], p[-mat_index]))
                Z_new <- aperm(Z_arr, match(1:n, c(mat_index, (1:n)[-mat_index])))

                L_k_mult <- prod(diag(L_k)) ^ (1 / p[mat_index])
                L_k <- L_k / L_k_mult

                Z_new_fnorm <- fnorm(Z_new)
                Z_new <- Z_new / Z_new_fnorm

                sig <- L_k_mult * Z_new_fnorm

                A[[mat_index]] <- L_k
            } else if (any(mat_index == mode_ldu)) {
                ## assume unit lower triangular along a mode
                lq_z <- lq(mat(Z_new, mat_index))
                L_k_temp <- lq_z$L
                L_k_diag <- diag(L_k_temp)

                L_k <- L_k_temp %*% diag(1 / L_k_diag)

                Z_arr <- array(diag(L_k_diag) %*% t(lq_z$Q), dim = c(p[mat_index], p[-mat_index]))
                Z_new <- aperm(Z_arr, match(1:n, c(mat_index, (1:n)[-mat_index])))

                Z_new_fnorm <- fnorm(Z_new)
                Z_new <- Z_new / Z_new_fnorm

                sig <- sig * Z_new_fnorm

                A[[mat_index]] <- A[[mat_index]] %*% L_k
            } else {
                ## assume diagonal along a mode
                lq_z <- lq(mat(Z_new, mat_index))
                L_k <- lq_z$L

                Z_arr <- array(t(lq_z$Q), dim = c(p[mat_index], p[-mat_index]))
                Z_new <- aperm(Z_arr, match(1:n, c(mat_index, (1:n)[-mat_index])))

                L_k_mult <- prod(diag(L_k)) ^ (1 / p[mat_index])
                L_k <- L_k / L_k_mult

                Z_new_fnorm <- fnorm(Z_new)
                Z_new <- Z_new / Z_new_fnorm

                sig <- sig * L_k_mult * Z_new_fnorm

                A[[mat_index]] <- A[[mat_index]] %*% L_k
            }
        }
        z_diff <- fnorm(Z_old - Z_new)
        sig_diff <- abs(sig_old / sig - 1)
        crit_vec <- c(z_diff, sig_diff)
        if (print_diff) {
            cat("Scale Diff = ", sig_diff, "\n")
            cat("     Scale = ", sig, "\n\n")
        }
        iter_index <- iter_index + 1
    }
    return(list(Z = Z_new, A = A, sig = sig))
}

#' The incredible higher-order polar decomposition (IHOP).
#'
#' Mmm, pancakes.
#'
#' This function will calculate the higher-order polar decomposition, a
#' generalization of the polar decomposition to tensors. It generalizes a
#' minimization formulation of the polar decomposition.
#'
#' Given an array \code{X}, \code{ihop} will output \code{L} a list of lower
#' triangular matrices with positive diagonal elements and unit Frobenius norm,
#' \code{R} a core array with certain orthogonality properties, and \code{sig} a
#' total variation parameter. We have that \code{X} is equal to \code{sig *
#' atrans(R, L)} up to numerical precision.
#'
#' \code{t(solve(L[[i]])) \%*\% mat(R, i)} will have orthonormal rows for all
#' \code{i}.
#'
#' For more details on the IHOP, see
#' \href{http://arxiv.org/abs/1410.1094}{Gerard and Hoff (2014)}.
#'
#' @param X An array of numerics.
#' @param itermax An integer. The maximum number of iterations to perform during
#'   the optimization procedure.
#' @param tol A numeric. The algorithm will stop when the Frobenius norm of the
#'   difference of core arrays between subsequent iterations is below \code{tol}
#'   (for \code{use_sig = FALSE}) or when the absolute difference between the
#'   ratio of subsequent values of \code{sig} and 1 is less than \code{tol} (for
#'   \code{use_sig = TRUE}).
#' @param print_diff A logical. Should we print the updates of the algorithm?
#' @param mode_rep A vector. Which component matrices should be set to be the
#'   identity?
#' @param use_sig A logical. See \code{tol}.
#' @return \code{R} A core array which, in combination with \code{L}, has
#'   certain orthogonality properties.
#'
#'   \code{L} A list of lower triangular matrices with unit Frobenius norm.
#'
#'   \code{sig} A numeric.
#'
#' @export
#'
#' @keywords decompositions
#'
#' @references Gerard, D. C., & Hoff, P. D. (2014).
#'   \href{http://arxiv.org/abs/1410.1094}{A higher-order LQ decomposition for
#'   separable covariance models}. \emph{arXiv preprint arXiv:1410.1094.}
#'
#' @author David Gerard.
#'
#' @examples
#' #Generate random data.
#' p <- c(2, 3, 4)
#' X <- array(stats::rnorm(prod(p)), dim = p)
#'
#' #Calculate IHOP.
#' ihop_x <- ihop(X)
#' R <- ihop_x$R
#' L <- ihop_x$L
#' sig <- ihop_x$sig
#'
#' #Reconstruct X
#' trim(X - sig * atrans(R, L))
#'
#' #Orthogonality properties
#' ortho_1 <- t(solve(L[[1]])) %*% mat(R, 1)
#' trim(ortho_1 %*% t(ortho_1))
#'
#' ortho_2 <- t(solve(L[[2]])) %*% mat(R, 2)
#' trim(ortho_2 %*% t(ortho_2))
#'
#' ortho_3 <- t(solve(L[[3]])) %*% mat(R, 3)
#' trim(ortho_3 %*% t(ortho_3))
ihop <- function(X, itermax = 100, tol = 10 ^ -9, print_diff = TRUE, mode_rep = NULL, use_sig = TRUE) {
    p <- dim(X)
    n <- length(p)

    ### set starting values to identity
    A <- list()
    for (index in 1:n) {
        A[[index]] <- diag(p[index]) / sqrt(p[index])
    }
    X_fnorm <- fnorm(X)
    Z_new <- X / X_fnorm
    sig <- X_fnorm * prod(sqrt(p))

    iter_index <- 1

    if (is.null(mode_rep)) {
        mode_rep <- n + 1
    }

    crit_vec <- c(tol + 1, tol + 1)  ### use z_diff or sig_diff?

    while (iter_index <= itermax & crit_vec[use_sig + 1] > tol) {

        Z_old <- Z_new
        sig_old <- sig
        for (mode_index in 1:n) {
            if (any(mode_index == mode_rep)) {
                ## Do nothing
            } else {
                polar_z <- polar(A[[mode_index]] %*% mat(Z_new, mode_index))
                chol_p <- t(chol(polar_z$P))

                Z_arr <- array(t(chol_p) %*% polar_z$Z, dim = c(p[mode_index], p[-mode_index]))
                Z_new <- aperm(Z_arr, match(1:n, c(mode_index, (1:n)[-mode_index])))

                L_k_mult <- fnorm(chol_p)
                L_k <- chol_p / L_k_mult

                Z_new_fnorm <- fnorm(Z_new)
                Z_new <- Z_new / Z_new_fnorm
                sig <- sig * L_k_mult * Z_new_fnorm
                A[[mode_index]] <- L_k
                iter_index <- iter_index + 1
            }
        }
        z_diff <- fnorm(Z_old - Z_new)
        sig_diff <- abs(sig_old / sig - 1)
        crit_vec <- c(z_diff, sig_diff)
        if (print_diff) {
            cat("Scale Diff = ", sig_diff, "\n")
            cat("     Scale = ", sig, "\n\n")
        }
    }
    return(list(R = Z_new, L = A, sig = sig))
}

#' Calculate the incredible SVD (ISVD).
#'
#' The ISVD is a generalization of the SVD to tensors. It is derived from the
#' incredible HOLQ.
#'
#' Let \code{sig * atrans(Z, L)} be the HOLQ of \code{X}. Then the ISVD
#' calculates the SVD of each \code{L[[i]]}, call it \code{U[[i]] \%*\% D[[i]]
#' \%*\% t(W[[i]])}. It then returns \code{l = sig}, \code{U}, \code{D}, and
#' \code{V = atrans(Z, W)}. These values have the property that \code{X} is
#' equal to \code{l * atrans(atrans(V, D), U)}, up to numerical precision.
#' \code{V} is also scaled all-orthonormal.
#'
#' For more details on the ISVD, see
#' \href{http://arxiv.org/abs/1410.1094}{Gerard and Hoff (2014)}.
#'
#' @param x_holq The output from \code{\link{holq}}.
#'
#' @return l A numeric.
#' @return U A list of orthogonal matrices.
#' @return D A list of diagonal matrices with positive diagonal entries and unit
#'   determinant. The diagonal entries are in descending order.
#' @return V A scaled all-orthonormal array.
#'
#' @export
#'
#' @author David Gerard.
#'
#' @keywords decompositions
#'
#' @references Gerard, D. C., & Hoff, P. D. (2014).
#'   \href{http://arxiv.org/abs/1410.1094}{A higher-order LQ decomposition for
#'   separable covariance models}. \emph{arXiv preprint arXiv:1410.1094.}
#'
#' @examples
#' #Generate random data.
#' p <- c(4,4,4)
#' X <- array(stats::rnorm(prod(p)), dim = p)
#'
#' #Calculate HOLQ, then ISVD
#' holq_x <- holq(X)
#' isvd_x <- get_isvd(holq_x)
#' l <- isvd_x$l
#' U <- isvd_x$U
#' D <- isvd_x$D
#' V <- isvd_x$V
#'
#' #Recover X
#' trim(X - l * atrans(atrans(V, D), U))
#'
#' #V is scaled all-orthonormal
#' trim(mat(V, 1) %*% t(mat(V, 1)), epsilon = 10^-5)
#'
#' trim(mat(V, 2) %*% t(mat(V, 2)), epsilon = 10^-5)
#'
#' trim(mat(V, 3) %*% t(mat(V, 3)), epsilon = 10^-5)
get_isvd <- function(x_holq) {
    ## x_holq is the output of the holq
    p <- dim(x_holq$Z)
    n <- length(p)
    A_svd <- lapply(x_holq$A, svd)

    D <- vector(length = n, mode = "list")
    U <- vector(length = n, mode = "list")
    V_k <- vector(length = n, mode = "list")
    for (mode_index in 1:n) {
        U[[mode_index]] <- A_svd[[mode_index]]$u
        D[[mode_index]] <- diag(A_svd[[mode_index]]$d)
        V_k[[mode_index]] <- A_svd[[mode_index]]$v
    }
    V <- atrans(x_holq$Z, lapply(V_k, t))
    l <- x_holq$sig
    return(list(l = l, U = U, D = D, V = V))
}


#' Calculate the higher-order orthogonal iteration (HOOI).
#'
#' This function will calculate the best rank \code{r} (where \code{r} is a
#' vector) approximation (in terms of sum of squared differences) to a given
#' data array.
#'
#' Given an array \code{X}, this code will find a core array \code{G} and a list
#' of matrices with orthonormal columns \code{U} that minimizes \code{fnorm(X -
#' atrans(G, U))}. If \code{r} is equal to the dimension of \code{X}, then it
#' returns the HOSVD (see \code{\link{hosvd}}).
#'
#' For details on the HOOI see
#' \href{http://epubs.siam.org/doi/abs/10.1137/S0895479898346995}{Lathauwer et
#' al (2000)}.
#'
#' @param X An array of numerics.
#' @param r A vector of integers. This is the given low multilinear rank of the
#'   approximation.
#' @param tol A numeric. Stopping criterion.
#' @param print_fnorm Should updates of the optimization procedure be printed?
#'   This number should get larger during the optimizaton procedure.
#' @return \code{G} An all-orthogonal core array.
#'
#'   \code{U} A vector of matrices with orthonormal columns.
#'
#' @export
#'
#' @author David Gerard.
#'
#' @keywords decompositions
#'
#' @references De Lathauwer, L., De Moor, B., & Vandewalle, J. (2000).
#'   \href{http://epubs.siam.org/doi/abs/10.1137/S0895479898346995}{On the best
#'   rank-1 and rank-(\eqn{r_1, r_2,..., r_n}) approximation of higher-order tensors}.
#'   \emph{SIAM Journal on Matrix Analysis and Applications}, 21(4), 1324-1342.
#'
#' @examples
#' #Generate random data.
#' p <- c(2, 3, 4)
#' X <- array(stats::rnorm(prod(p)), dim = p)
#'
#' #Calculate HOOI
#' r <- c(2, 2, 2)
#' hooi_x <- hooi(X, r = r)
#' G <- hooi_x$G
#' U <- hooi_x$U
#'
#' #Reconstruct the hooi approximation.
#' X_approx <- atrans(G, U)
#' fnorm(X - X_approx)
hooi <- function(X, r, tol = 10 ^ -6, print_fnorm = FALSE) {
    p <- dim(X)
    n <- length(p)
    U_new <- hosvd(X, r = r)$U
    epsilon <- tol + 1  ## how far U_old is from U_new
    while (epsilon > tol) {
        U_old <- U_new
        ## do svd along each mode
        for (mode_index in 1:n) {
            U_mult <- U_new
            U_mult[[mode_index]] <- diag(p[mode_index])
            Y <- atrans(X, lapply(U_mult, t))
            U_new[[mode_index]] <- svd(mat(Y, mode_index))$u[, 1:r[mode_index]]
        }

        ## calculate how far the matrices are from one another
        dist_current <- rep(NA, length = n)
        for (mode_index in 1:n) {
            dist_current[mode_index] <- fnorm(U_new[[mode_index]] - U_old[[mode_index]])
        }
        epsilon <- max(dist_current)

        if (print_fnorm == TRUE) {
            cat(fnorm(atrans(X, lapply(U_new, t))), "\n")
        }
    }
    G <- atrans(X, lapply(U_new, t))
    return(list(G = G, U = U_new))
}
