#' T3Clusf: Tucker3 Fuzzy Cluster Analysis
#' 
#' This is an implementation of the T3Clusf algorithm of Rocci & Vichi (2005).
#' 
#' @param X Three-way data array, with no missing values.
#' @param Q Integer giving the number of dimensions required for mode B (variables).
#' This is the first mode of the array, excluding the mode clustered over (see \code{margin}).
#' @param R Integer giving the number of dimensions required for mode C (occasions). 
#' This is the second mode of the array, excluding the mode clustered over (see \code{margin}).
#' @param G Integer giving the number of clusters required.
#' @param margin Integer giving the margin of the array to cluster over. The remaining two
#' modes, in the original order, corresponds to \code{Q} and \code{R}.
#' @param alpha Numeric value giving the fuzziness parameter.
#' @param eps Small numeric value giving the empirical convergence threshold.
#' @param maxit Integer giving the maximum number of iterations allowed.
#' @param verbose Integer giving the number of iterations after which the loss valus is printed.
#' @param nr.starts Integer giving the number of random starts required.
#' @param parallel Logical indicating whether to parallelize over random starts if 
#' \code{nr.starts > 1}.
#' @param mc.cores Argument passed to \code{\link{mclapply}}.
#' @export
#' @references 
#' Rocci, R., & Vichi, M. (2005). \emph{Three-mode component analysis with crisp or fuzzy partition of units}. 
#' Psychometrika, 70(4), 715-736.
#' @examples 
#' data("dcars")
#' set.seed(13)
#' res <- T3Clusf(X = dcars, Q = 3, R = 2, G = 3, alpha = 2)
#' 
T3Clusf <- function(X, Q, R = Q, G = 2, margin = 3L, alpha = 1, eps = 1e-8, maxit = 100L,
                    verbose = 1, nr.starts = 1L, parallel = TRUE, 
                    mc.cores = detectCores() - 1) {
  
  ## Recurse if multiple starts needed
  if (nr.starts > 1L) {
    if (parallel) {
      out <- mclapply(seq_len(nr.starts), function(a) T3Clusf(X = X, Q = Q, R = R, G = G, 
                      margin = margin, alpha = alpha, eps = eps, maxit = maxit, 
                      verbose = verbose, nr.starts = 1L), mc.cores = mc.cores)
    } else {
      out <- replicate(nr.starts, T3Clusf(X = X, Q = Q, R = R, G = G, margin = margin, 
                                          alpha = alpha, eps = eps, maxit = maxit, 
                                          verbose = verbose, nr.starts = 1L))
    }
    return(out)
  }
  
  ## Call
  cll <- match.call()
  
  ## Permute array so that clustering is over FIRST dimension
  if (!is(X, "array")) stop("'X' must be an array.")
  if (length(dim(X)) != 3) stop("'X' must be a three-dimensional array.")
  if (!(margin %in% seq_len(3)))  stop("'margin' must be either 1, 2, or 3.")
  X <- aperm.default(X, perm = c(margin, seq_len(3)[-margin]))
  dims <- dim(X)
  Xmat <- matrix(X, nrow = dims[1], ncol = prod(dims[-1]))
  
  ## Parameter checks
  if (any(is.na(X))) stop("'X' contains missing values.")
  if (alpha < 1) stop("'alpha' must be larger than or equal to 1.")
  
  ## Starting values
  ## Orthogonal B
  B <- matrix(rnorm(dims[2] * Q), nrow = dims[2], ncol = Q)
  B <- qr(B)
  B <- qr.Q(B)
  ## Orthogonal C
  C <- matrix(rnorm(dims[3] * R), nrow = dims[3], ncol = R)
  C <- qr(C)
  C <- qr.Q(C)
  ## Membership matrix U
  if (alpha == 1) {
    U <- sample.int(n = G, size = dims[1], replace = TRUE)
    U <- diag(G)[U, ]
  } else {
    U <- matrix(runif(dims[1] * G), nrow = dims[1], ncol = G)
    U <- diag(1 / rowSums(U)) %*% U
  }
  ## Staring values of Xmns (KJ x G) and Y (G x QR)
  Xmns <- apply(U^alpha, 2, function(z) apply(z * X, 2:3, sum) / sum(z))
  Y <- crossprod(Xmns, kronecker(C, B))
  
  ## Monitor loss
  loss <- rep(NA, maxit)
  iter <- 0L
  
  ## Iterate
  while (iter < maxit) {
    
    iter <- iter + 1L
    
    ## Update U
    CBY <- tcrossprod(kronecker(C, B), Y)
    dists <- apply(Xmat, 1, function(x, y) colSums((x - CBY)^2))
    if (alpha == 1) {
      U <- diag(G)[apply(dists, 2, which.min), ]
    } else {
      dists <- dists^(-1/(alpha - 1))
      U <- t(dists / matrix(colSums(dists), nrow = G, ncol = dims[1], byrow = TRUE))
    }
    
    ## Update Xmns and Y
    Xmns <- apply(U^alpha, 2, function(z) apply(z * X, 2:3, sum) / sum(z))
    Y <- crossprod(Xmns, kronecker(C, B))

    ## Check for empty clusters
    
    ## Update B
    omega <- colSums(U^alpha)
    CC <- tcrossprod(C)
    bmat <- matrix(0, nrow = dims[2], ncol = dims[2])
    for (i in seq_len(G)) {
      matX <- matrix(Xmns[, i], nrow = dims[2], ncol = dims[3])
      bmat <- bmat + omega[i] * tcrossprod(matX %*% CC, matX)
    }
    B <- eigen(bmat)$vectors[, seq_len(Q)]
    Y <- crossprod(Xmns, kronecker(C, B))
    
    ## Update C
    BB <- tcrossprod(B)
    cmat <- matrix(0, nrow = dims[3], ncol = dims[3])
    for (i in seq_len(G)) {
      matX <- matrix(Xmns[, i], nrow = dims[2], ncol = dims[3])
      cmat <- cmat + omega[i] * crossprod(matX, BB) %*% matX
    }
    C <- eigen(cmat)$vectors[, seq_len(R)]
    Y <- crossprod(Xmns, kronecker(C, B))
    
    ## Calculate the loss
    CC <- tcrossprod(C)
    CCBB <- kronecker(CC, BB)
    target <- CCBB %*% Xmns
    losscomps <- matrix(nrow = dims[1], ncol = G)
    for (i in seq_len(dims[1])) {
      losscomps[i, ] <- colSums((matrix(Xmat[i, ], nrow = prod(dims[-1]), ncol = G) - target)^2)
    }
    indloss <- rowSums(U^alpha * losscomps)
    loss[iter] <- sum(indloss)
    
    ## Monitor convergence
    if (iter %% verbose == 0)
      cat(sprintf(paste0("%", nchar(as.character(maxit)), "d"), iter), 
          "| Loss =", sprintf("%5f", loss[iter]), "\n")
    
    ## Check convergence
    if (iter > 2) {
      if (loss[iter] > loss[iter - 1]) 
        warning("The loss increased in iteration ", iter)
      if (1 - loss[iter] / loss[iter - 1] < eps) break
    }
  }
  
  ## Reshape Y into array
  Yarr <- array(t(Y), dim = c(Q, R, G))
  
  out <- list(U = U, B = B, C = C, Y = Yarr, iter = iter, loss = loss[seq_len(iter)], 
              minloss = loss[iter], sizes = colSums(U), call = cll)
  return(out)
}