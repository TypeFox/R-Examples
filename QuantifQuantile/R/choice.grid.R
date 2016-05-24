#' @export
#' @title Choice of the quantization grids
#' @name choice.grid
#' @description This function provides \code{ng} optimal quantization
#'  grids for \code{X}, with \code{N} fixed.
#' @param X vector or matrix that we want to quantize. 
#' @param N size of the quantization grids.
#' @param ng number of quantization grids needed.
#' @param p L_p norm optimal quantization.
#' @return An array of dimension \code{d}*\code{N}*\code{ng} that
#'  corresponds to \code{ng} quantization grids.
#' @seealso \code{\link{QuantifQuantile}}, \code{\link{QuantifQuantile.d2}} and 
#' \code{\link{QuantifQuantile.d}}
#' @details This function works for any dimension of \code{X}. If the covariate
#'  is univariate, \code{X} is a vector while \code{X} is a matrix with \code{d}
#'  rows when the covariate is \code{d}-dimensional.
#' @details These grids are constructed using a stochastic gradient algorithm,
#' called CLVQ when \code{p}=2.
#' @references Charlier, I. and Paindaveine, D. and Saracco, J.,
#' \emph{Conditional quantile estimator based on optimal 
#' quantization: from theory to practice}, Submitted.
#' @references Pages, G. (1998) \emph{A space quantization method for numerical 
#' integration}, Journal of Computational and Applied Mathematics, 89(1), 1-38
#' @examples
#' X <- runif(300,-2,2)
#' N <- 10
#' ng <- 20
#' choice.grid(X,N,ng)

choice.grid <- function(X, N, ng=1, p = 2) {
    if (!is.numeric(X)) 
        stop("X must be numeric")
    if (!is.vector(X) & !is.matrix(X)) 
        stop("X must be a matrix or a vector")
    if ((!(floor(N) == N)) | (N <= 0)) 
        stop("N must be entire and positive")
    if ((!(floor(ng) == ng)) | (ng <= 0)) 
        stop("B must be entire")
    if (p < 1) 
        stop("p must be at least 1")
    if (is.vector(X)) {
        d <- 1
        n <- length(X)
        primeX <- matrix(sample(X, n * ng, replace = TRUE), 
            nrow = ng)
        hatX <- replicate(ng, sample(unique(X), N, replace = FALSE))
        hatX0 <- hatX  #save the initial grids
        gammaX <- array(0, dim = c(ng, n + 1))
        # initialisation of the step parameter gamma
        a_gamma <- 4 * N
        b_gamma <- pi^2 * N^(-2)
        BtestX = array(Inf, dim = c(ng, 1))
        # choice of gamma0X
        projXbootinit <- array(0, dim = c(n, ng))
        # index of the grid on which X is projected
        iminx <- array(0, dim = c(n, ng))
        for (i in 1:n) {
            RepX <- matrix(rep(X[i], N * ng), ncol = ng, byrow = TRUE)
            Ax <- sqrt((RepX - hatX0)^2)
            iminx[i, ] <- apply(Ax, 2, which.min)
            mx <- matrix(c(iminx[i, ], c(1:ng)), nrow = ng)
            projXbootinit[i, ] <- hatX0[mx]
        }
        RepX <- matrix(rep(X, ng), ncol = ng)
        distortion <- apply((RepX - projXbootinit)^2, 2, sum)/n
        temp_gammaX <- which(distortion > 1)
        if (length(temp_gammaX) > 0) {
            distortion[temp_gammaX] <- array(1, dim = c(length(temp_gammaX), 
                1))
        }
        gamma0X <- distortion
        if(any(gamma0X < 0.005)){gamma0X[gamma0X<0.005] <- 1}
        gammaX = array(0, dim = c(ng, n + 1))
        for (i in 1:ng) {
            # calculation of the step parameter
            gammaX[i, ] <- gamma0X[i] * a_gamma/(a_gamma + gamma0X[i] * 
                b_gamma * c(1:(n + 1) - 1))
        }
        gammaX[, 1] <- gamma0X
        iminX <- array(0, dim = c(n, ng))  #index that will change 
        tildeX <- array(0, dim = c(N, ng))
        # updating of the grids, providing optimal grids
        for (i in 1:n) {
            for (j in 1:ng) {
                tildeX[, j] <- matrix(rep(primeX[j, i], N), nrow = N, 
                  byrow = TRUE)
            }
            Ax <- sqrt((tildeX - hatX)^2)
            # calculation of each distance to determine the point of
            # the grid the closer of the stimuli
            iminX[i, ] <- apply(Ax, 2, which.min)
            mX <- matrix(c(iminX[i, ], c(1:ng)), nrow = ng)
            if(sqrt(sum((hatX[mX] - primeX[, i])^2))==0){
              hatX[mX] <- hatX[mX] - gammaX[, i + 1] * (hatX[mX] - primeX[, i])
            }else{
              hatX[mX] <- hatX[mX] - gammaX[, i + 1] * (hatX[mX] - 
                  primeX[, i]) * (sqrt(sum((hatX[mX] - primeX[, 
                  i])^2)))^(p - 1)/sqrt(sum((hatX[mX] - primeX[,                                                                                                                          i])^2))
            }
        }
    } else {
        n <- ncol(X)
        d <- nrow(X)
        primeX <- array(X[, sample(c(1:n), n * ng, 
            replace = T)], dim = c(d, n, ng))
        hatX <- replicate(ng, unique(X)[, sample(c(1:n), 
            N, replace = F)])  # initial grids chosen randomly in the sample
        hatX <- array(hatX, dim = c(d, N, ng))
        hatX0 <- hatX  #save the initial grids
        # initialisation of the step parameter gamma
        a_gamma <- 4 * N^(1/d)
        b_gamma <- pi^2 * N^(-2/d)
        BtestX <- array(Inf, dim = c(1, ng))
        # choice of gamma0X
        for (i in 1:(N - 1)) {
            for (j in (i + 1):N) {
                Bx <- array(0, dim = c(ng, 1))
                Bx <- sqrt(apply((hatX[, i, , drop = FALSE] - 
                  hatX[, j, , drop = FALSE])^2, c(2, 3), sum))/2
                temp_gammaX <- which(Bx < BtestX)
                BtestX[temp_gammaX] <- Bx[temp_gammaX]
            }
        }
        temp_gammaX = which(BtestX > 1)
        if (length(temp_gammaX) > 0) {
            BtestX[temp_gammaX] <- array(1, dim = c(length(temp_gammaX), 
                1))
        }
        gamma0X <- BtestX
        if(any(gamma0X < 0.005)){gamma0X[gamma0X<0.005] <- 1}
        
        gammaX <- array(0, dim = c(ng, n + 1))
        for (i in 1:ng) {
            # calculation of the step parameter
            gammaX[i, ] <- gamma0X[i] * a_gamma/(a_gamma + gamma0X[i] * 
                b_gamma * c(1:(n + 1) - 1))
        }
        gammaX[, 1] <- gamma0X
        iminX <- array(0, dim = c(n, ng))  #index that will change 
        tildeX <- array(0, dim = c(d, N, ng))
        # updating of the grids, providing optimal grids
        for (i in 1:n) {
            for (j in 1:ng) {
                tildeX[, , j] <- matrix(rep(primeX[, i, j], N), 
                  nrow = N, byrow = FALSE)
            }
            Ax <- sqrt(apply((tildeX - hatX)^2, c(2, 3), sum))
            # calculation of each distance to determine the point of the grid 
            #the closer of the stimuli
            iminX[i, ] <- apply(Ax, 2, which.min)
            for (k in 1:d) {
                m <- matrix(c(rep(k, ng), iminX[i, ], c(1:ng)), ncol = 3)
                if(p==2){
                  hatX[m] = hatX[m] - gammaX[, i + 1] * (hatX[m] - primeX[k, i, ]) 
                }else{
                  hatX[m] = hatX[m] - gammaX[,i + 1] * (hatX[m] - primeX[k, i, ]) * 
                    (sqrt(sum((hatX[m] - primeX[k, i, ])^2)))^(p - 1)/sqrt(sum((hatX[m] - primeX[k, 
                    i, ])^2))
                } 
            }
        }
    }
    output <- list(init_grid=hatX0,opti_grid=hatX)
    output
} 
