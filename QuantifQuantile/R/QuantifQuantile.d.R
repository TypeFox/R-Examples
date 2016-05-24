#' @export
#' @title QuantifQuantile for general X
#' @name QuantifQuantile.d
#' @description Estimation of conditional quantiles using optimal quantization 
#' when \code{X} is d-dimensional.
#' 
#' @details \itemize{\item This function calculates estimated conditional 
#' quantiles with a method based on optimal quantization for any dimension for 
#' the covariate. The matrix of covariate \code{X} must have \code{d} rows 
#' (dimension). For particular cases of \code{d} =1 or 2, it is strongly 
#' recommended to use \code{\link{QuantifQuantile}} and 
#' \code{\link{QuantifQuantile.d2}} respectively (computationally faster). The 
#' argument \code{x} must also have \code{d} rows.
#' \item The criterion for selecting the number of quantizers is implemented in 
#' this function. The user has to choose a grid \code{testN} of possible values 
#' in which \code{N} will be selected. It actually minimizes some bootstrap 
#' estimated version of the ISE (Integrated Squared Error). More precisely, for 
#' \code{N} fixed, it calculates the sum according to \code{alpha} of 
#' \code{hatISE_N} and then minimizes the resulting vector to get \code{N_opt}.
#'  However, the user can choose to select a different value of \code{N_opt} for
#'  each \code{alpha} by setting \code{same_N=FALSE}. In this case, the vector 
#'  \code{N_opt} is obtained by minimizing each column of \code{hatISE_N} 
#'  separately. The reason why \code{same_N=TRUE} by default is that taking 
#'  \code{N_opt} according to \code{alpha} could provide crossing conditional 
#'  quantile curves (rarely observed for not too close values of \code{alpha}). 
#'  The function \code{\link{plot.QuantifQuantile}} 
#'  illustrates the selection of \code{N_opt}. If the graph is not decreasing 
#'  then increasing, the argument \code{testN} should be adapted.
#'  \item This function can use parallel computation to save time, by simply 
#'  increasing the parameter \code{ncores}. Parallel computation relies on 
#'  \code{\link{mclapply}} from \code{\link{parallel}} package, hence is not available
#'  on Windows unless \code{ncores}=1 (default value).
#'  }
#'  
#' @param X matrix of covariates.
#' @param Y vector of response variables.
#' @param x matrix of values for \code{x} in q_alpha(x).
#' @param alpha vector of order of the quantiles.
#' @param testN grid of values of \code{N} that will be tested.
#' @param p L_p norm optimal quantization.
#' @param B number of bootstrap replications for the bootstrap estimator.
#' @param tildeB number of bootstrap replications for the choice of \code{N}.
#' @param same_N whether to use the same value of \code{N} for each \code{alpha}
#' (\code{TRUE} by default).
#' @param ncores number of cores to use. Default is set to 1 (see Details below).
#' 
#' @return An object of class \code{QuantifQuantile} which is a list with the 
#' following components:
#' @return \item{hatq_opt}{A matrix containing the estimated conditional 
#' quantiles. The number of columns is the number of considered values for \code{x} 
#' and the number of rows the size of the order vector \code{alpha}. This object 
#'  can also be returned using the usual \code{fitted.values} function.}
#' @return \item{N_opt}{Optimal selected value for \code{N}. An integer if 
#' \code{same_N=TRUE} and a vector of integers of length \code{length(alpha)} 
#' otherwise.}
#' @return \item{hatISE_N}{The matrix of estimated ISE provided by our selection 
#' criterion for \code{N} before taking the mean according to \code{alpha}. The 
#' number of columns is then \code{length(testN)} and the number of rows 
#' \code{length(alpha)}.}
#' @return \item{hatq_N}{A 3-dimensional array containing the estimated 
#' conditional quantiles for each considered value for \code{alpha}, \code{x} and \code{N}.}
#' @return \item{X}{The matrix of covariates.}
#' @return \item{Y}{The vector of response variables.}
#' @return \item{x}{The considered vector of values for \code{x} in q_alpha(x).}
#' @return \item{alpha}{The considered vector of order for the quantiles.}
#' @return \item{testN}{The considered grid of values for \code{N} that were tested.}
#'  
#' @references Charlier, I. and Paindaveine, D. and Saracco, J.,
#' \emph{Conditional quantile estimation through optimal quantization}, 
#' Journal of Statistical Planning and Inference, 2015 (156), 14-30.
#' @references Charlier, I. and Paindaveine, D. and Saracco, J.,
#' \emph{Conditional quantile estimator based on optimal 
#' quantization: from theory to practice}, Submitted.
#' 
#' @seealso \code{\link{QuantifQuantile}} and \code{\link{QuantifQuantile.d2}} 
#' for particular dimensions one and two.
#' @seealso \code{\link{plot.QuantifQuantile}}, 
#' \code{\link{print.QuantifQuantile}}, \code{\link{summary.QuantifQuantile}}
#'
#' @examples
#' \dontrun{
#' set.seed(644925)
#' n <- 500
#' X <- runif(n,-2,2)
#' Y <- X^2+rnorm(n)
#' x <- seq(min(X),max(X),length=100)
#' res <- QuantifQuantile.d(X,Y,x,testN=seq(15,35,by=5))
#' }
#' \dontrun{
#' set.seed(272422)
#' n <- 1000
#' X <- matrix(runif(n*2,-2,2),ncol=n)
#' Y <- apply(X^2,2,sum)+rnorm(n)
#' x1 <- seq(min(X[1,]),max(X[1,]),length=20)
#' x2 <- seq(min(X[2,]),max(X[2,]),length=20)
#' x <- matrix(c(rep(x1,20),sort(rep(x2,20))),nrow=nrow(X),byrow=TRUE)
#' res <- QuantifQuantile.d(X,Y,x,testN=seq(90,140,by=10),B=20,tildeB=15)
#' }
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
#' @importFrom stats quantile


QuantifQuantile.d <- function(X, Y, x, alpha = c(0.05, 0.25, 
    0.5, 0.75, 0.95), testN = c(35, 40, 45, 50, 55), p = 2, B = 50, 
    tildeB = 20, same_N=TRUE,ncores=1) {
    if (!is.numeric(X)) 
        stop("Y must be numeric")
    if (!is.numeric(x)) 
        stop("x must be numeric")
    if (!is.vector(Y)) 
        stop("Y must be a vector")
    if (!all(floor(testN) == testN & testN > 0)) 
        stop("testN must have entire positive entries")
    if (!all(alpha > 0 & alpha < 1)) 
        stop("alpha must be strictly between 0 and 1")
    if ((!(floor(B) == B)) | (B <= 0)) 
        stop("B must be a positive entire")
    if ((!(floor(tildeB) == tildeB)) | (tildeB <= 0)) 
        stop("tildeB must be a positive entire")
    if (p < 1) 
        stop("p must be at least 1")
    if (!is.logical(same_N))
      stop("same_N must be logical")
    
    if (is.vector(X)) {
        d <- 1
        n <- length(X)
        X <- matrix(X, nrow = d)
        if (!is.vector(x)) 
            stop("x must have same dimension as X")
        x <- matrix(x, nrow = d)
    } else {
        if (!is.matrix(X)) 
            stop("X must be a matrix with d rows")
        n <- ncol(X)
        d <- nrow(X)
        if (!is.matrix(x)) 
            stop("x must be a matrix with d rows")
        if (nrow(x) != d) 
            stop("x must be a matrix with d rows")
    }
    
    m <- length(x)/d
    hatISE_N <- array(0, dim = c(length(alpha), length(testN)))
    hatq_N <- array(0, dim = c(length(alpha), m, length(testN)))
    
    primeX <- array(X[, sample(c(1:n), n * (B + tildeB), replace = T)], 
        dim = c(d, n, (B + tildeB)))
    
    #estimation for different values of N
    calc_hatq_N <- function(N){  
        hatX <- choice.grid(X, N, ng = (B + tildeB))$opti_grid
        # projection of the sample X on the B+tildeB optimal grids
        projXboot <- array(0, dim = c(d, n, B + tildeB))
        # index of the grid on which X is projected
        iminx <- array(0, dim = c(n, B + tildeB))
        for (i in 1:n) {
            RepX <- array(rep(X[, i], N * (B + tildeB)), dim = c(d, 
                N, B + tildeB))
            Ax <- sqrt(apply((RepX - hatX)^2, c(2, 3), sum))
            iminx[i, ] <- apply(Ax, 2, which.min)
            mx <- array(0, dim = c(d, B + tildeB, 3))
            for (k in 1:d) {
                mx[k, , ] <- matrix(c(rep(k, B + tildeB), iminx[i, 
                  ], c(1:(B + tildeB))), ncol = 3)
                projXboot[k, i, ] <- hatX[mx[k, , ]]
            }
        }
        
        # estimation of q_alpha(x) for N fixed 
        # save the B+tildeB estimation of q_alpha(x)
        Hatq <- array(0, dim = c(m, length(alpha), B + tildeB))
        # save by Voronoi cell
        Hatq_cell <- array(0, dim = c(N, length(alpha), (B + 
            tildeB)))
        
        proj_gridx_boot <- function(z) {
            proj <- array(0, dim = c(d, B + tildeB))
            Repz <- array(rep(z, N * (B + tildeB)), dim = c(d, 
                N, B + tildeB))
            A <- sqrt(apply((Repz - hatX)^2, c(2, 3), sum))
            i_min <- apply(A, 2, which.min)
            mx <- array(0, dim = c(d, B + tildeB, 3))
            for (k in 1:d) {
                mx[k, , ] <- matrix(c(rep(k, B + tildeB), i_min, 
                  c(1:(B + tildeB))), ncol = 3)
                proj[k, ] <- hatX[mx[k, , ]]
            }
            proj
        }
        
        projection_x <- apply(x, FUN = proj_gridx_boot, MARGIN = 2)
        projection_x <- array(projection_x, dim = c(d, B + tildeB, 
            m))
        
        repY <- matrix(rep(Y, B + tildeB), nrow = n)
        
        # calculation of the conditional quantile for each cell. 
        # Since any point of a cell is projected on the center of this cell, 
        # the corresponding conditional quantiles are equal
        
        for (i in 1:N) {
            for (j in 1:(B + tildeB)) {
                a <- apply(projXboot[, , j, drop = F], FUN = function(z) {
                  all(z == hatX[, i, j])
                }, MARGIN = 2)
                a <- c(1:n)[a]
                if (length(a) > 0) {
                  Hatq_cell[i, , j] <- quantile(repY[a, j], probs = alpha)
                }
            }
        }
        
        # we now identify the cell in which belongs each x to associate the
        # corresponding value of conditional quantiles
        
        identification <- function(z) {
            identification <- array(0, dim = c(B + tildeB, 1))
            i <- apply(x, FUN = function(r) {
                all(r == z)
            }, MARGIN = 2)
            i <- c(1:m)[i]
            for (j in 1:(B + tildeB)) {
                identif <- apply(hatX[, , j, drop = F], FUN = function(z) {
                  all(z == projection_x[, j, i])
                }, MARGIN = 2)
                identification[j, ] <- c(1:N)[identif]
            }
            identification
        }
        
        identification_projection_x <- apply(x, FUN = identification, 
            2)
        
        for (i in 1:length(alpha)) {
            for (j in 1:m) {
                r <- matrix(c(identification_projection_x[, j], 
                  rep(i, B + tildeB), c(1:(B + tildeB))), ncol = 3)
                Hatq[j, i, ] <- Hatq_cell[r]
            }
        }
        
        # the final estimation is the mean of the B estimations
        hatq <- array(0, dim = c(m, length(alpha)))
        hatq <- apply(Hatq[, , c(1:B), drop = FALSE], c(1, 2), 
            mean)
        
        # the last tilde B are used to estimate the ISE
        HATq <- array(rep(hatq, tildeB), dim = c(m, length(alpha), 
            tildeB))
        hatISE <- (HATq - Hatq[, , c((1 + B):(B + tildeB)), drop = FALSE])^2
        hatISE <- apply(hatISE, 2, sum)/(m * tildeB)
        
        print(N)
        list(hatq=hatq,hatISE=hatISE)
    }
    
    parallel_hatq_hatISE <- mclapply(testN,calc_hatq_N,mc.cores = ncores,mc.set.seed=F)
    
    for(i in 1:length(testN)){
      hatq_N[,,i] <- t(parallel_hatq_hatISE[[i]]$hatq)
      hatISE_N[,i] <- parallel_hatq_hatISE[[i]]$hatISE
    }
    
    if(same_N){
      #choice of optimal N
      hatISEmean_N <- apply(hatISE_N, 2, mean)
      i_opt <- which.min(hatISEmean_N)
      #optimal value for N chosen as minimizing the sum of hatISE for the 
      # different alpha's
      N_opt <- testN[i_opt]
      
      # table of the associated estimated conditional quantiles
      hatq_opt <- hatq_N[, , i_opt, drop = F]
      hatq_opt <- matrix(hatq_opt, nrow = length(alpha))
    }else{
      #choice of optimal N
      i_opt <- apply(hatISE_N, 1, which.min)
      #optimal value for N chosen as minimizing the sum of hatISE for the 
      # different alpha's
      N_opt <- testN[i_opt]
      # table of the associated estimated conditional quantiles
      hatq_opt <- array(0, dim = c(length(alpha), length(x)/d))
      for(i in 1:length(alpha)){
        hatq_opt[i, ] <- hatq_N[i, , i_opt[i]]
      }
    }
    
    if(length(N_opt)==1){
      if(N_opt==min(testN)){warning("N_opt is on the left boundary of testN")}
      if(N_opt==max(testN)){warning("N_opt is on the right boundary of testN")}
    }else{
      if(any(N_opt==min(testN))){warning(
        "N_opt is on the left boundary of testN for at least one value of alpha")}
      if(any(N_opt==max(testN))){warning(
        "N_opt is on the right boundary of testN for at least one value of alpha")}
    }
    
    output <- list(hatq_opt = hatq_opt,fitted.values = hatq_opt, N_opt = N_opt, 
        hatISE_N = hatISE_N, hatq_N = hatq_N, X = X, Y = Y, x = x, 
        alpha = alpha, testN = testN)
    class(output) <- "QuantifQuantile"
    output
    ############################ 
    
}  # fin de la fonction 
