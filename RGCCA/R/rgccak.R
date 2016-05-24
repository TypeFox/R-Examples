#' The function rgccak() is called by rgcca() and does not have to be used by the user. 
#' The function rgccak() computes the RGCCA block components, outer weight vectors, etc., 
#' for each block and each dimension. Depending on the dimensionality of each block \eqn{\mathbf{X}_j , j = 1, \ldots, J}, 
#' the primal (when \eqn{n > p_j}) or the dual (when \eqn{n < p_j}) algorithm is used (see Tenenhaus et al. 2013) 
#' @param A  A list that contains the \eqn{J} blocks of variables. Either the blocks (\eqn{\mathbf{X}_1, \mathbf{X}_2, \ldots, \mathbf{X}_J}) or the residual matrices (\eqn{\mathbf{X}_{h1}, \mathbf{X}_{h2}, \ldots, \mathbf{X}_{hJ}}).
#' @param C  A design matrix that describes the relationships between blocks. (Default: complete design).
#' @param tau A \eqn{1 \times J} vector that contains the values of the shrinkage parameters \eqn{\tau_j}, \eqn{ j=1, \ldots J}. (Default: \eqn{\tau_j = 1}, \eqn{ j=1, \ldots, J}).
#' If tau = "optimal" the shrinkage intensity paramaters are estimated using the Schafer and Strimmer (2005) 
#' analytical formula. 
#' @param scheme Either "horst", "factorial" or "centroid" (default: centroid).
#' @param scale  if scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
#' @param verbose  Will report progress while computing if verbose = TRUE (default: TRUE).
#' @param init The mode of initialization to use in the RGCCA algorithm. The alternatives are either by Singular Value Decompostion or random (default : "svd").
#' @param bias A logical value for either a biaised or unbiaised estimator of the var/cov.
#' @param tol Stopping value for convergence.
#' @return \item{Y}{A \eqn{n \times J} matrix of RGCCA outer components}
#' @return \item{Z}{A \eqn{n \times J} matrix of RGCCA inner components}
#' @return \item{a}{A list of outer weight vectors}
#' @return \item{crit}{The values of the objective function to be optimized in each iteration of the iterative procedure.}
#' @return \item{converg}{Speed of convergence of the algorithm to reach the tolerance.}
#' @return \item{AVE}{Indicators of model quality based on the Average Variance Explained (AVE): 
#' AVE(for one block), AVE(outer model), AVE(inner model).}
#' @return \item{C}{A design matrix that describes the relationships between blocks (user specified).}
#' @return \item{tau}{\eqn{1 \times J} vector containing the value for the tau penalties applied to each of the \eqn{J} blocks of data (user specified)}
#' @return \item{scheme}{The scheme chosen by the user (user specified).}
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Tenenhaus A. et al., (2013), Kernel Generalized Canonical Correlation Analysis, submitted.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @title Internal function for computing the RGCCA parameters (RGCCA block components, outer weight vectors, etc.).
#' @export rgccak

rgccak <- function (A, C, tau = "optimal", scheme = "centroid", scale = FALSE, 
                    verbose = FALSE, init="svd", bias = TRUE, tol = .Machine$double.eps) {
  A <- lapply(A, as.matrix)
  J <- length(A)
  n <- NROW(A[[1]])
  pjs <- sapply(A,NCOL)
  Y <- matrix(0, n, J)
  
  if (scale == TRUE) A <- lapply(A, function(x) scale2(x, bias = bias))
  if (!is.numeric(tau)) tau = sapply(A, tau.estimate)

  
  a <- alpha <- M <- Minv <- K <- list()
  which.primal <- which((n>=pjs)==1)
  which.dual <- which((n<pjs)==1)
    
    
  ####################################
  # Initialisation and normalisation #
  ####################################
  
  if (init=="svd") {
    #SVD Initialisation of a_j or \alpha_j
    for (j in which.primal){
      a[[j]] <- initsvd(A[[j]])
    }
    for (j in which.dual){
      alpha[[j]] <- initsvd(A[[j]])
      K[[j]] <- A[[j]]%*%t(A[[j]])
    }
  } else if (init=="random") {
    #Random Initialisation of a_j or \alpha_j
    for (j in which.primal) {
      a[[j]] <- rnorm(pjs[j])
    }
    for (j in which.dual) {
      alpha[[j]] <- rnorm(n)
      K[[j]] <- A[[j]]%*%t(A[[j]])
    }
  } else {
    stop("init should be either random or by SVD.")
  }
  
  N = ifelse(bias, n, n-1)  
    
  # Normalisation of a_j or \alpha_j
  for (j in which.primal) {
    ifelse( tau[j] == 1, {
      a[[j]] <- drop(1/sqrt(t(a[[j]]) %*% a[[j]])) * a[[j]]
      Y[, j] <- A[[j]] %*% a[[j]] } , {
      M[[j]] <- ginv(tau[j] * diag(pjs[j]) + ((1 - tau[j])/(N)) * (t(A[[j]]) %*% A[[j]]))
      a[[j]] <- drop(1/sqrt(t(a[[j]]) %*% M[[j]] %*% a[[j]])) * M[[j]] %*% a[[j]]
      Y[, j] <- A[[j]] %*% a[[j]]}
    )
  }
    
  for (j in which.dual) {
      ifelse(tau[j] == 1, {alpha[[j]] = drop(1/sqrt(t(alpha[[j]])%*%K[[j]]%*% alpha[[j]]))*alpha[[j]]
                           a[[j]] = t(A[[j]])%*%alpha[[j]]
                           Y[, j] = A[[j]] %*% a[[j]]}
                        , {M[[j]] = tau[j]*diag(n)+(1-tau[j])/(N)*K[[j]]
                          Minv[[j]] = ginv(M[[j]])
                          alpha[[j]] = drop(1/sqrt(t(alpha[[j]]) %*% M[[j]]%*%K[[j]]%*% alpha[[j]]))*alpha[[j]]
                          a[[j]] = t(A[[j]])%*%alpha[[j]]
                          Y[, j] = A[[j]] %*% a[[j]]}
            )
  }
    
  iter = 1
  crit = numeric()
  converg = numeric()
  Z = matrix(0, NROW(A[[1]]), J)
  a_temp = a
    
  g <- function(x) switch(scheme,horst=x,factorial=x**2,centroid=abs(x))
  
  repeat {
    Yold <- Y
            
    ################
    # Horst Scheme #
    ################
    if (scheme == "horst") {
                                
        for (j in which.primal) {
          ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * Y)
                             a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% t(A[[j]]) %*% Z[, j])) * t(A[[j]]) %*% Z[, j]
                             Y[, j] = A[[j]] %*% a[[j]]}
                          , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * Y)
                             a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * M[[j]] %*% t(A[[j]]) %*% Z[, j]
                             Y[, j] = A[[j]] %*% a[[j]]}
                )
        }
             
              
        for (j in which.dual) {
          ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * Y) 
                             alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j]))* Z[, j]
                             a[[j]] = t(A[[j]])%*% alpha[[j]]
                             Y[, j] = A[[j]] %*% a[[j]]}
          
                          , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * Y) 
                             alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j]))*(Minv[[j]] %*% Z[, j])
                             a[[j]] = t(A[[j]])%*% alpha[[j]]
                             Y[, j] = A[[j]] %*% a[[j]]}
                )
        }
    }

    ####################
    # Factorial Scheme #
    ####################
    
    if (scheme == "factorial") {
                                              
        for (j in which.primal){
          ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(cov(Y[, j], Y), n), n, J, byrow = TRUE) * Y)
                             a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% t(A[[j]]) %*% Z[, j])) * (t(A[[j]]) %*% Z[, j])
                             Y[, j] = A[[j]] %*% a[[j]]}
                          , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(cov(Y[, j], Y), n), n, J, byrow = TRUE) * Y)
                             a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * (M[[j]] %*% t(A[[j]]) %*% Z[, j])
                             Y[, j] = A[[j]] %*% a[[j]]}
                )
        }
          
        for (j in which.dual) {
          ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(cov(Y[, j], Y), n), n, J, byrow = TRUE) * Y)
                             alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j])) * Z[, j]
                             a[[j]] = t(A[[j]])%*% alpha[[j]]
                             Y[, j] = A[[j]] %*% a[[j]]}
                          , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(cov(Y[, j], Y), n), n, J, byrow = TRUE) * Y)
                            alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])) * (Minv[[j]] %*% Z[, j])
                            a[[j]] = t(A[[j]])%*% alpha[[j]]
                            Y[, j] = A[[j]] %*% a[[j]]}
                )
        }
    }

    ###################
    # Centroid Scheme #
    ###################
    
    if (scheme == "centroid") {
      for (j in which.primal) {
        ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * sign(matrix(rep(cov(Y[, j], Y), n), n, J, byrow = TRUE)) * Y)
                           a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% t(A[[j]]) %*% Z[, j])) * (t(A[[j]]) %*% Z[, j])
                           Y[, j] = A[[j]] %*% a[[j]]}
                        , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * sign(matrix(rep(cov(Y[, j], Y), n), n, J, byrow = TRUE)) * Y)
                           a[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * (M[[j]] %*% t(A[[j]]) %*% Z[, j])
                           Y[, j] = A[[j]] %*% a[[j]]}
              )
      }
      
      for (j in which.dual) {
        ifelse(tau[j]==1, {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * sign(matrix(rep(cov(Y[, j], Y), n), n, J, byrow = TRUE)) * Y)
                           alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j])) * Z[, j]
                           a[[j]] = t(A[[j]])%*% alpha[[j]]
                           Y[, j] = A[[j]] %*% a[[j]]}
                        , {Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * sign(matrix(rep(cov(Y[, j], Y), n), n, J, byrow = TRUE)) * Y)
                           alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])) * (Minv[[j]] %*% Z[, j])
                           a[[j]] = t(A[[j]])%*% alpha[[j]]
                           Y[, j] = A[[j]] %*% a[[j]]}
        )
      }
    }
  
    num_converg <- sum((rowSums(Yold) - rowSums(Y))^2)
    den_converg <- sum(rowSums(Yold)^2)
    converg[iter] <- num_converg/den_converg
    stationnary_point = rep(FALSE, length(A))
    for (j in 1:J) stationnary_point[j] = sum(round(abs(a_temp[[j]] - a[[j]]), 8) < tol) == NCOL(A[[j]])
  
    a_temp <- a
    crit[iter] <- sum(C*g(cov2(Y, bias = bias)))
  
    if (iter > 1000) warning("The RGCCA algorithm did not converge after 1000 iterations.")
    
    if ((converg[iter] < tol & sum(stationnary_point) == J | iter > 1000)) 
          break
    iter <- iter + 1
  }
  
  if(sum(stationnary_point) == J & verbose) cat("The RGCCA algorithm converged to a fixed point of the stationary equations after", iter-1, "iterations \n")
  if (verbose) plot(crit, xlab = "iteration", ylab = "criteria")
         

  AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)
  result <- list(Y = Y, a = a, crit = crit[which(crit != 0)], 
                converg = converg[which(converg != 0)],
                AVE_inner = AVEinner, C = C, tau = tau, scheme = scheme)
  
  return(result)
}
