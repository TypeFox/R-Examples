#' Generalized Procrustes Rotation
#' 
#' This function finds K orthogonal rotation matrices so that the rotated versions of the input 
#' configurations match each other optimally in the least-squares sense. The algorithm depends on the
#' starting values for the rotation matrices. At present identity matrices are used as starting values.
#' Only rotations / reflections are considered -- no scaling or translation factors are included.
#' 
#' @param configs A list of original configuration matrices
#' @param maxit The maximum number of iterations allowed
#' @param reltol The relative error tolerance for determining numeric convergence.
#' @param random Logical indicating whether or not to use random starts (only applicable when the
#' dimensionality is two).
#' @references
#' Gower, J. C., & Dijksterhuis, G. B. (2004). Procrustes problems (Vol. 3). Oxford: Oxford University Press.
#' @export
genproc <- function(configs, maxit = 50L, reltol = 1e-6, random = FALSE) {
  
  ## Sanity checks
  if (!is.list(configs)) stop("Arguments 'configs' must be a list.")
  if (length(configs) == 1) stop("Only a single configuration given.")
  if (length(unique(sapply(configs, length))) != 1) stop("All configurations must be of the same dimension.")
  
  ## Preliminaries
  K <- length(configs)
  p <- ncol(configs[[1]])
  crit <- rep(NA, maxit)
  iter <- 0
  
  ## Initialize rotation matrices randomly or to identity
  if (p == 2 && random) {
    rot2 <- function(angle) matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), nrow = 2, ncol = 2)
    angles <- runif(K, min = 0, max = 2*pi)
    Qmats <- lapply(angles, rot2)
  } else {
    Qmats <- replicate(K, diag(p), simplify = FALSE)
  }
  
  ## Initialize current configurations
  curr <- Map("%*%", configs, Qmats)
  
  ## Procrustes algorithm
  repeat {
    
    ## Increase counter
    iter <- iter + 1
    
    for (k in seq_len(K)) {
      
      ## Get Gk
      Gk <- Reduce('+', curr[-k]) / (K - 1)
      
      ## Get svd and update rotation matrix
      svdGk <- svd(crossprod(Gk, curr[[k]]))
      Qcur <- tcrossprod(svdGk$v, svdGk$u)
      Qmats[[k]] <- Qmats[[k]] %*% Qcur
      
      ## Calculate rotation
      curr[[k]] <- curr[[k]] %*% Qcur
    }
    
    ## Calculate criterion and check convergence
    G <- Reduce('+', curr) / K
    crit[iter] <- K * sum(sapply(Map("-", curr, list(G)), function(x) sum(x^2)))
    if (iter > 1) {
      if (crit[iter] / crit[iter - 1] > 1 - reltol) {
        break
      }
    }
  }
  return(list(original = configs, rotated = curr, rotations = Qmats, iter = iter, 
              crit = crit[1:iter], mincrit = crit[iter]))
}