#' Reconciliate individual predictions using GTOP 
#'
#' Uses a Game Theory approach to reconciliate hierarchical 
#' time series predicitons 
#'   
#' In hierarchical time series forecasts, one predicts individuals
#' quantities and a global quantity. There exists a contraint that
#' matches the sum of the individual quantities to the global quantity.
#' However, forecasting models don't take into account this constraint.
#' With GTOP you can reconciliate the individual and global quantities 
#' in order to match the aggregate consistency contraint.
#'
#'       
#' @param preds_indiv vector contains the individual predictions
#' @param pred_total  prediction for the sum of individuals
#' @param weights_indiv  vector, contains the weights of the individuals
#' @param weight_total   weight of the total
#' @param bounds_indiv vector, contains the bounds of the individuals
#' @param solver string, use quadratic programming (\code{quad}) or
#'               Lasso-like solvers (\code{lasso}) 
#' @return A list with \itemize{
#'   \item{pred_indivs} the reconciliated predictions for the 
#'         individuals and the total, 
#'   \item{solution} the solution to the associate minimisation 
#'         problem.
#'         }
#' @examples
#' K <- 5
#' indiv <- rep(0, K)
#' total <- 1
#' gtop(preds_indiv   = indiv, 
#'      pred_total = total, 
#'      weights_indiv = rep(1, K), 
#'      weight_total = 2,
#'      bounds_indiv  = rep(1 / K, K))


gtop <- function(preds_indiv, pred_total,
                 weights_indiv, weight_total, 
                 bounds_indiv, solver = "quad") {
  if(solver == "quad") 
    res <- gtopQuad(preds_indiv,   pred_total,
                    weights_indiv, weight_total, 
                    bounds_indiv)
  else if(solver == "lasso") 
    res <- gtopLasso(preds_indiv,   pred_total,
                    weights_indiv, weight_total, 
                    bounds_indiv)
  else stop("Invalid option for solver")
  return(res)
} 



gtopQuad <- function(preds_indiv,   pred_total,
                     weights_indiv, weight_total, 
                     bounds_indiv) {
  
  if(any(weights_indiv <= 0)) stop("weights_indiv must all be positive")
  
  K <- length(preds_indiv)
  z <- pred_total - sum(preds_indiv)
  y <- c(preds_indiv, pred_total)
  
  A2 <- diag(c(weights_indiv, weight_total)) # gtop's A^2 matrix & 
  # solve.QP's Dmat matrix
  Rinv <- diag(1 / sqrt(c(weights_indiv, weight_total))) # A2 = R^T %*% R 
  
  
  # solve.QP's Dmat matrix
  
  Amat <- matrix(0, K + 1, 2 * K + 1) # solve.QP's A matrix (constraints matrix)
  Amat[, 1] <- c(rep(1, K), -1)                    # sum of indivs == total
  Amat[cbind(1:K, 1 + seq(1, 2 * K, by = 2))] <-  1    # lower boundaries 
  Amat[cbind(1:K, 1 + seq(2, 2 * K, by = 2))] <- -1    # upper boundaries
  
  # bvec : solve.QP's independent coefficients on the contraint
  # dvec : linear part of the objective function  
  bvec <- c(0, t(Amat[, -1]) %*% y - rep(bounds_indiv, each = 2))
  dvec <- A2 %*% y
  
  res <- solve.QP(Rinv, dvec, Amat, bvec, meq = 1, factorized = TRUE) 
  
  v <- - (t(res$solution - y) %*% A2 %*% (res$solution - y))
  
  return(list(preds_indiv = res$solution, solution = v))
} 



gtopLasso <- function(preds_indiv,   pred_total,
                      weights_indiv, weight_total, 
                      bounds_indiv) {
  
  if(any(weights_indiv <= 0)) stop("weights_indiv must all be positive")
  if(any(bounds_indiv  <= 0)) stop("bounds_indiv must all be positive")
  
  K <- length(preds_indiv)
  z <- pred_total - sum(preds_indiv)
  
  A <- rbind(diag(sqrt(weights_indiv)), rep(sqrt(weight_total), K))
  
  U <- matrix(rep(bounds_indiv * weight_total, K), nrow = K, byrow=FALSE) + 
    diag(bounds_indiv * weights_indiv)
  
  v <- weight_total * z * bounds_indiv
  
  Uinv <- diag(1 / (weights_indiv * bounds_indiv)) - 
    1 / (1 / weight_total + sum(1 / weights_indiv)) *
    (1 / weights_indiv) %*% t(1 / (weights_indiv * bounds_indiv))
  
  D    <-   A %*% Uinv
  b    <- - D %*% v
  
  # Lasso shooting divides the least square part by two which
  #  is equivalent to multipliying lambda by two so we had to
  #  divide lambda by two in order to compensate
  t <- lassoshooting(x = D, y = b, lambda = 1, thr = 1e-8)$coefficients
  s <- Uinv %*% (t + v)
  
  W <- sum( (D %*% t - b)^2 ) + 2 * sum(abs(t))
  V <- W - weight_total * z^2
  
  return(list(preds_indiv = preds_indiv + s, solution = V))
}
