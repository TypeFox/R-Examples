#'Simultaneous generation of multivariate data with arbitrary marginals
#'
#'This package implements a specific method for generating n-dimensional random
#'vectors with given marginal distributions and correlation matrix. The method
#'uses the NORTA(NORmal To Anything) approach which generates a standard normal
#'random vector and then transforms it into a random vector with specified
#'marginal distributions and the RA(Retrospective Approximation) algorithm which
#'is a generic stochastic root-finding algorithm. Data generation is
#'accomplished by first using \code{\link{BoundingRA}} to calculate an
#'intermediate multivariate normal correlation matrix, then the matrix is used
#'to generate samples from multivariate normal distribution. The engine function
#'\code{\link{genNORTARA}} will transforms the normal samples to the wanted data
#'set with specified inputmarginals from users. The function
#'\code{\link{valid_input_cormat}} returns the lower and upper bounds of the
#'mixture pre-specified marginal distributions. The function
#'\code{\link{check_input_cormat}} checks the input target correlation matrix
#'whether it is in the lower and upper bounds, if in the bounds, then the
#'function will return \code{TRUE} it means the input target correaltion matrix
#'is feasible, otherwise, it will print the  elements' positions which are out
#'of bounds and give an error message.
#'
#'@details \tabular{ll}{ Package: \tab nortaRA\cr Type: \tab Package\cr Version:
#'  \tab 1.0.0\cr Date: \tab 2014-12-06\cr License: \tab MIT + file LICENSE }
#'@author Po Su \cr Maintainer: Po Su \email{desolator@@sjtu.edu.cn}
#'@name nortaRA-package
#'@aliases nortaRA
NULL

# Described in Appendeix:NORTA RVG Methods step 5.
v_i_j_estimator <- function(i, j, m_size_history, k,
                             r_solution_mat_history,
                             r_estimator_mat_history){
  if (class(r_solution_mat_history) != "list" || class(r_estimator_mat_history) != "list") {
    stop("r_solution_mat_history and r_estimator_mat_history must be list !")
  }
  if (length(r_solution_mat_history) < k || length(r_estimator_mat_history) < k
       || length(m_size_history) < k) {
    stop( paste0("r_solution_mat_history and r_estimator_mat_history
                 and m_size_history must be all no less than ",k,"!"))
  }
  if (k == 1) return(0)
  r_solu_i_j_vec <- NULL
  for (t in 1:k) {
   r_solu_i_j_vec <- c(r_solu_i_j_vec, r_solution_mat_history[[t]][i,j])
  }
  a <- sum(m_size_history[1:k] *  r_solu_i_j_vec^2)
  b <- sum(m_size_history[1:k]) *  r_estimator_mat_history[[k]][i,j]^2
  res <- ((a - b) / (k - 1))^(1 / 2)
  res
}

# Described in Appendeix:NORTA RVG Methods step 6.
delta_k_i_j <- function(i, j, delta1, k, mk, c1,
                        c2, m_size_history,
                        r_solution_mat_history,
                        r_estimator_mat_history) {

  if (k == 1 || k == 2) return(delta1)
  a <- v_i_j_estimator(i, j, m_size_history,
                      k-1, r_solution_mat_history, r_estimator_mat_history)
  b <- 1 / sum(m_size_history[1:(k - 1)])
  d <- 1 / (c1 * mk)
  res <- c2 * a * ( b + d )^(1/2)
  res
}

# Described in Appendeix:NORTA RVG Methods step 4.
r_k_bar_matrix_i_j <- function(i, j, k, r_k_solution_i_j, mk, m_size_history, r_solution_mat_history){
  if (class(r_solution_mat_history) != "list") {
    stop("r_solution_mat_history  must be list ! ")
  }
  if (length(r_solution_mat_history) < (k - 1) || length(m_size_history) < (k - 1)) {    stop( paste0("r_solution_mat_history and  m_size_history must be all no less than ",k - 1,"!"))
  }
  if (k == 1) return(r_k_solution_i_j)
  r_solu_i_j_vec <- NULL
  for( t in 1:(k - 1)) {
    r_solu_i_j_vec <- c(r_solu_i_j_vec, r_solution_mat_history[[t]][i,j])
  }
  a <-sum(  m_size_history[1:(k - 1)] *  r_solu_i_j_vec ) + mk * r_k_solution_i_j
  b <-sum(  m_size_history[1:(k - 1)] ) + mk
  res <- a / b
  res
}
