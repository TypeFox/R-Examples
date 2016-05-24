#'Computes an intermediate multivariate normal correlation matrix before using
#'the NORTA approach.
#'
#'The function computes an intermediate multivariate correlation matrix with a
#'specific RA(Retrospective Approximation) algorithm called bounding RA.
#'
#'
#'The function computes an intermediate multivariate correlation matrix with a
#'specific RA(Retrospective Approximation) algorithm called bounding RA.Then the
#'result matrix will be used to applying the NORTA approach to get \code{n}
#'observations which have specified marginals and target correlation matrix
#'\code{cor_matrix}.
#'@param cor_matrix The specified input correlation matrix.
#'@inheritParams valid_input_cormat
#'@param m1 The initial sample size.
#'@param c1 The sample-size multiplier(c1>1).
#'@param c2 The step-size multiplier(c2>0).
#'@param delta1 The initial step size(detla1>0).
#'@param sigma0 The standard error tolerance.
#'@param epsilon The initial error tolerance.
#'@param maxit The maximum number of numerical searches.
#'@return An intermediate normal correalation matrix of the same size as
#'  \code{cor_matrix}.
#'
#'@author Po Su
#'@seealso \code{\link{genNORTARA}}, \code{\link{valid_input_cormat}},
#'  \code{\link{check_input_cormat}}
#'@references  Huifen Chen, (2001) Initialization for NORTA: Generation of
#'  Random Vectors with Specified Marginals and Correlations. INFORMS Journal on
#'  Computing \bold{13(4):312-331.}
#'@note In this implementation of the NORTA approach with bounding RA algorithm
#'  according to the reference paper, the initial sample  size is adjusted from
#'  40 to 60, the other parameters remain the same. And the third choice of
#'  random seeds \code{w_bar} is used for efficiency while the Appendix of the paper
#'  uses the first choice.
#'
#'@examples
#'\dontrun{
#'invcdfnames <- c("qt","qpois","qnorm")
#'paramslists <- list(
#'                m1 = list(df = 3),
#'                m2 = list(lambda = 5),
#'                m3 = NULL  # it means list(mean = 0, sd = 1)
#'                  )
#'cor_matrix <- matrix(c(1,0.5,-0.3,0.5,1,0.4,-0.3,0.4,1), 3)
#'Sigma <- BoundingRA(cor_matrix,invcdfnames,paramslists)
#'Sigma
#'invcdfnames <- c("qunif","qunif","qunif")
#'paramslists <- list(
#'                m1 = NULL, #it means list(min = 0, max = 1)
#'                m2 = NULL,
#'                m3 = NULL
#'               )
#'cor_matrix <- matrix(c(1,0.5,-0.3,0.5,1,0.4,-0.3,0.4,1), 3)
#'Sigma <- BoundingRA(cor_matrix,invcdfnames,paramslists)
#'Sigma
#'#NB:For element 0.5 in cor_matrix, the true root should be around 2*sin(0.5*3.14/6).
#'res <- sapply(c(0.5,-0.3,0.4), function(x){2*sin(x*pi/6)})
#'trueroots <- diag(1/2,3,3)
#'trueroots[upper.tri(trueroots)] <- res
#'trueroots <- trueroots + t(trueroots)
#'trueroots
#'abserrors <- abs(Sigma - trueroots)
#'abserrors
#'}
#'@export
BoundingRA <- function(cor_matrix, invcdfnames,
                       paramslists,  m1 = 60,
                        c1 = 2, c2 = 1, delta1 = 0.0001,
                        sigma0 = 0.01, epsilon = 1e+50, maxit = 1000){
  if( length(invcdfnames) != length(paramslists)) {
    stop("inversecdfs should have the same length as paramslists!")
  }
  ndim <- ncol(cor_matrix)
  m_size_history <- NULL
  r_solution_mat_history <- list()
  r_estimator_mat_history <- list()
  var_r_estimator <- matrix(rep(0, ndim * ndim), ndim)
  # initialize to make the while  statement can continue for the first time
  var_r_estimator[1,ndim] <- sigma0 + 0.1
  mk <- m1
  epsilonk <- epsilon
  k <- 0
  while (k < 4 || any(var_r_estimator > sigma0)) {
    k <- k + 1
    # Here choose w_bar the third choice for better effciency while the reference paper
    # choose the first choice for saving computer's storage capacity
    w_k_bar_initial <- matrix(rnorm(mk*ndim, mean = 0, sd = 1), nrow = mk)
    mk_initial <- mk
    add_sample_size_k <-NULL
    fail_search <- TRUE
    while (fail_search){
    if ( k == 1) r_matrix <- cor_matrix
    else   r_matrix <- r_estimator_mat_history[[k-1]]
    r_k_bar_matrix <- matrix(rep(0, ndim * ndim ), ndim)
    r_solution_matrix <- matrix(rep(0, ndim * ndim ), ndim)
    msave <- mk
    # start with new mk in k iteration  with the same initial random number seeds
    w_k_bar <- rbind( w_k_bar_initial, add_sample_size_k )
    y_estimator <- Perform_chol(r_matrix, w_k_bar, invcdfnames,
                                paramslists)
    check_result <- check_y_estimator(y_estimator, r_matrix, msave, mk,
                                      invcdfnames, paramslists)
    y_estimator <- check_result$y_estimator
    mk <- check_result$mk
    l <- r_matrix
    u <- r_matrix
    y_bar_l <- cor_matrix + 1
    y_bar_u <- cor_matrix - 1
    for (i in 1:(ndim - 1)) {
      Breakouter <- FALSE
      for (j in (i + 1):ndim) {
        p <- 1
        delta_k <- delta_k_i_j(i, j, delta1, k, mk, c1,
                              c2, m_size_history,
                              r_solution_mat_history,
                              r_estimator_mat_history)
        fail_search <- FALSE
        while (p < maxit){
        if (y_estimator[i,j] < cor_matrix[i,j]) {
          if (-1 <= r_matrix[i,j] && r_matrix[i,j] < 1) {
            l[i,j] <- r_matrix[i,j]
            y_bar_l[i,j] <- y_estimator[i,j]
          } else {
            r_solution_matrix[i,j] <- 1
            r_k_bar_matrix[i,j] <- r_k_bar_matrix_i_j(i,j,k,r_solution_matrix[i,j],mk,m_size_history, r_solution_mat_history)
          }
        } else {
          if (-1 < r_matrix[i,j] && r_matrix[i,j] <= 1) {
            u[i,j] <- r_matrix[i,j]
            y_bar_u[i,j] <- y_estimator[i,j]

          } else {
            r_solution_matrix[i,j] <- -1
            r_k_bar_matrix[i,j] <- r_k_bar_matrix_i_j(i,j,k,r_solution_matrix[i,j],mk,m_size_history, r_solution_mat_history)
          }
        }
        if ((y_bar_l[i,j] - cor_matrix[i,j]) * (y_bar_u[i,j] - cor_matrix[i,j]) > 0)
        {
          if (y_estimator[i,j] < cor_matrix[i,j]) {

            r_matrix[i,j] <- min(r_matrix[i,j] + delta_k,  1)

          } else {

            r_matrix[i,j] <- max(r_matrix[i,j] - delta_k, -1)

          }
          delta_k <- 2 * delta_k
          r_mat <- matrix(c(1,r_matrix[i,j], r_matrix[i,j],1),2)
          y_estimator[i,j] <- Perform_chol(r_mat, w_k_bar[ ,c(i,j)],
                                           invcdfnames[c(i,j)],
                                           paramslists[c(i,j)] )[1,2]

        } else
          {
            while (p < maxit) {
            r_regula_falsi <- l[i,j] + (u[i,j] - l[i,j]) * (cor_matrix[i,j] - y_bar_l[i,j]) /
              (y_bar_u[i,j] - y_bar_l[i,j])

            if (u[i,j] - l[i,j] < epsilonk)
            {
              r_solution_matrix[i,j] <- r_regula_falsi
              r_k_bar_matrix[i,j] <- r_k_bar_matrix_i_j(i,j,k,r_solution_matrix[i,j],mk,m_size_history, r_solution_mat_history)
              break
            } else {
              r_mat <- matrix(c(1,r_regula_falsi,r_regula_falsi,1),2)
              y_estimator[i,j] <- Perform_chol(r_mat, w_k_bar[ ,c(i,j)],
                                               invcdfnames[c(i,j)],
                                               paramslists[c(i,j)])[1,2]
              if (y_estimator[i,j] < cor_matrix[i,j]) {
                if (-1 <= r_regula_falsi && r_regula_falsi < 1) {
                  l[i,j] <- r_regula_falsi
                  y_bar_l[i,j] <- y_estimator[i,j]

                } else {
                  r_solution_matrix[i,j] <- 1
                  r_k_bar_matrix[i,j] <- r_k_bar_matrix_i_j(i,j,k,r_solution_matrix[i,j],mk,m_size_history, r_solution_mat_history)
                }
              } else {
                if (-1 < r_regula_falsi && r_regula_falsi <= 1) {
                  u[i,j] <- r_regula_falsi
                  y_bar_u[i,j] <- y_estimator[i,j]

                } else {
                  r_solution_matrix[i,j] <- -1
                  r_k_bar_matrix[i,j] <- r_k_bar_matrix_i_j(i,j,k,r_solution_matrix[i,j],mk,m_size_history, r_solution_mat_history)
                }
              }

            }
            p <- p + 1
          }
          break
          }
        p <- p + 1
        }
        if (p >= maxit) {
          fail_search <- TRUE
          Breakouter <-TRUE
          break
         }
       }
      if(Breakouter) break
      }
     if (fail_search) {
     mk <- c1 * mk
     add_sample_size_k <-  matrix(rnorm((mk - mk_initial) * ndim, mean = 0, sd = 1), nrow = mk - mk_initial)
     }
     }
  m_size_history <-  c(m_size_history, mk)
  r_solution_mat_history[[k]] <- r_solution_matrix
  r_estimator_mat_history[[k]] <- r_k_bar_matrix
  for ( i in 1:(ndim - 1))
    for (j in (i + 1):ndim) {
       tmp <- v_i_j_estimator(i, j, m_size_history,k, r_solution_mat_history,
                                                    r_estimator_mat_history)
       var_r_estimator[i,j] <- tmp / sum(m_size_history[1:k])
    }
  epsilonk <- epsilonk / sqrt(c1)
  mk <- c1 * mk
  }
  rou_estimator <- r_estimator_mat_history[[k]]
  rou_estimator <- rou_estimator+t(rou_estimator)
  diag(rou_estimator) <- 1
  if (!corpcor::is.positive.definite(rou_estimator)) {
    warning( "The estimator of the target Normal correlation matrix is not positive definite. It was replaced by  Nearest positive definite matrix by using function Matrix::nearPD !")
    rou_estimator <- as.matrix(Matrix::nearPD(rou_estimator, corr = TRUE, keepDiag = TRUE)$mat)
  }
  rou_estimator
}
