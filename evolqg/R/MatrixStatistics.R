#' Calculate mean values for various matrix statistics
#'
#' Calculates: Mean Squared Correlation, ICV, Autonomy, ConditionalEvolvability, Constraints, Evolvability, Flexibility, Pc1Percent, Respondability.
#' @aliases Autonomy ConditionalEvolvability Constraints Evolvability Flexibility Pc1Percent Respondability
#' @param cov.matrix A covariance matrix
#' @param iterations Number of random vectors to be used in calculating the stochastic statistics
#' @param full.results If TRUE, full distribution of statistics will be returned.
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @return dist Full distribution of stochastic statistics, only if full.resuts == TRUE
#' @return mean Mean value for all statistics
#' @export
#' @references Hansen, T. F., and Houle, D. (2008). Measuring and comparing evolvability
#' and constraint in multivariate characters. Journal of evolutionary
#' biology, 21(5), 1201-19. doi:10.1111/j.1420-9101.2008.01573.x
#' @author Diogo Melo Guilherme Garcia
#' @importFrom Matrix Matrix solve chol
#' @examples
#' cov.matrix <- cov(iris[,1:4])
#' MeanMatrixStatistics(cov.matrix)
#' 
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' #library(doParallel)
#' ##Windows:
#' #cl <- makeCluster(2)
#' #registerDoParallel(cl)
#' ##Mac and Linux:
#' #registerDoParallel(cores = 2)
#' #MeanMatrixStatistics(cov.matrix, parallel = TRUE)
#' @keywords Autonomy
#' @keywords ConditionalEvolvability
#' @keywords Constraints
#' @keywords Evolvability
#' @keywords Flexibility
#' @keywords Pc1Percent
#' @keywords Respondability
MeanMatrixStatistics <- function (cov.matrix, iterations = 1000, full.results = F, parallel = FALSE) {
  matrix.stat.functions = list ('respondability' = Respondability,
                                'evolvability' = Evolvability,
                                'conditional.evolvability' = ConditionalEvolvability,
                                'autonomy' = Autonomy,
                                'flexibility' = Flexibility,
                                'constraints' = Constraints)
  num.traits <- dim (cov.matrix) [1]
  beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
  beta.mat <- apply (beta.mat, 2, Normalize)
  iso.vec <- Normalize (rep(1, num.traits))
  null.dist <- abs (t (iso.vec) %*% beta.mat)
  null.dist <- sort (null.dist)
  crit.value <- null.dist [round (0.95 * iterations)]
  if(full.results) cat ('critical value: ', crit.value, '\n')
  MatrixStatisticsMap <- function (CurrentFunc) CurrentFunc(cov.matrix = cov.matrix, beta.mat = beta.mat)
  stat.dist <- t(laply (matrix.stat.functions, MatrixStatisticsMap, .parallel = parallel))
  stat.dist <- cbind (stat.dist, null.dist)
  colnames (stat.dist) <- c('respondability',
                            'evolvability',
                            'conditional.evolvability',
                            'autonomy',
                            'flexibility',
                            'constraints',
                            'null.dist')
  stat.mean <- colMeans (stat.dist[,-7])
  integration <- c (CalcR2 (cov2cor(cov.matrix)), Pc1Percent (cov.matrix), CalcICV(cov.matrix))
  names (integration) <- c ('MeanSquaredCorrelation', 'pc1%', 'ICV')
  stat.mean <- c (integration, stat.mean)
  if(full.results)
    return (list ('dist' = stat.dist, 'mean' = stat.mean))
  else
    return (stat.mean)
}

#' @export
#' @importFrom Matrix Matrix solve chol nearPD diag
Autonomy <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  cov.matrix <- Matrix(cov.matrix)
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  tryCatch({cv = chol(cov.matrix)}, error = function(cond){
    warning("matrix is singular, can't compute autonomy directly. Using nearPD, results could be wrong")
  })
   (1/diag(t (beta.mat) %*% solve (cov.matrix, beta.mat))) / diag(t (beta.mat) %*% cov.matrix %*% beta.mat)
} 
#' @export
#' @importFrom Matrix Matrix solve chol nearPD diag
ConditionalEvolvability <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  cov.matrix <- Matrix(cov.matrix)
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  tryCatch({chol(cov.matrix)}, error = function(cond){
    warning("matrix is singular, can't compute conditional evolvability directly. Using nearPD, results could be wrong")
    cov.matrix <- nearPD(cov.matrix)[[1]]
    chol(cov.matrix)
  })
  return (1/diag(t (beta.mat) %*% solve (cov.matrix, beta.mat)))
}
#' @export
Constraints <- function  (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  PC1 <- eigen(cov.matrix)$vectors[,1]
  Cb <- apply(cov.matrix %*% beta.mat, 2, Normalize)
  abs(apply(Cb, 2, `%*%`, PC1))
}
#' @export
Evolvability <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  } 
  diag(t(beta.mat) %*% cov.matrix %*% beta.mat)
}
#' @export
Flexibility <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  Cb <- apply(cov.matrix %*% beta.mat, 2, Normalize)
  diag(t (beta.mat) %*% Cb)
  }
#' @export
Pc1Percent <- function (cov.matrix) return (eigen (cov.matrix)$values [1] / sum (eigen (cov.matrix)$values))
#' @export
Respondability <- function (cov.matrix, beta.mat = NULL, iterations = 1000) {
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  apply(cov.matrix %*% beta.mat, 2, Norm)
}
