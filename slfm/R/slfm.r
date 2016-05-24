#' Sparse Latent Factor Model
#'
#' This function is used to fit a Bayesian sparse
#' latent factor model.
#' 
#' @references
#' 1. Duarte, J. D. N. and Mayrink, V. D. (2015). Factor analysis with mixture modeling to evaluate coherent patterns in microarray data. In Interdisciplinary Bayesian Statistics, volume 118 of Springer Proceedings in Mathematics & Statistics, pages 185-195. Springer International Publishing.
#'
#' @param x matrix with the pre-processed data
#' @param a prior shape parameter for Gamma distribution 
#' @param b prior scale parameter for Gamma distribution
#' @param gamma_a prior parameter for Beta distribution
#' @param gamma_b prior parameter for Beta distribution
#' @param omega_0 prior variance of the spike component
#' @param omega_1 prior variance of the slab component
#' @param sample sample size after burn-in
#' @param burnin burn-in size
#' @param lag lag for MCMC
#' @param degenerate use the degenerate version of mixture
#' @return x: data matrix
#' @return q_star: matrix of MCMC chains for q_star parameter
#' @return alpha: summary table of MCMC chains for alpha parameter
#' @return lambda: summary table of MCMC chains for lambda parameter
#' @return sigma: summary table of MCMC chains for sigma parameter
#' @return classification: classification of each alpha (`present`, `marginal`, `absent`)
#' @export
#' @importFrom coda as.mcmc
#' @importFrom stats window
#' @importFrom Rcpp evalCpp
#' @useDynLib slfm
#' @examples
#' mat <- matrix(rnorm(2000), nrow = 20)
#' slfm(mat, sample = 1000)
slfm <- function(
  x, a = 2.1, b = 1.1, gamma_a = 1, gamma_b = 1,
  omega_0 = 0.01, omega_1 = 10, sample = 1000, burnin = round(0.25*sample), lag = 1, degenerate = FALSE) {
  
  ite <- sample + burnin

  # Convert the x input to numeric matrix
  x <- data.matrix(x)

  if(degenerate) {
    res <- slfm_MDN(x, a, b, gamma_a, gamma_b, omega_1, burnin, lag, sample)
  } else {
    res <- slfm_MNN(x, a, b, gamma_a, gamma_b, omega_0, omega_1, burnin, lag, sample)
  }

  after_burnin <- (burnin + 1):ite

  q_star_matrix <- coda::as.mcmc(res[["qstar"]][after_burnin,])
  q_star_matrix <- window(q_star_matrix, thin = lag)
  hpds_q_star <- coda::HPDinterval(q_star_matrix)
  stats_q_star <- summary(q_star_matrix)$statistics
  alpha_clas <- format_classification(class_interval(hpds_q_star))
  alpha_clas_mean <- class_mean(stats_q_star)

  z_matrix <- coda::as.mcmc(res[["z"]][after_burnin,])
  z_matrix <- window(z_matrix, thin = lag)

  alpha_matrix <- coda::as.mcmc(res[["alpha"]][after_burnin,])
  alpha_matrix <- window(alpha_matrix, thin = lag)

  table_alpha <- alpha_estimation(res[["alpha"]][after_burnin,], alpha_clas_mean, res[["z"]][after_burnin,])

  lambda_matrix <- coda::as.mcmc(res[["lambda"]][after_burnin,])
  lambda_matrix <- window(lambda_matrix, thin = lag)
  stats_lambda <- summary(lambda_matrix)$statistics
  hpds_lambda <- coda::HPDinterval(lambda_matrix)
  table_lambda <- cbind(stats_lambda, hpds_lambda)[,-4]
  colnames(table_lambda)[4:5] = c("Upper HPD", "Lower HPD")

  sigma_matrix <- coda::as.mcmc(res[["sigma2"]][after_burnin,])
  sigma_matrix <- window(sigma_matrix, thin = lag)
  stats_sigma <- summary(sigma_matrix)$statistics
  hpds_sigma <- coda::HPDinterval(sigma_matrix)
  table_sigma <- cbind(stats_sigma, hpds_sigma)[,-4]
  colnames(table_sigma)[4:5] = c("Upper HPD", "Lower HPD")


  obj <- list(
    x = x,
    alpha = table_alpha,
    lambda = table_lambda,
    sigma2 = table_sigma,
    alpha_matrix = alpha_matrix,
    lambda_matrix = lambda_matrix,
    sigma2_matrix = sigma_matrix,
    q_star = q_star_matrix,
    z_matrix = z_matrix,
    classification = alpha_clas)
  class(obj) <- "slfm"
  obj
}

print.slfm <- function(x) {
  cat("SLFM object", "\n")
  cat("\n")
  cat("Dimensions","\n")
  cat("- alpha:", nrow(x$alpha),"\n")
  cat("- lambda:", nrow(x$lambda),"\n")
  cat("\n")
  cat("Classification:","\n")
  print(x$classification)
}

alpha_estimation <- function(x, alpha_clas, z_matrix) {
  table_list <- lapply(1:length(alpha_clas), function(i) {
    chain_indicator <- z_matrix[, i] == 1
    if(alpha_clas[i]) {
      chain <- x[chain_indicator, i]
    } else {
      chain <- x[!chain_indicator, i]
    }
    chain.mcmc <- coda::as.mcmc(chain)
    stats <- summary(chain.mcmc)$statistics
    hpds <- coda::HPDinterval(chain.mcmc)
    table <- c(stats, hpds)[-4]
  })
  table <- do.call(rbind, table_list)
  colnames(table)[4:5] = c("Upper HPD", "Lower HPD")
  table
}