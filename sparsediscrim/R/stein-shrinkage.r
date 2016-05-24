#' Bias correction function from Pang et al. (2009).
#'
#' This function computes the function \eqn{h_{nu, p}(t)} on page 1023 of Pang
#' et al. (2009).
#'
#' @param nu a specified constant (nu = N - K)
#' @param p the feature space dimension.
#' @param t a constant specified by the user that indicates the exponent to use
#' with the variance estimator. By default, t = -1 as in Pang et al. See the
#' paper for more details.
#'
#' @references Pang, H., Tong, T., & Zhao, H. (2009). "Shrinkage-based Diagonal
#' Discriminant Analysis and Its Applications in High-Dimensional Data,"
#' Biometrics, 65, 4, 1021-1029.
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2009.01200.x/abstract}
#' @return the bias correction value
h <- function(nu, p, t = -1) {
  if (nu <= 0) {
    stop("The value for 'nu' must be positive.")
  }
	
	if(nu / 2 + t / p > 0) {
		(nu / 2)^t * (gamma(nu / 2) / gamma(nu / 2 + t / p))^p
	} else {
		warning("The bias-correction has resulted in a NaN. Incrementing nu by 1...")
		(nu / 2)^t * (gamma(nu / 2) / gamma((nu + 1) / 2 + t / p))^p
	}
}

#' Stein Risk function from Pang et al. (2009).
#'
#' This function finds the value for \eqn{alpha \in [0,1]} that empirically
#' minimizes the average risk under a Stein loss function, which is given on
#' page 1023 of Pang et al. (2009).
#'
#' @param N the sample size.
#' @param K the number of classes.
#' @param var_feature a vector of the sample variances for each dimension.
#' @param num_alphas The number of values used to find the optimal amount of
#' shrinkage.
#' @param t a constant specified by the user that indicates the exponent to use
#' with the variance estimator. By default, t = -1 as in Pang et al. See the
#' paper for more details.
#' @references Pang, H., Tong, T., & Zhao, H. (2009). "Shrinkage-based Diagonal
#' Discriminant Analysis and Its Applications in High-Dimensional Data,"
#' Biometrics, 65, 4, 1021-1029.
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2009.01200.x/abstract}
#' @return list with
#' \itemize{
#'   \item \code{alpha}: the alpha that minimizes the average risk under a Stein
#'    loss function. If the minimum is not unique, we randomly select an
#'    \code{alpha} from the minimizers.
#'   \item \code{risk}: the minimum average risk attained.
#' }
risk_stein <- function(N, K, var_feature, num_alphas = 101, t = -1) {
	nu <- N - K
	p <- length(var_feature)
	alphas <- seq(0, 1, length = num_alphas)
	
	# The pooled variance is defined in Pang et al. (2009) as the geometric mean
	# of the sample variances of each feature.
	var_pool <- prod(var_feature)^(t / p)
	
	# Here we compute the average risk for the Stein loss function on page 1023
	# for all values of alpha.
	risk_alphas <- sapply(alphas, function(alpha) {
		risk <- h(nu = nu, p = p)^alpha * h(nu = nu, p = 1)^(1 - alpha)
		risk <- risk / (h(nu = nu, p = 1, t = alpha * t / p))^(p - 1)
		risk <- risk / h(nu = nu, p = 1, t = (1 - alpha + alpha / p) * t)
		risk <- risk * (var_pool)^(alpha * t)
		risk <- risk * mean(var_feature^(-alpha * t))
		risk <- risk - log(h(nu = nu, p = p)^alpha * h(nu = nu, p = 1)^(1 - alpha))
		risk <- risk - t * digamma(nu / 2)
		risk <- risk + t * log(nu / 2) - 1
		risk
	})

	# Which of the alphas empirically minimize this risk?
	# If there are ties in the minimum risk, we randomly select
	# the value of alpha from the minimizers.
	alpha_min_risk <- alphas[which(min(risk_alphas) == risk_alphas)]
	alpha_star <- sample(alpha_min_risk, 1)
	
	list(alpha = alpha_star, var_pool = var_pool)
}

#' Shrinkage-based estimator of variances for each feature from Pang et al.
#' (2009).
#'
#' This function computes the shrinkage-based estimator of variance of each
#' feature (variable) from Pang et al. (2009) for the SDLDA classifier.
#'
#' @param N the sample size.
#' @param K the number of classes.
#' @param var_feature a vector of the sample variances for each feature.
#' @param num_alphas The number of values used to find the optimal amount of
#' shrinkage.
#' @param t a constant specified by the user that indicates the exponent to use
#' with the variance estimator. By default, t = -1 as in Pang et al. See the
#' paper for more details.
#'
#' @references Pang, H., Tong, T., & Zhao, H. (2009). "Shrinkage-based Diagonal
#' Discriminant Analysis and Its Applications in High-Dimensional Data,"
#' Biometrics, 65, 4, 1021-1029.
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2009.01200.x/abstract}
#' @return a vector of the shrunken variances for each feature.
var_shrinkage <- function(N, K, var_feature, num_alphas = 101, t = -1) {
	nu <- N - K
	p <- length(var_feature)
	
	risk_stein_out <- risk_stein(N = N, K = K, var_feature = var_feature,
                               num_alphas = num_alphas, t = t)

	var_pool <- risk_stein_out$var_pool
	alpha <- risk_stein_out$alpha
	
	var_feature_shrink <- (h(nu = nu, p = p, t = t) * var_pool)^alpha *
    (h(nu = nu, p = 1, t = t) * var_feature)^(1 - alpha)
	var_feature_shrink
}
