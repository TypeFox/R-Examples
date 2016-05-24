#' @name FB
#' @aliases FB
#' @title Fuller-Burmeister estimator
#' 
#' @description Produces estimates for population totals and means using the Fuller - Burmeister estimator from survey data obtained
#'  from a dual frame sampling desing. Confidence intervals are also computed, if required.
#' 
#' @usage FB(ysA, ysB, pi_A, pi_B, domains_A, domains_B, conf_level = NULL)
#' @param ysA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{c} containing information about variable of interest from \eqn{s_A}.
#' @param ysB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{c} containing information about variable of interest from \eqn{s_B}.
#' @param pi_A A numeric vector of length \eqn{n_A} or a square numeric matrix of dimension \eqn{n_A} containing first order or first and second order inclusion probabilities for units included in \eqn{s_A}.
#' @param pi_B A numeric vector of length \eqn{n_B} or a square numeric matrix of dimension \eqn{n_B} containing first order or first and second order inclusion probabilities for units included in \eqn{s_B}.
#' @param domains_A A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domains_B A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".
#' @param conf_level (Optional) A numeric value indicating the confidence level for the confidence intervals.
#' @details Fuller-Burmeister estimator of population total is given by
#'  \deqn{\hat{Y}_{FB} = \hat{Y}_a^A + \hat{\beta_1}\hat{Y}_{ab}^A + (1 - \hat{\beta_1})\hat{Y}_{ab}^B + \hat{Y}_b^B + \hat{\beta_2}(\hat{N}_{ab}^A - \hat{N}_{ab}^B)}
#'  where optimal values for \eqn{\hat{\beta}} to minimize variance of the estimator are:
#'  \deqn{
#'    \left( \begin{array}{c}
#'        \hat{\beta}_1\\
#'        \hat{\beta}_2
#'    \end{array} \right)
#'        = -
#'    \left( \begin{array}{cc}
#'        \hat{V}(\hat{Y}_{ab}^A - \hat{Y}_{ab}^B) & \widehat{Cov}(\hat{Y}_{ab}^A - \hat{Y}_{ab}^B, \hat{N}_{ab}^A - \hat{N}_{ab}^B)\\
#'        \widehat{Cov}(\hat{Y}_{ab}^A - \hat{Y}_{ab}^B, \hat{N}_{ab}^A - \hat{N}_{ab}^B) & \hat{V}(\hat{N}_{ab}^A - \hat{N}_{ab}^B)
#'    \end{array} \right)^{-1}
#'        \times}
#'  \deqn{
#'    \left( \begin{array}{c}
#'        \widehat{Cov}(\hat{Y}_a^A + \hat{Y}_b^B + \hat{Y}_{ab}^B, \hat{Y}_{ab}^A - \hat{Y}_{ab}^B)\\
#'        \widehat{Cov}(\hat{Y}_a^A + \hat{Y}_b^B + \hat{Y}_{ab}^B, \hat{N}_{ab}^A - \hat{N}_{ab}^B)
#'    \end{array} \right)
#'    }
#'  Due to Fuller-Burmeister estimator is not defined for estimating population sizes, estimation of the mean is computed as \eqn{\hat{Y}_{FB} / \hat{N}_H}, where \eqn{\hat{N}_H}
#'  is the estimation of the population size using Hartley estimator.
#'  Estimated variance for the Fuller-Burmeister estimator can be obtained through expression
#'  \deqn{\hat{V}(\hat{Y}_{FB}) = \hat{V}(\hat{Y}_a^A) + \hat{V}(\hat{Y}^B) + 
#'     \hat{\beta}_1[\widehat{Cov}(\hat{Y}_a^A, \hat{Y}_{ab}^A) - \widehat{Cov}(\hat{Y}^B, \hat{Y}_{ab}^B)]} 
#'  \deqn{ + \hat{\beta}_2[\widehat{Cov}(\hat{Y}_a^A, \hat{N}_{ab}^A) - \widehat{Cov}(\hat{Y}^B, \hat{N}_{ab}^B)]
#'  }
#'  If both first and second order probabilities are known, variances and covariances involved in calculation of \eqn{\hat{\beta}} and \eqn{\hat{V}(\hat{Y}_{FB})} are estimated using functions \code{VarHT} and \code{CovHT}, respectively. If
#'  only first order probabilities are known, variances are estimated using Deville's method and covariances are estimated using following expression
#'  \deqn{\widehat{Cov}(\hat{X}, \hat{Y}) = \frac{\hat{V}(X + Y) - \hat{V}(X) - \hat{V}(Y)}{2}}
#' @return \code{FB} returns an object of class "EstimatorDF" which is a list with, at least, the following components:
#'  \item{Call}{the matched call.}
#'  \item{Est}{total and mean estimation for main variable(s).}
#'  \item{VarEst}{variance estimation for main variable(s).}
#'  If parameter \code{conf_level} is different from \code{NULL}, object includes component
#'  \item{ConfInt}{total and mean estimation and confidence intervals for main variables(s).}
#'  In addition, components \code{TotDomEst} and \code{MeanDomEst} are available when estimator is based on estimators of the domains. Component \code{Param} shows value of parameters involded in calculation of the estimator (if any).
#'  By default, only \code{Est} component (or \code{ConfInt} component, if parameter \code{conf_level} is different from \code{NULL}) is shown. It is possible to access to all the components of the objects by using function \code{summary}.
#' @references Fuller, W.A. and Burmeister, L.F. (1972).
#'  \emph{Estimation for Samples Selected From Two Overlapping Frames} ASA Proceedings of the Social Statistics Sections, 245 - 249.
#' @seealso \code{\link{Hartley}} \code{\link{JackFB}}
#' @examples
#' data(DatA)
#' data(DatB)
#' data(PiklA)
#' data(PiklB)
#'
#' #Let calculate Fuller-Burmeister estimator for variable Clothing
#' FB(DatA$Clo, DatB$Clo, PiklA, PiklB, DatA$Domain, DatB$Domain)
#' 
#' #Now, let calculate Fuller-Burmeister estimator and a 90% confidence interval
#' #for variable Leisure, considering only first order inclusion probabilities
#' FB(DatA$Lei, DatB$Lei, DatA$ProbA, DatB$ProbB, DatA$Domain, 
#' DatB$Domain, 0.90)
#' @export
FB = function (ysA, ysB, pi_A, pi_B, domains_A, domains_B, conf_level = NULL)
{
	cnames <- names(ysA)
	ysA <- as.matrix(ysA)
	ysB <- as.matrix(ysB)
	pi_A <- as.matrix(pi_A)
	pi_B <- as.matrix(pi_B)

	if (any(is.na(ysA)))
		stop("There are missing values in sample from frame A.")
	if (any(is.na(ysB)))
		stop("There are missing values in sample from frame B.")
	if (any(is.na(pi_A)))
		stop("There are missing values in pikl from frame A.")
	if (any(is.na(pi_B)))
		stop("There are missing values in pikl from frame B.")
	if (any(is.na(domains_A)))
		stop("There are missing values in domains from frame A.")
	if (any(is.na(domains_B)))
		stop("There are missing values in domains from frame B.")
	if (ncol(ysA) != ncol(ysB))
		stop("Number of variables does not match.")
	if (nrow(ysA) != nrow(pi_A) | nrow(ysA) != length(domains_A) | length(domains_A) != nrow(pi_A))
		stop("Arguments from frame A have different sizes.")
	if (nrow(ysB) != nrow(pi_B) | nrow(ysB) != length(domains_B) | length(domains_B) != nrow(pi_B))
		stop("Arguments from frame B have different sizes.")
	if (length(which(domains_A == "a")) + length(which(domains_A == "ab")) != length(domains_A))
		stop("Domains from frame A are not correct.")
	if (length(which(domains_B == "b")) + length(which(domains_B == "ba")) != length(domains_B))
		stop("Domains from frame B are not correct.")

	cl <- match.call()

	n_A <- nrow(ysA)
  	n_B <- nrow(ysB)
	c <- ncol(ysA)

	ones_ab_A <- Domains (rep (1, n_A), domains_A, "ab")
	ones_ab_B <- Domains (rep (1, n_B), domains_B, "ba")

	est <- matrix(, 2, c, dimnames = list(c("Total", "Mean"), cnames))
	varest <- matrix(, 2, c, dimnames = list(c("Var. Total", "Var. Mean"), cnames))
	totdom <- matrix(, 4, c, dimnames = list(c("Total dom. a", "Total dom. ab", "Total dom. b", "Total dom. ba"), cnames))
	meandom <- matrix(, 4, c, dimnames = list(c("Mean dom. a", "Mean dom. ab", "Mean dom. b", "Mean dom. ba"), cnames))
	par <- 	matrix(, 2, c, dimnames = list(c("beta1", "beta2"), cnames))
	if (is.null(conf_level))
		interv <- NULL
	else
		interv <- matrix(, 6, c, dimnames = list(c("Total", "Lower Bound", "Upper Bound", "Mean", "Lower Bound", "Upper Bound"), cnames))

	H <- Hartley (rep(1, nrow(ysA)), rep(1, nrow(ysB)), pi_A, pi_B, domains_A, domains_B)
	size_estimation <- H$Est[1,1]
	domain_size_estimation <- H$TotDomEst[,1]

	if (!is.null(dim(drop(pi_A))) & !is.null(dim(drop(pi_B)))) {

		Nhat_ab_A <- HT (ones_ab_A, diag(pi_A))
		Nhat_ab_B <- HT (ones_ab_B, diag(pi_B))
		Vhat_Nhat_ab_A <- VarHT (ones_ab_A, pi_A)
		Vhat_Nhat_ab_B <- VarHT (ones_ab_B, pi_B)

		for (k in 1:c) {

			if (nrow(pi_A) != ncol(pi_A))
				stop("Pikl from frame A is not a square matrix.")
			if (nrow(pi_B) != ncol(pi_B))
				stop("Pikl from frame B is not a square matrix.")

			data_a_A <- Domains (ysA[,k], domains_A, "a")
			data_ab_A <- Domains (ysA[,k], domains_A, "ab")
			data_b_B <- Domains (ysB[,k], domains_B, "b")
			data_ab_B <- Domains (ysB[,k], domains_B, "ba")	

			Yhat_a_A <- HT (data_a_A, diag(pi_A))
			Yhat_ab_A <- HT (data_ab_A, diag(pi_A))
			Yhat_b_B <- HT (data_b_B, diag(pi_B))
			Yhat_ab_B <- HT (data_ab_B, diag(pi_B))
			Vhat_Yhat_a_A <- VarHT(data_a_A, pi_A)
			Vhat_Yhat_ab_A <- VarHT (data_ab_A, pi_A)
			Vhat_Yhat_b_B <- VarHT(data_b_B, pi_B)
			Vhat_Yhat_ab_B <- VarHT (data_ab_B, pi_B)

			Covhat_Yhat_a_A_Yhat_ab_A <- CovHT(data_a_A, data_ab_A, pi_A)
			Covhat_Yhat_b_B_Yhat_ab_B <- CovHT(data_b_B, data_ab_B, pi_B)
			Covhat_Yhat_a_A_Nhat_ab_A <- CovHT(data_a_A, ones_ab_A, pi_A)
			Covhat_Yhat_b_B_Nhat_ab_B <- CovHT(data_b_B, ones_ab_B, pi_B)
			Covhat_Yhat_ab_A_Nhat_ab_A <- CovHT(data_ab_A, ones_ab_A, pi_A)
			Covhat_Yhat_ab_B_Nhat_ab_B <- CovHT(data_ab_B, ones_ab_B, pi_B)
		
			mat <- matrix (0, nrow = 2, ncol = 2)
			mat[1,1] <- Vhat_Yhat_ab_A + Vhat_Yhat_ab_B
			mat[1,2] <- Covhat_Yhat_ab_A_Nhat_ab_A + Covhat_Yhat_ab_B_Nhat_ab_B
			mat[2,1] <- Covhat_Yhat_ab_A_Nhat_ab_A + Covhat_Yhat_ab_B_Nhat_ab_B
			mat[2,2] <- Vhat_Nhat_ab_A + Vhat_Nhat_ab_B
		
			vec <- matrix (0, nrow = 2, ncol = 1)
			vec[1,1] <- Covhat_Yhat_a_A_Yhat_ab_A - Covhat_Yhat_b_B_Yhat_ab_B - Vhat_Yhat_ab_B
			vec[2,1] <- Covhat_Yhat_a_A_Nhat_ab_A - Covhat_Yhat_b_B_Nhat_ab_B - Covhat_Yhat_ab_B_Nhat_ab_B

			beta <- -solve(mat,vec)

			totdom[,k] <- c(Yhat_a_A, Yhat_ab_A, Yhat_b_B, Yhat_ab_B)
			meandom[,k] <- totdom[,k]/domain_size_estimation
			par[,k] <- beta
	
			total_estimation <- Yhat_a_A + Yhat_b_B + beta[1] * Yhat_ab_A + (1 - beta[1]) * Yhat_ab_B + beta[2] * (Nhat_ab_A - Nhat_ab_B)
			mean_estimation <- total_estimation / size_estimation
			est[,k] <- c(total_estimation, mean_estimation)

			Vhat_Yhat_a_A <- VarHT (data_a_A, pi_A)
			Vhat_Yhat_b_B <- VarHT (data_b_B, pi_B)
			Vhat_Yhat_B <- Vhat_Yhat_b_B + Vhat_Yhat_ab_B + 2 * Covhat_Yhat_b_B_Yhat_ab_B
			Vhat_Yhat_FB <- Vhat_Yhat_a_A + Vhat_Yhat_B + beta[1] * (Covhat_Yhat_a_A_Yhat_ab_A - Covhat_Yhat_b_B_Yhat_ab_B - Vhat_Yhat_ab_B) + beta[2] * (Covhat_Yhat_a_A_Nhat_ab_A - Covhat_Yhat_b_B_Nhat_ab_B - Covhat_Yhat_ab_B_Nhat_ab_B)
			Vhat_Ymeanhat_FB <- 1/size_estimation^2 * Vhat_Yhat_FB
			varest[,k] <- c(Vhat_Yhat_FB, Vhat_Ymeanhat_FB)

			if (!is.null(conf_level)) {

				total_upper <- total_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_FB)
				total_lower <- total_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_FB)
				mean_upper <- mean_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_FB)
				mean_lower <- mean_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_FB)
				interv[,k] <- c(total_estimation, total_lower, total_upper, mean_estimation, mean_lower, mean_upper)	
			}
		}
	}
	else {

		if (is.null(dim(drop(pi_A))) & is.null(dim(drop(pi_B)))){
			
			Nhat_ab_A <- HT (ones_ab_A, pi_A)
			Nhat_ab_B <- HT (ones_ab_B, pi_B)
			Vhat_Nhat_ab_A <- varest (ones_ab_A, pik = pi_A)
			Vhat_Nhat_ab_B <- varest (ones_ab_B, pik = pi_B)

			for (k in 1:c) {

				data_a_A <- Domains (ysA[,k], domains_A, "a")
				data_ab_A <- Domains (ysA[,k], domains_A, "ab")
				data_b_B <- Domains (ysB[,k], domains_B, "b")
				data_ab_B <- Domains (ysB[,k], domains_B, "ba")	

				Yhat_a_A <- HT (data_a_A, pi_A)
				Yhat_ab_A <- HT (data_ab_A, pi_A)
				Yhat_b_B <- HT (data_b_B, pi_B)
				Yhat_ab_B <- HT (data_ab_B, pi_B)
				Vhat_Yhat_a_A <- varest(data_a_A, pik = pi_A)
				Vhat_Yhat_ab_A <- varest (data_ab_A, pik = pi_A)
				Vhat_Yhat_b_B <- varest(data_b_B, pik = pi_B)
				Vhat_Yhat_ab_B <- varest (data_ab_B, pik = pi_B)

				Covhat_Yhat_ab_A_Nhat_ab_A <- (varest (Ys = data_ab_A + ones_ab_A, pik = pi_A) - Vhat_Yhat_ab_A - Vhat_Nhat_ab_A) / 2
				Covhat_Yhat_ab_B_Nhat_ab_B <- (varest (Ys = data_ab_B + ones_ab_B, pik = pi_B) - Vhat_Yhat_ab_B - Vhat_Nhat_ab_B) / 2
				Covhat_Yhat_a_A_Yhat_ab_A <- (varest (Ys = ysA[,k], pik = pi_A) - Vhat_Yhat_a_A - Vhat_Yhat_ab_A) / 2
				Covhat_Yhat_b_B_Yhat_ab_B <- (varest (Ys = ysB[,k], pik = pi_B) - Vhat_Yhat_b_B - Vhat_Yhat_ab_B) / 2
				dat <- ysA[,k]; dat[domains_A == "ab"] <- 1
				Covhat_Yhat_a_A_Nhat_ab_A <- (varest (Ys = dat, pik = pi_A) - Vhat_Yhat_a_A - Vhat_Nhat_ab_A) / 2
				dat <- ysB[,k]; dat[domains_B == "ba"] <- 1
				Covhat_Yhat_b_B_Nhat_ab_B <- (varest (Ys = dat, pik = pi_B) - Vhat_Yhat_b_B - Vhat_Nhat_ab_B) / 2
		
				mat <- matrix (0, nrow = 2, ncol = 2)
				mat[1,1] <- Vhat_Yhat_ab_A + Vhat_Yhat_ab_B
				mat[1,2] <- Covhat_Yhat_ab_A_Nhat_ab_A + Covhat_Yhat_ab_B_Nhat_ab_B
				mat[2,1] <- mat [1,2]
				mat[2,2] <- Vhat_Nhat_ab_A + Vhat_Nhat_ab_B
		
				vec <- matrix (0, nrow = 2, ncol = 1)
				vec[1,1] <- Covhat_Yhat_a_A_Yhat_ab_A - Covhat_Yhat_b_B_Yhat_ab_B - Vhat_Yhat_ab_B
				vec[2,1] <- Covhat_Yhat_a_A_Nhat_ab_A - Covhat_Yhat_b_B_Nhat_ab_B - Covhat_Yhat_ab_B_Nhat_ab_B

				beta <- -solve(mat, vec)

				totdom[,k] <- c(Yhat_a_A, Yhat_ab_A, Yhat_b_B, Yhat_ab_B)
				meandom[,k] <- totdom[,k]/domain_size_estimation
				par[,k] <- beta

				total_estimation <- Yhat_a_A + Yhat_b_B + beta[1] * Yhat_ab_A + (1 - beta[1]) * Yhat_ab_B + beta[2] * (Nhat_ab_A - Nhat_ab_B)
				mean_estimation <- total_estimation / size_estimation
				est[,k] <- c(total_estimation, mean_estimation)

      				Vhat_Yhat_a_A <- varest (data_a_A, pik = pi_A)
				Vhat_Yhat_b_B <- varest (data_b_B, pik = pi_B)
				Vhat_Yhat_B <- Vhat_Yhat_b_B + Vhat_Yhat_ab_B + 2 * Covhat_Yhat_b_B_Yhat_ab_B
				Vhat_Yhat_FB <- Vhat_Yhat_a_A + Vhat_Yhat_B + beta[1] * (Covhat_Yhat_a_A_Yhat_ab_A - Covhat_Yhat_b_B_Yhat_ab_B - Vhat_Yhat_ab_B) + beta[2] * (Covhat_Yhat_a_A_Nhat_ab_A - Covhat_Yhat_b_B_Nhat_ab_B - Covhat_Yhat_ab_B_Nhat_ab_B)
				Vhat_Ymeanhat_FB <- 1/size_estimation^2 * Vhat_Yhat_FB
				varest[,k] <- c(Vhat_Yhat_FB, Vhat_Ymeanhat_FB)

				if (!is.null(conf_level)) {

					total_upper <- total_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_FB)
					total_lower <- total_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_FB)
					mean_upper <- mean_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_FB)
					mean_lower <- mean_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_FB)
					interv[,k] <- c(total_estimation, total_lower, total_upper, mean_estimation, mean_lower, mean_upper)			
				}
			}
		}
		else
			stop("Invalid option: Probability vector in one frame and probability matrix in the other frame. Type of probabilities structures must match.")
	}
	results = list(Call = cl, Est = est, VarEst = varest, TotDomEst = totdom, MeanDomEst = meandom, Param = par, ConfInt = interv)
   	class(results) = "EstimatorDF"
   	attr(results, "attributesDF") = conf_level
   	return(results)
}