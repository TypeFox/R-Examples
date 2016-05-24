#' @name BKA
#' @aliases BKA
#' @title Bankier-Kalton-Anderson estimator
#' 
#' @description Produces estimates for population total and mean using the Bankier-Kalton-Anderson estimator from survey data obtained
#'  from a dual frame sampling design. Confidence intervals are also computed, if required.
#' 
#' @usage BKA(ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, 
#' conf_level = NULL)
#' @param ysA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{c} containing information about variable(s) of interest from \eqn{s_A}.
#' @param ysB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{c} containing information about variable(s) of interest from \eqn{s_B}.
#' @param pi_A A numeric vector of length \eqn{n_A} or a square numeric matrix of dimension \eqn{n_A} containing first order or first and second order inclusion probabilities for units included in \eqn{s_A}.
#' @param pi_B A numeric vector of length \eqn{n_B} or a square numeric matrix of dimension \eqn{n_B} containing first order or first and second order inclusion probabilities for units included in \eqn{s_B}.
#' @param pik_ab_B A numeric vector of size \eqn{n_A} containing first order inclusion probabilities according to sampling design in frame B for units belonging 
#'  to overlap domain that have been selected in \eqn{s_A}.
#' @param pik_ba_A A numeric vector of size \eqn{n_B} containing first order inclusion probabilities according to sampling design in frame A for units belonging 
#'  to overlap domain that have been selected in \eqn{s_B}.
#' @param domains_A A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domains_B A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".
#' @param conf_level (Optional) A numeric value indicating the confidence level for the confidence intervals, if desired.
#' @details BKA estimator of population total is given by
#'  \deqn{\hat{Y}_{BKA} = \sum_{i \in s_A}\tilde{d}_i^Ay_i + \sum_{i \in s_B}\tilde{d}_i^By_i}
#'  where
#'  \eqn{\tilde{d}_i^A =\left\{\begin{array}{lcc}
#'  d_i^A & \textrm{if } i \in a\\
#'  (1/d_i^A + 1/d_i^B)^{-1} & \textrm{if } i \in ab
#'  \end{array}
#'  \right.}
#'  and
#'  \eqn{\tilde{d}_i^B =\left\{\begin{array}{lcc}
#'  d_i^B & \textrm{if } i \in b\\
#'  (1/d_i^A + 1/d_i^B)^{-1} & \textrm{if } i \in ba
#'  \end{array}
#'  \right.}
#'  being \eqn{d_i^A} and \eqn{d_i^B} the design weights, obtained as the inverse of the first order inclusion probabilities, that is, \eqn{d_i^A = 1/\pi_i^A} and \eqn{d_i^B = 1/\pi_i^B}.
#'  
#'  To estimate variance of this estimator, one uses following approach proposed by Rao and Skinner (1996)
#'  \deqn{\hat{V}(\hat{Y}_{BKA}) = \hat{V}(\sum_{i \in s_A}\tilde{z}_i^A) + \hat{V}(\sum_{i \in s_B}\tilde{z}_i^B)}
#'  with \eqn{\tilde{z}_i^A = \delta_i(a)y_i + (1 - \delta_i(a))y_i\pi_i^A/(\pi_i^A + \pi_i^B)} and \eqn{\tilde{z}_i^B = \delta_i(b)y_i + (1 - \delta_i(b))y_i\pi_i^B/(\pi_i^A + \pi_i^B)},
#'  being \eqn{\delta_i(a)} and \eqn{\delta_i(b)} the indicator variables for domain \eqn{a} and domain \eqn{b}, respectively.
#'  If both first and second order probabilities are known, variances and covariances involved in calculation of \eqn{\hat{\beta}} and \eqn{\hat{V}(\hat{Y}_{FB})} are estimated using functions \code{VarHT} and \code{CovHT}, respectively. If
#'  only first order probabilities are known, variances are estimated using Deville's method and covariances are estimated using following expression
#'  \deqn{\widehat{Cov}(\hat{X}, \hat{Y}) = \frac{\hat{V}(X + Y) - \hat{V}(X) - \hat{V}(Y)}{2}}
#' @return \code{BKA} returns an object of class "EstimatorDF" which is a list with, at least, the following components:
#'  \item{Call}{the matched call.}
#'  \item{Est}{total and mean estimation for main variable(s).}
#'  \item{VarEst}{variance estimation for main variable(s).}
#'  If parameter \code{conf_level} is different from \code{NULL}, object includes component
#'  \item{ConfInt}{total and mean estimation and confidence intervals for main variables(s).}
#'  In addition, components \code{TotDomEst} and \code{MeanDomEst} are available when estimator is based on estimators of the domains. Component \code{Param} shows value of parameters involded in calculation of the estimator (if any).
#'  By default, only \code{Est} component (or \code{ConfInt} component, if parameter \code{conf_level} is different from \code{NULL}) is shown. It is possible to access to all the components of the objects by using function \code{summary}.
#' @references Bankier, M. D. (1986)
#'  \emph{Estimators Based on Several Stratified Samples With Applications to Multiple Frame Surveys}. 
#'  Journal of the American Statistical Association, Vol. 81, 1074 - 1079.
#' @references Kalton, G. and Anderson, D. W. (1986) 
#'  \emph{Sampling Rare Populations}. 
#'  Journal of the Royal Statistical Society, Ser. A, Vol. 149, 65 - 82.
#' @references Rao, J. N. K. and Skinner, C. J. (1996)
#' \emph{Estimation in Dual Frame Surveys with Complex Designs}.
#' Proceedings of the Survey Method Section, Statistical Society of Canada, 63 - 68.
#' @references Skinner, C. J. and Rao, J. N. K. (1996) 
#'  \emph{Estimation in Dual Frame Surveys with Complex Designs}.
#' Journal of the American Statistical Association, Vol. 91, 433, 349 - 356.
#' @seealso \code{\link{JackBKA}}
#' @examples
#' data(DatA)
#' data(DatB)
#' data(PiklA)
#' data(PiklB)
#' 
#' #Let calculate BKA estimator for population total for variable Leisure
#' BKA(DatA$Lei, DatB$Lei, PiklA, PiklB, DatA$ProbB, DatB$ProbA, 
#' DatA$Domain, DatB$Domain)
#' 
#' #Now, let calculate BKA estimator and a 90% confidence interval for population 
#' #total for variable Feeding considering only first order inclusion probabilities
#' BKA(DatA$Feed, DatB$Feed, DatA$ProbA, DatB$ProbB, DatA$ProbB, 
#' DatB$ProbA, DatA$Domain, DatB$Domain, 0.90)
#' @export
BKA = function (ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, conf_level = NULL)
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
	if (any(is.na(pik_ab_B[domains_A == "ab"])))
		stop("Some values in pik_ab_B are 0 when they should not.")
	if (any(is.na(pik_ba_A[domains_B == "ba"])))
		stop("Some values in pik_ba_A are 0 when they should not.")
	if (nrow(ysA) != nrow(pi_A) | nrow(ysA) != length(domains_A) | length(domains_A) != nrow(pi_A) | nrow(ysA) != length(pik_ab_B))
		stop("Arguments from frame A have different sizes.")
	if (nrow(ysB) != nrow(pi_B) | nrow(ysB) != length(domains_B) | length(domains_B) != nrow(pi_B) | nrow(ysB) != length(pik_ba_A))
		stop("Arguments from frame B have different sizes.")
	if (ncol(ysA) != ncol(ysB))
		stop("Number of variables does not match.")
	if (length(which(domains_A == "a")) + length(which(domains_A == "ab")) != length(domains_A))
		stop("Domains from frame A are not correct.")
	if (length(which(domains_B == "b")) + length(which(domains_B == "ba")) != length(domains_B))
		stop("Domains from frame B are not correct.")

	cl <- match.call()

	n_A <- nrow(ysA)
 	n_B <- nrow(ysB)
	c <- ncol(ysA)

	ysA <- cbind(rep(1, n_A), ysA)
	ysB <- cbind(rep(1, n_B), ysB)

	delta_a_A <- Domains (rep (1, n_A), domains_A, "a")
	delta_ab_A <- Domains (rep (1, n_A), domains_A, "ab")
	delta_b_B <- Domains (rep (1, n_B), domains_B, "b")
	delta_ab_B <- Domains (rep (1, n_B), domains_B, "ba")

	est <- matrix(, 2, c, dimnames = list(c("Total", "Mean"), cnames))
	varest <- matrix(, 2, c, dimnames = list(c("Var. Total", "Var. Mean"), cnames))
	totdom <- NULL
	meandom <- NULL
	par <- NULL
	if (is.null(conf_level))
		interv <- NULL
	else
		interv <- matrix(, 6, c, dimnames = list(c("Total", "Lower Bound", "Upper Bound", "Mean", "Lower Bound", "Upper Bound"), cnames))

	for (k in 1:(c+1)) {

		if (!is.null(dim(drop(pi_A))) & !is.null(dim(drop(pi_B)))) {
		
			if (nrow(pi_A) != ncol(pi_A))
				stop("Pikl from frame A is not a square matrix.")
			if (nrow(pi_B) != ncol(pi_B))
				stop("Pikl from frame B is not a square matrix.")

			w_tilde_iS_A <- (1 / diag(pi_A)) * (domains_A == "a") + (1 / (diag(pi_A) + pik_ab_B)) * (domains_A == "ab")
			w_tilde_iS_B <- (1 / diag(pi_B)) * (domains_B == "b") + (1 / (diag(pi_B) + pik_ba_A)) * (domains_B == "ba")

			if (k == 1){
				Nhat_a_A <- sum(ysA[,k] * w_tilde_iS_A * delta_a_A)
				Nhat_ab_A <- sum(ysA[,k] * w_tilde_iS_A * delta_ab_A)
				Nhat_b_B <- sum(ysB[,k] * w_tilde_iS_B * delta_b_B)
				Nhat_ab_B <- sum(ysB[,k] * w_tilde_iS_B * delta_ab_B)
				domain_size_estimation <- c(Nhat_a_A, Nhat_ab_A, Nhat_b_B, Nhat_ab_B)
				size_estimation <- sum(ysA[,k] * w_tilde_iS_A) + sum(ysB[,k] * w_tilde_iS_B)
			}
			else{
				Yhat_a_A <- sum(ysA[,k] * w_tilde_iS_A * delta_a_A)
				Yhat_ab_A <- sum(ysA[,k] * w_tilde_iS_A * delta_ab_A)
				Yhat_b_B <- sum(ysB[,k] * w_tilde_iS_B * delta_b_B)
				Yhat_ab_B <- sum(ysB[,k] * w_tilde_iS_B * delta_ab_B)
				total_estimation <- sum(ysA[,k] * w_tilde_iS_A) + sum(ysB[,k] * w_tilde_iS_B)
			}

			if (k > 1) {

				mean_estimation <- total_estimation / size_estimation 
				est[,k-1] <- c(total_estimation, mean_estimation)

				p <- diag(pi_A) / (diag(pi_A) + pik_ab_B)
				zA <- delta_a_A * ysA[,k] + (1 - delta_a_A) * ysA[,k] * p			

				q <- diag(pi_B) / (diag(pi_B) + pik_ba_A)
				zB <- delta_b_B * ysB[,k] + (1 - delta_b_B) * ysB[,k] * q

				Vhat_Yhat_SF <- VarHT (zA, pi_A) + VarHT (zB, pi_B)
				Vhat_Ymeanhat_SF <- 1/size_estimation^2 * Vhat_Yhat_SF
				varest[,k-1] <- c(Vhat_Yhat_SF, Vhat_Ymeanhat_SF)

				if (!is.null(conf_level)) {

					total_upper <- total_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_SF)
					total_lower <- total_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_SF)
					mean_upper <- mean_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_SF)
					mean_lower <- mean_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_SF)
					interv[,k-1] <- c(total_estimation, total_lower, total_upper, mean_estimation, mean_lower, mean_upper)			
				}
			}
		}
		else{
		
			if (is.null(dim(drop(pi_A))) & is.null(dim(drop(pi_B)))){
	
				w_tilde_iS_A <- (1 / pi_A) * (domains_A == "a") + (1 / (pi_A + pik_ab_B)) * (domains_A == "ab")
				w_tilde_iS_B <- (1 / pi_B) * (domains_B == "b") + (1 / (pi_B + pik_ba_A)) * (domains_B == "ba")

				if (k == 1){
					Nhat_a_A <- sum(ysA[,k] * w_tilde_iS_A * delta_a_A)
					Nhat_ab_A <- sum(ysA[,k] * w_tilde_iS_A * delta_ab_A)
					Nhat_b_B <- sum(ysB[,k] * w_tilde_iS_B * delta_b_B)
					Nhat_ab_B <- sum(ysB[,k] * w_tilde_iS_B * delta_ab_B)
					domain_size_estimation <- c(Nhat_a_A, Nhat_ab_A, Nhat_b_B, Nhat_ab_B)
					size_estimation <- sum(ysA[,k] * w_tilde_iS_A) + sum(ysB[,k] * w_tilde_iS_B)
				}
				else{
					Yhat_a_A <- sum(ysA[,k] * w_tilde_iS_A * delta_a_A)
					Yhat_ab_A <- sum(ysA[,k] * w_tilde_iS_A * delta_ab_A)
					Yhat_b_B <- sum(ysB[,k] * w_tilde_iS_B * delta_b_B)
					Yhat_ab_B <- sum(ysB[,k] * w_tilde_iS_B * delta_ab_B)
					total_estimation <- sum(ysA[,k] * w_tilde_iS_A) + sum(ysB[,k] * w_tilde_iS_B)
				}

				if (k > 1) {

					mean_estimation <- total_estimation / size_estimation
					est[,k-1] <- c(total_estimation, mean_estimation)
					
					p <- pi_A / (pi_A + pik_ab_B)
					zA <- delta_a_A * ysA[,k] + (1 - delta_a_A) * ysA[,k] * p			

					q <- pi_B / (pi_B + pik_ba_A)
					zB <- delta_b_B * ysB[,k] + (1 - delta_b_B) * ysB[,k] * q

					Vhat_Yhat_SF <- varest (zA, pik = pi_A) + varest (zB, pik = pi_B)
					Vhat_Ymeanhat_SF <- 1/size_estimation^2 * Vhat_Yhat_SF	
					varest[,k-1] <- c(Vhat_Yhat_SF, Vhat_Ymeanhat_SF)

					if (!is.null(conf_level)) {
						
						total_upper <- total_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_SF)
						total_lower <- total_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_SF)
						mean_upper <- mean_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_SF)
						mean_lower <- mean_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_SF)
						interv[,k-1] <- c(total_estimation, total_lower, total_upper, mean_estimation, mean_lower, mean_upper)			
					}
				}
			}
			else

				stop("Invalid option: Probability vector in one frame and probability matrix in the other frame. Type of both structures must match.")
		}
	}
   	results = list(Call = cl, Est = est, VarEst = varest, TotDomEst = totdom, MeanDomEst = meandom, Param = par, ConfInt = interv)
   	class(results) = "EstimatorDF"
   	attr(results, "attributesDF") = conf_level
   	return(results)			
}