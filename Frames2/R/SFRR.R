#' @name SFRR
#' @aliases SFRR
#' @title Raking ratio estimator
#' 
#' @description Produces estimates for population total and mean using the raking ratio estimator from survey data obtained
#'  from a dual frame sampling desing. Confidence intervals are also computed, if required.
#' 
#' @usage SFRR(ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, N_A, N_B, 
#' conf_level = NULL)
#' @param ysA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{c} containing information about variable of interest from \eqn{s_A}.
#' @param ysB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{c} containing information about variable of interest from \eqn{s_B}.
#' @param pi_A A numeric vector of length \eqn{n_A} or a square numeric matrix of dimension \eqn{n_A} containing first order or first and second order inclusion probabilities for units included in \eqn{s_A}.
#' @param pi_B A numeric vector of length \eqn{n_B} or a square numeric matrix of dimension \eqn{n_B} containing first order or first and second order inclusion probabilities for units included in \eqn{s_B}.
#' @param pik_ab_B A numeric vector of size \eqn{n_A} containing first order inclusion probabilities according to sampling desing in frame B for units belonging 
#'  to overlap domain that have been selected in \eqn{s_A}.
#' @param pik_ba_A A numeric vector of size \eqn{n_B} containing first order inclusion probabilities according to sampling desing in frame A for units belonging 
#'  to overlap domain that have been selected in \eqn{s_A}.
#' @param domains_A A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domains_B A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "b" and "ba".
#' @param N_A A numeric value indicating the size of frame A
#' @param N_B A numeric value indicating the size of frame B
#' @param conf_level (Optional) A numeric value indicating the confidence level for the confidence intervals, if desired.
#' @details Raking ratio estimator of population total is given by
#'  \deqn{\hat{Y}_{SFRR} = \frac{N_A - \hat{N}_{ab,rake}}{\hat{N}_a^A}\hat{Y}_a^A + \frac{N_B - \hat{N}_{ab,rake}}{\hat{N}_b^B}\hat{Y}_b^B + \frac{\hat{N}_{ab,rake}}{\hat{N}_{abS}}\hat{Y}_{abS}}
#'  where \eqn{\hat{Y}_{abS} = \sum_{i \in s_{ab}^A}\tilde{d}_i^Ay_i + \sum_{i \in s_{ab}^B}\tilde{d}_i^By_i, \hat{N}_{abS} = \sum_{i \in s_{ab}^A}\tilde{d}_i^A + \sum_{i \in s_{ab}^B}\tilde{d}_i^B} and
#'  \eqn{\hat{N}_{ab,rake}} is the smallest root of the quadratic equation \eqn{\hat{N}_{ab,rake}x^2 - [\hat{N}_{ab,rake}(N_A + N_B) + \hat{N}_{aS}\hat{N}_{bS}]x + \hat{N}_{ab,rake}N_AN_B = 0},
#'  with \eqn{\hat{N}_{aS} = \sum_{s_a^A}\tilde{d}_i^B} and \eqn{\hat{N}_{bS} = \sum_{s_b^B}\tilde{d}_i^B}. Weights \eqn{\tilde{d}_i^A} and \eqn{\tilde{d}_i^B} are obtained as follows
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
#'  being \eqn{d_i^A} and \eqn{d_i^B} the design weights, obtained as the inverse of the first order inclusion probabilities, that is \eqn{d_i^A = 1/\pi_i^A} and \eqn{d_i^B = 1/\pi_i^B}.
#'  
#'  To obtain an estimator of the variance for this estimator, one has taken into account that raking ratio estimator coincides with SF calibration estimator when frame sizes are known and "raking"
#'  method is used. So, one can use here Deville's expression to calculate an estimator for the variance of the raking ratio estimator
#'  \deqn{\hat{V}(\hat{Y}_{SFRR}) = \frac{1}{1-\sum_{k\in s} a_k^2}\sum_{k\in s}(1-\pi_k)\left(\frac{e_k}{\pi_k} - \sum_{l\in s} a_{l} \frac{e_l}{\pi_l}\right)^2}
#'  where \eqn{a_k=(1-\pi_k)/\sum_{l\in s} (1-\pi_l)} and \eqn{e_k} are the residuals of the regression with auxiliary variables as regressors.
#' @return \code{SFRR} returns an object of class "EstimatorDF" which is a list with, at least, the following components:
#'  \item{Call}{the matched call.}
#'  \item{Est}{total and mean estimation for main variable(s).}
#'  \item{VarEst}{variance estimation for main variable(s).}
#'  If parameter \code{conf_level} is different from \code{NULL}, object includes component
#'  \item{ConfInt}{total and mean estimation and confidence intervals for main variables(s).}
#'  In addition, components \code{TotDomEst} and \code{MeanDomEst} are available when estimator is based on estimators of the domains. Component \code{Param} shows value of parameters involded in calculation of the estimator (if any).
#'  By default, only \code{Est} component (or \code{ConfInt} component, if parameter \code{conf_level} is different from \code{NULL}) is shown. It is possible to access to all the components of the objects by using function \code{summary}.
#' @references Lohr, S. and Rao, J.N.K. (2000).
#' \emph{Inference in Dual Frame Surveys}. Journal of the American Statistical Association, Vol. 95, 271 - 280.
#' @references Rao, J.N.K. and Skinner, C.J. (1996).
#' \emph{Estimation in Dual Frame Surveys with Complex Designs}. Proceedings of the Survey Method Section, Statistical Society of Canada, 63 - 68.
#' @references Skinner, C.J. and Rao J.N.K. (1996).
#' \emph{Estimation in Dual Frame Surveys with Complex Designs}. Journal of the American Statistical Association, Vol. 91, 443, 349 - 356.
#' @references Skinner, C.J. (1991).
#' \emph{On the Efficiency of Raking Ratio Estimation for Multiple Frame Surveys}. Journal of the American Statistical Association, Vol. 86, 779 - 784.
#' @seealso \code{\link{JackSFRR}}
#' @examples
#' data(DatA)
#' data(DatB)
#' data(PiklA)
#' data(PiklB)
#' 
#' #Let calculate raking ratio estimator for population total for variable Clothing
#' SFRR(DatA$Clo, DatB$Clo, PiklA, PiklB, DatA$ProbB, DatB$ProbA, DatA$Domain, 
#' DatB$Domain, 1735, 1191)
#' 
#' #Now, let calculate raking ratio estimator and a 90% confidence interval for 
#' #population total for variable Feeding, considering only first order inclusion probabilities
#' SFRR(DatA$Feed, DatB$Feed, DatA$ProbA, DatB$ProbB, DatA$ProbB, DatB$ProbA, 
#' DatA$Domain, DatB$Domain, 1735, 1191, 0.90)
#' @export
SFRR = function (ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, N_A, N_B, conf_level = NULL)
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
	if (ncol(ysA) != ncol(ysB))
		stop("Number of variables does not match.")
	if (nrow(ysA) != nrow(pi_A) | nrow(ysA) != length(domains_A) | length(domains_A) != nrow(pi_A) | nrow(ysA) != length(pik_ab_B))
		stop("Arguments from frame A have different sizes.")
	if (nrow(ysB) != nrow(pi_B) | nrow(ysB) != length(domains_B) | length(domains_B) != nrow(pi_B) | nrow(ysB) != length(pik_ba_A))
		stop("Arguments from frame B have different sizes.")
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

	ones_a_A <- Domains (rep (1, n_A), domains_A, "a")
	ones_b_B <- Domains (rep (1, n_B), domains_B, "b")
	ones_ab_A <- Domains (rep (1, n_A), domains_A, "ab")
	ones_ab_B <- Domains (rep (1, n_B), domains_B, "ba")

	est <- matrix(, 2, c, dimnames = list(c("Total", "Mean"), cnames))
	varest <- matrix(, 2, c, dimnames = list(c("Var. Total", "Var. Mean"), cnames))
	totdom <- matrix(, 3, c, dimnames = list(c("Total dom. a", "Total abS", "Total dom. b"), cnames))
	meandom <- matrix(, 3, c, dimnames = list(c("Mean dom. a", "Mean abS", "Mean dom. b"), cnames))
	par <- 	NULL
	if (is.null(conf_level))
		interv <- NULL
	else
		interv <- matrix(, 6, c, dimnames = list(c("Total", "Lower Bound", "Upper Bound", "Mean", "Lower Bound", "Upper Bound"), cnames))


	if (!is.null(dim(drop(pi_A))) & !is.null(dim(drop(pi_B)))) {

		if (nrow(pi_A) != ncol(pi_A))
			stop("Pikl from frame A is not a square matrix.")
		if (nrow(pi_B) != ncol(pi_B))
			stop("Pikl from frame B is not a square matrix.")

		Nhat_a_A <- HT (ones_a_A, diag(pi_A))
		Nhat_b_B <- HT (ones_b_B, diag(pi_B))

		wi_A <- (1 / diag(pi_A)) * (domains_A == "a") + (1 / (diag(pi_A) + pik_ab_B)) * (domains_A == "ab")
		wi_B <- (1 / diag(pi_B)) * (domains_B == "b") + (1 / (diag(pi_B) + pik_ba_A)) * (domains_B == "ba")

		Nhat_aS_A <- sum(wi_A * ones_a_A)
		Nhat_bS_B <- sum(wi_B * ones_b_B)
		Nhat_abS <- sum (wi_A * ones_ab_A) + sum(wi_B * ones_ab_B)

		term_a <- Nhat_abS
		term_b <- -(Nhat_abS * (N_A + N_B) + Nhat_aS_A * Nhat_bS_B)
		term_c <- Nhat_abS * N_A * N_B

		Nhat_ab_rake <- (- term_b - sqrt(term_b * term_b - 4 * term_a * term_c)) / (2 * term_a)
 
		for (k in 1:(c+1)) {

			data_a_A <- Domains (ysA[,k], domains_A, "a")
			data_b_B <- Domains (ysB[,k], domains_B, "b")
			Yhat_a_A <- HT (data_a_A, diag(pi_A))
			Yhat_b_B <- HT (data_b_B, diag(pi_B))
			Yhat_abS <- sum (wi_A * ysA[,k] * ones_ab_A) + sum(wi_B * ysB[,k] * ones_ab_B)

			if (k == 1){
				domain_size_estimation <- c(Yhat_a_A, Yhat_abS, Yhat_b_B)
				size_estimation <- sum((N_A - Nhat_ab_rake) * Yhat_a_A / Nhat_a_A, (N_B - Nhat_ab_rake) * Yhat_b_B / Nhat_b_B, Nhat_ab_rake / Nhat_abS * Yhat_abS, na.rm = TRUE)
			}
			else
				total_estimation <- sum((N_A - Nhat_ab_rake) * Yhat_a_A / Nhat_a_A, (N_B - Nhat_ab_rake) * Yhat_b_B / Nhat_b_B, Nhat_ab_rake / Nhat_abS * Yhat_abS, na.rm = TRUE) 
			
			if (k > 1) {

				totdom[,k-1] <- c(Yhat_a_A, Yhat_abS, Yhat_b_B)
				meandom[,k-1] <- totdom[,k-1]/domain_size_estimation

				mean_estimation <- total_estimation / size_estimation
				est[,k-1] <- c(total_estimation, mean_estimation)

				d <- c(wi_A, wi_B)
				n <- n_A + n_B
				sample <- c(ysA[,k], ysB[,k])
				domains <- factor(c(as.character(domains_A), as.character(domains_B)))
				delta_a <- Domains (rep (1, n), domains, "a")
				delta_ab <- Domains (rep (1, n), domains, "ab")
				delta_b <- Domains (rep (1, n), domains, "b")
				delta_ba <- Domains (rep (1, n), domains, "ba")
				Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b)
				total <- c(N_A - Nhat_ab_rake, Nhat_ab_rake, N_B - Nhat_ab_rake)
				g <- calib (Xs, d, total, method = "raking")

				Vhat_Yhat_SFRR <- varest(sample, Xs, 1/d, g)
				Vhat_Ymeanhat_SFRR <- 1/size_estimation^2 * Vhat_Yhat_SFRR
				varest[,k-1] <- c(Vhat_Yhat_SFRR, Vhat_Ymeanhat_SFRR)

				if (!is.null(conf_level)) {

					total_upper <- total_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_SFRR)
					total_lower <- total_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_SFRR)
					mean_upper <- mean_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_SFRR)
					mean_lower <- mean_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_SFRR)
					interv[,k-1] <- c(total_estimation, total_lower, total_upper, mean_estimation, mean_lower, mean_upper)
		
				}
			}
		}
	}
	else {

		if (is.null(dim(drop(pi_A))) & is.null(dim(drop(pi_B)))){

			Nhat_a_A <- HT (ones_a_A, pi_A)
			Nhat_b_B <- HT (ones_b_B, pi_B)

			wi_A <- (1 / pi_A) * (domains_A == "a") + (1 / (pi_A + pik_ab_B)) * (domains_A == "ab")
			wi_B <- (1 / pi_B) * (domains_B == "b") + (1 / (pi_B + pik_ba_A)) * (domains_B == "ba")

			Nhat_aS_A <- sum(wi_A * ones_a_A)
			Nhat_bS_B <- sum(wi_B * ones_b_B)
			Nhat_abS <- sum (wi_A * ones_ab_A) + sum(wi_B * ones_ab_B)

			term_a <- Nhat_abS
			term_b <- -(Nhat_abS * (N_A + N_B) + Nhat_aS_A * Nhat_bS_B)
			term_c <- Nhat_abS * N_A * N_B

			Nhat_ab_rake <- (- term_b - sqrt(term_b * term_b - 4 * term_a * term_c)) / (2 * term_a)
 
			for (k in 1:(c+1)) {

				data_a_A <- Domains (ysA[,k], domains_A, "a")
				data_b_B <- Domains (ysB[,k], domains_B, "b")
				Yhat_a_A <- HT (data_a_A, pi_A)
				Yhat_b_B <- HT (data_b_B, pi_B)
				Yhat_abS <- sum (wi_A * ysA[,k] * ones_ab_A) + sum(wi_B * ysB[,k] * ones_ab_B)

				if (k == 1){
					domain_size_estimation <- c(Yhat_a_A, Yhat_abS, Yhat_b_B)
					size_estimation <- sum((N_A - Nhat_ab_rake) * Yhat_a_A / Nhat_a_A, (N_B - Nhat_ab_rake) * Yhat_b_B / Nhat_b_B, Nhat_ab_rake / Nhat_abS * Yhat_abS, na.rm = TRUE)
				}
				else
					total_estimation <- sum((N_A - Nhat_ab_rake) * Yhat_a_A / Nhat_a_A, (N_B - Nhat_ab_rake) * Yhat_b_B / Nhat_b_B, Nhat_ab_rake / Nhat_abS * Yhat_abS, na.rm = TRUE)
				
				if (k > 1) {

					totdom[,k-1] <- c(Yhat_a_A, Yhat_abS, Yhat_b_B)
					meandom[,k-1] <- totdom[,k-1]/domain_size_estimation

					mean_estimation <- total_estimation / size_estimation
					est[,k-1] <- c(total_estimation, mean_estimation)

					d <- c(wi_A, wi_B)
					n <- n_A + n_B
					sample <- c(ysA[,k], ysB[,k])
					domains <- factor(c(as.character(domains_A), as.character(domains_B)))
					delta_a <- Domains (rep (1, n), domains, "a")
					delta_ab <- Domains (rep (1, n), domains, "ab")
					delta_b <- Domains (rep (1, n), domains, "b")
					delta_ba <- Domains (rep (1, n), domains, "ba")
					Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b)
					total <- c(N_A - Nhat_ab_rake, Nhat_ab_rake, N_B - Nhat_ab_rake)
					g <- calib (Xs, d, total, method = "raking")

					Vhat_Yhat_SFRR <- varest(sample, Xs, 1/d, g)
					Vhat_Ymeanhat_SFRR <- 1/size_estimation^2 * Vhat_Yhat_SFRR
					varest[,k-1] <- c(Vhat_Yhat_SFRR, Vhat_Ymeanhat_SFRR)

					if (!is.null(conf_level)) {

						total_upper <- total_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_SFRR)
						total_lower <- total_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_SFRR)
						mean_upper <- mean_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_SFRR)
						mean_lower <- mean_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_SFRR)
						interv[,k-1] <- c(total_estimation, total_lower, total_upper, mean_estimation, mean_lower, mean_upper)	
					}
				}
			}
		}
		else
			stop("Invalid option: Probability vector in one frame and probability matrix in the other frame. Type of both structures must match.")
		
	}
   	results = list(Call = cl, Est = est, VarEst = varest, TotDomEst = totdom, MeanDomEst = meandom, Param = par, ConfInt = interv)
   	class(results) = "EstimatorDF"
   	attr(results, "attributesDF") = conf_level
   	return(results)			
}