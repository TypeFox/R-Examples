#' @name CalSF
#' @aliases CalSF
#' @title SF calibration estimator
#' 
#' @description Produces estimates for population totals and means using the SF calibration estimator from survey data obtained
#'  from a dual frame sampling design. Confidence intervals are also computed, if required. 
#' 
#' @usage CalSF(ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, N_A = NULL,
#' N_B = NULL, N_ab = NULL, xsAFrameA = NULL, xsBFrameA = NULL, xsAFrameB = NULL, 
#' xsBFrameB = NULL, xsT = NULL, XA = NULL, XB = NULL, X = NULL, met = "linear", 
#' conf_level = NULL)
#' @param ysA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{c} containing information about variable(s) of interest from \eqn{s_A}.
#' @param ysB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{c} containing information about variable(s) of interest from \eqn{s_B}.
#' @param pi_A A numeric vector of length \eqn{n_A} or a square numeric matrix of dimension \eqn{n_A} containing first order or first and second order inclusion probabilities for units included in \eqn{s_A}.
#' @param pi_B A numeric vector of length \eqn{n_B} or a square numeric matrix of dimension \eqn{n_B} containing first order or first and second order inclusion probabilities for units included in \eqn{s_B}.
#' @param pik_ab_B A numeric vector of size \eqn{n_A} containing first order inclusion probabilities according to sampling desing in frame B for units belonging 
#'  to overlap domain that have been selected in \eqn{s_A}.
#' @param pik_ba_A A numeric vector of size \eqn{n_B} containing first order inclusion probabilities according to sampling desing in frame A for units belonging 
#'  to overlap domain that have been selected in \eqn{s_B}.
#' @param domains_A A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domains_B A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".
#' @param N_A (Optional) A numeric value indicating the size of frame A
#' @param N_B (Optional) A numeric value indicating the size of frame B
#' @param N_ab (Optional) A numeric value indicating the size of the overlap domain
#' @param xsAFrameA (Optional) A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, containing auxiliary information in frame A for units included in \eqn{s_A}.
#' @param xsBFrameA (Optional) A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, containing auxiliary information in frame A for units included in \eqn{s_B}. For units in domain \eqn{b}, these values are 0.
#' @param xsAFrameB (Optional) A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, containing auxiliary information in frame B for units included in \eqn{s_A}. For units in domain \eqn{a}, these values are 0.
#' @param xsBFrameB (Optional) A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, containing auxiliary information in frame B for units included in \eqn{s_B}.
#' @param xsT (Optional) A numeric vector of length \eqn{n} or a numeric matrix or data frame of dimensions \eqn{n} x \eqn{m_T}, with \eqn{m_T} the number of auxiliary variables in both frames, containing auxiliary information for all units in the entire sample \eqn{s = s_A \cup s_B}.
#' @param XA (Optional) A numeric value or vector of length \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, indicating the population totals for the auxiliary variables considered in frame A.
#' @param XB (Optional) A numeric value or vector of length \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, indicating the population totals for the auxiliary variables considered in frame B.
#' @param X (Optional) A numeric value or vector of length \eqn{m_T}, with \eqn{m_T} the number of auxiliary variables in both frames, indicating the population totals for the auxiliary variables considered in both frames.
#' @param met (Optional) A character vector indicating the distance that must be used in calibration process. Possible values are "linear", "raking" and "logit". Default is "linear".
#' @param conf_level (Optional) A numeric value indicating the confidence level for the confidence intervals, if desired.
#' @details SF calibration estimator of population total is given by
#'  \deqn{\hat{Y}_{CalSF} = \hat{Y}_a + \hat{Y}_{ab} + \hat{Y}_b}
#'  where \eqn{\hat{Y}_a = \sum_{i \in s_a}\tilde{d}_i y_i, \hat{Y}_{ab} = \sum_{i \in (s_{ab} \cup s_{ba})}\tilde{d}_i y_i} 
#'  and \eqn{\hat{Y}_b = \sum_{i \in s_b} \tilde{d}_i y_i}, with \eqn{\tilde{d}_i} calibration weights which are calculated
#'  having into account a different set of constraints, depending on the case. For instance, if \eqn{N_A, N_B} and \eqn{N_{ab}} are known and no other auxiliary information is available, calibration constraints are
#'  \deqn{\sum_{i \in s_a}\tilde{d}_i = N_a, \sum_{i \in s_{ab} \cup s_{ba}}\tilde{d}_i = N_{ab}, \sum_{i \in s_{ba}}\tilde{d}_i = N_{ba}}
#'
#'  Function covers following scenarios:
#'  \itemize{
#'	\item There is not any additional auxiliary variable 
#'      \itemize{
#'      	\item \eqn{N_A, N_B} and \eqn{N_{ab}} unknown
#'		\item \eqn{N_{ab}} known and \eqn{N_A} and \eqn{N_B} unknown 
#'              \item \eqn{N_A} and \eqn{N_B} known and \eqn{N_{ab}} unknown
#'              \item \eqn{N_A, N_B} and \eqn{N_{ab}} known
#'      }
#'      \item At least, information about one additional auxiliary variable is available 
#'      \itemize{
#'		\item \eqn{N_{ab}} known and \eqn{N_A} and \eqn{N_B} unknown
#'              \item \eqn{N_A} and \eqn{N_B} known and \eqn{N_{ab}} unknown
#'              \item \eqn{N_A, N_B} and \eqn{N_{ab}} known
#'      }
#'  }
#'
#'  To obtain an estimator of the variance for this estimator, one can use Deville's expression
#'  \deqn{\hat{V}(\hat{Y}_{CalSF}) = \frac{1}{1-\sum_{k\in s} a_k^2}\sum_{k\in s}(1-\pi_k)\left(\frac{e_k}{\pi_k} - \sum_{l\in s} a_{l} \frac{e_l}{\pi_l}\right)^2}
#'  where \eqn{a_k=(1-\pi_k)/\sum_{l\in s} (1-\pi_l)} and \eqn{e_k} are the residuals of the regression with auxiliary variables as regressors.
#' @return \code{CalSF} returns an object of class "EstimatorDF" which is a list with, at least, the following components:
#'  \item{Call}{the matched call.}
#'  \item{Est}{total and mean estimation for main variable(s).}
#'  \item{VarEst}{variance estimation for main variable(s).}
#'  If parameter \code{conf_level} is different from \code{NULL}, object includes component
#'  \item{ConfInt}{total and mean estimation and confidence intervals for main variables(s).}
#'  In addition, components \code{TotDomEst} and \code{MeanDomEst} are available when estimator is based on estimators of the domains. Component \code{Param} shows value of parameters involded in calculation of the estimator (if any).
#'  By default, only \code{Est} component (or \code{ConfInt} component, if parameter \code{conf_level} is different from \code{NULL}) is shown. It is possible to access to all the components of the objects by using function \code{summary}.
#' @references Ranalli, M. G., Arcos, A., Rueda, M. and Teodoro, A. (2013)
#'  \emph{Calibration estimation in dual frame surveys}. arXiv:1312.0761 [stat.ME]
#' @references Deville, J. C., Sarndal, C. E. (1992)
#'  \emph{Calibration estimators in survey sampling.}
#'  Journal of the American Statistical Association, 87, 376 - 382
#' @seealso \code{\link{JackCalSF}}
#' @examples
#' data(DatA)
#' data(DatB)
#' data(PiklA)
#' data(PiklB)
#' 
#' #Let calculate SF calibration estimator for variable Clothing, without
#' #considering any auxiliary information
#' CalSF(DatA$Clo, DatB$Clo, PiklA, PiklB, DatA$ProbB, DatB$ProbA, 
#' DatA$Domain, DatB$Domain)
#' 
#' #Now, let calculate SF calibration estimator for variable Leisure when the frame
#' #sizes and the overlap domain size are known
#' CalSF(DatA$Lei, DatB$Lei, PiklA, PiklB, DatA$ProbB, DatB$ProbA, DatA$Domain, 
#' DatB$Domain, N_A = 1735, N_B = 1191, N_ab = 601)
#' 
#' #Finally, let calculate SF calibration estimator and a 90% confidence interval
#' #for population total for variable Feeding, considering Income and Metres2 as auxiliary 
#' #variables and with frame sizes and overlap domain size known.
#' CalSF(DatA$Feed, DatB$Feed, PiklA, PiklB, DatA$ProbB, DatB$ProbA, DatA$Domain, 
#' DatB$Domain, N_A = 1735, N_B =  1191, N_ab = 601, xsAFrameA = DatA$Inc, 
#' xsBFrameA = DatB$Inc, xsAFrameB = DatA$M2, xsBFrameB = DatB$M2, 
#' XA = 4300260, XB = 176553, conf_level = 0.90)
#' @export
CalSF = function (ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, N_A = NULL, N_B = NULL, N_ab = NULL, xsAFrameA = NULL, xsBFrameA = NULL, xsAFrameB = NULL, xsBFrameB = NULL, xsT = NULL, XA = NULL, XB = NULL, X = NULL, met = "linear", conf_level = NULL)
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
	if (nrow(ysA) != nrow(pi_A) | nrow(ysA) != length(domains_A) | length(domains_A) != nrow(pi_A))
		stop("Arguments from frame A have different sizes.")
	if (nrow(ysB) != nrow(pi_B) | nrow(ysB) != length(domains_B) | length(domains_B) != nrow(pi_B))
		stop("Arguments from frame B have different sizes.")
	if (length(which(domains_A == "a")) + length(which(domains_A == "ab")) != length(domains_A))
		stop("Domains from frame A are not correct.")
	if (length(which(domains_B == "b")) + length(which(domains_B == "ba")) != length(domains_B))
		stop("Domains from frame B are not correct.")
	if ((is.null (N_A) & !is.null (N_B)) | (!is.null (N_A) & is.null (N_B)))
		stop("Only one value has been indicated for N_A and N_B. This is not valid.") 
	if ((is.null (xsAFrameA) & !is.null (xsBFrameA)) | (!is.null (xsAFrameA) & is.null (xsBFrameA)))
		stop("Auxiliary information from frame A is available only in one frame. This is not a possible option.")
	if ((is.null (xsAFrameB) & !is.null (xsBFrameB)) | (!is.null (xsAFrameB) & is.null (xsBFrameB)))
		stop("Auxiliary information from frame B is available only in one frame. This is not a possible option.")

	cl <- match.call()

	sample <- rbind(ysA, ysB)
	n_A <- nrow(ysA)
 	n_B <- nrow(ysB)
	n <- n_A + n_B
	c <- ncol(ysA)

	domains <- factor(c(as.character(domains_A), as.character(domains_B)))  
	ysA <- cbind(rep(1, n_A), ysA)
	ysB <- cbind(rep(1, n_B), ysB)
	sample <- rbind(ysA, ysB)

	delta_a <- Domains (rep (1, n), domains, "a")
	delta_ab <- Domains (rep (1, n), domains, "ab")
	delta_b <- Domains (rep (1, n), domains, "b")
	delta_ba <- Domains (rep (1, n), domains, "ba")

	est <- matrix(, 2, c, dimnames = list(c("Total", "Mean"), cnames))
	varest <- matrix(, 2, c, dimnames = list(c("Var. Total", "Var. Mean"), cnames))
	totdom <- NULL
	meandom <- NULL
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
		
		pik_A <- diag(pi_A); pik_B <- diag(pi_B)
		w_tilde_iS_A <- (1 / diag(pi_A)) * (domains_A == "a") + (1 / (diag(pi_A) + pik_ab_B)) * (domains_A == "ab")
		w_tilde_iS_B <- (1 / diag(pi_B)) * (domains_B == "b") + (1 / (diag(pi_B) + pik_ba_A)) * (domains_B == "ba")
		d <- c(w_tilde_iS_A, w_tilde_iS_B)

		for (k in 1:(c+1)) {

			if (is.null(xsAFrameA) & is.null(xsBFrameB) & is.null(xsT)) {

				if (is.null(N_ab)) {
				
					if (is.null(N_A) & is.null(N_B)) {

						Nhat_abS <- sum((1 / (pik_A + pik_ab_B)) * (domains_A == "ab")) + sum((1 / (pik_B + pik_ba_A)) * (domains_B == "ba"))

						Nhat_A <- HT (rep(1, n_A), pik_A)
						Nhat_B <- HT (rep(1, n_B), pik_B)					

						Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b)
						total <- c(Nhat_A - Nhat_abS, Nhat_abS, Nhat_B - Nhat_abS)
					}
					else {
	
						Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba)
						total <- c(N_A, N_B)
					}
				}
				else {

					if (is.null(N_A) & is.null(N_B)) {

						Xs <- cbind(delta_ab + delta_ba)
						total <- c(N_ab)
					}
					else{

						Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b)
						total <- c(N_A - N_ab, N_ab, N_B - N_ab)
					}
				}
			}
			else {

				if (is.null(N_ab)) {

					if (is.null(xsAFrameA)){

						if (is.null(xsBFrameB)){

							xsT <- as.matrix(xsT)
							Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, xsT)
							total <- c(N_A, N_B, X)

						}
						else{

							xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
							XFrameB <- rbind(xsAFrameB, xsBFrameB)

							if (is.null(xsT)){

								Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, (delta_b + delta_ab + delta_ba) * XFrameB)
								total <- c(N_A, N_B, XB)
							}
							else {

								xsT <- as.matrix(xsT)
								Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + + delta_ab + delta_ba, (delta_b + delta_ab + delta_ba) * XFrameB, xsT)
								total <- c(N_A, N_B, XB, X)
							}
						}
					}
					else {

						xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
						XFrameA <- rbind(xsAFrameA, xsBFrameA)

						if (is.null(xsBFrameB)){
						
							if (is.null(xsT)){

								Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA)
								total <- c(N_A, N_B, XA) 
							}
							else{

								xsT <- as.matrix(xsT)
								Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + + delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA, xsT)
								total <- c(N_A, N_B, XA, X) 
							}
						}
						else{

							xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
							XFrameB <- rbind(xsAFrameB, xsBFrameB)

							if (is.null(xsT)){

								Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA, (delta_b + delta_ab + delta_ba) * XFrameB)
								total <- c(N_A, N_B, XA, XB) 
							}
							else{

								xsT <- as.matrix(xsT)
								Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA, (delta_b + delta_ab + delta_ba) * XFrameB, xsT)
								total <- c(N_A, N_B, XA, XB, X) 
							}
						}

					}
				
				}
				else {

					if (is.null(N_A) & is.null(N_B)) {

						if (is.null(xsAFrameA)){

							if (is.null(xsBFrameB)){

								Xs <- cbind(delta_ab + delta_ba, xsT)
								total <- c(N_ab, X)
							}
							else{

								xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
								XFrameB <- rbind(xsAFrameB, xsBFrameB)

								if (is.null(xsT)){

									Xs <- cbind(delta_ab + delta_ba, (delta_b + delta_ab + delta_ba) * XFrameB)
									total <- c(N_ab, XB)
								}
								else {
									Xs <- cbind(delta_ab + delta_ba, (delta_b + delta_ab + delta_ba) * XFrameB, xsT)
									total <- c(N_ab, XB, X)
								}
							}
						}
						else {

							xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
							XFrameA <- rbind(xsAFrameA, xsBFrameA)

							if (is.null(xsBFrameB)){
						
								if (is.null(xsT)){

									Xs <- cbind(delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA)
									total <- c(N_ab, XA) 
								}

								else{
									Xs <- cbind(delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA, xsT)
									total <- c(N_ab, XA, X) 
								}
							}
							else{

								xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
								XFrameB <- rbind(xsAFrameB, xsBFrameB)

								if (is.null(xsT)){

									Xs <- cbind(delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA, (delta_b + delta_ab + delta_ba) * XFrameB)                
									total <- c(N_ab, XA, XB)
								}
								else{
									Xs <- cbind(delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA, (delta_b + delta_ab + delta_ba) * XFrameB, xsT)
									total <- c(N_ab, XA, XB, X) 
								}
							}

						}
					}
					else{

						if (is.null(xsAFrameA)){

							if (is.null(xsBFrameB)){

								Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, xsT)
								total <- c(N_A - N_ab, N_ab, N_B - N_ab, X)
							}
							else{

								xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
								XFrameB <- rbind(xsAFrameB, xsBFrameB)

								if (is.null(xsT)){

									Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, (delta_b + delta_ab + delta_ba) * XFrameB)
									total <- c(N_A - N_ab, N_ab, N_B - N_ab, XB)
								}
								else {
									Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, (delta_b + delta_ab + delta_ba) * XFrameB, xsT)
									total <- c(N_A - N_ab, N_ab, N_B - N_ab, XB, X)
								}
							}
						}
						else {

							xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
							XFrameA <- rbind(xsAFrameA, xsBFrameA)

							if (is.null(xsBFrameB)){
						
								if (is.null(xsT)){

									Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, (delta_a + delta_ab + delta_ba) * XFrameA)
									total <- c(N_A - N_ab, N_ab, N_B - N_ab, XA) 
								}
	
								else{
									Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, (delta_a + delta_ab + delta_ba) * XFrameA, xsT)
									total <- c(N_A - N_ab, N_ab, N_B - N_ab, XA, X) 
								}
							}
							else{

								xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
								XFrameB <- rbind(xsAFrameB, xsBFrameB)

								if (is.null(xsT)){

									Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, (delta_a + delta_ab + delta_ba) * XFrameA, (delta_b + delta_ab + delta_ba) * XFrameB)                
									total <- c(N_A - N_ab, N_ab, N_B - N_ab, XA, XB)
								}
								else{

									Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, (delta_a + delta_ab + delta_ba) * XFrameA, (delta_b + delta_ab + delta_ba) * XFrameB, xsT)
									total <- c(N_A - N_ab, N_ab, N_B - N_ab, XA, XB, X) 
								}
							}

						}
					}
				}
			}
			g <- calib (Xs, d, total, method = met)

			if (k == 1)
				size_estimation <- sum (g * d * sample[,k])
			else
				total_estimation <- sum (g * d * sample[,k])

			if (k > 1) {

				mean_estimation <- total_estimation / size_estimation
				est[,k-1] <- c(total_estimation, mean_estimation)
				Vhat_Yhat_CalSF <- varest(sample[,k], Xs, 1/d, g)
				Vhat_Ymeanhat_CalSF <- 1/size_estimation^2 * Vhat_Yhat_CalSF
				varest[,k-1] <- c(Vhat_Yhat_CalSF, Vhat_Ymeanhat_CalSF)

				if (!is.null(conf_level)) {

					total_upper <- total_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_CalSF)
					total_lower <- total_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_CalSF)
					mean_upper <- mean_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_CalSF)
					mean_lower <- mean_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_CalSF)
					interv[,k-1] <- c(total_estimation, total_lower, total_upper, mean_estimation, mean_lower, mean_upper)	
				}
			}
		}
	}
	else {

		if (is.null(dim(drop(pi_A))) & is.null(dim(drop(pi_B)))){

			w_tilde_iS_A <- (1 / pi_A) * (domains_A == "a") + (1 / (pi_A + pik_ab_B)) * (domains_A == "ab")
			w_tilde_iS_B <- (1 / pi_B) * (domains_B == "b") + (1 / (pi_B + pik_ba_A)) * (domains_B == "ba")
			d <- c(w_tilde_iS_A, w_tilde_iS_B)

			for (k in 1:(c+1)) {

				if (is.null(xsAFrameA) & is.null(xsBFrameB) & is.null(xsT)) {

					if (is.null(N_ab)) {
				
						if (is.null(N_A) & is.null(N_B)) {

							Nhat_abS <- sum((1 / (pi_A + pik_ab_B)) * (domains_A == "ab")) + sum((1 / (pi_B + pik_ba_A)) * (domains_B == "ba"))

							Nhat_A <- HT (rep(1, n_A), pi_A)
							Nhat_B <- HT (rep(1, n_B), pi_B)
							Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b)
							total <- c(Nhat_A - Nhat_abS, Nhat_abS, Nhat_B - Nhat_abS)
						}
						else {
	
							Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba)
							total <- c(N_A, N_B)
						}
					}
					else {

						if (is.null(N_A) & is.null(N_B)) {
				
							Xs <- cbind(delta_ab + delta_ba)
							total <- c(N_ab)
						}
						else{

							Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b)
							total <- c(N_A - N_ab, N_ab, N_B - N_ab)
						}
					}
				}
				else {

					if (is.null(N_ab)) {

						if (is.null(xsAFrameA)){

							if (is.null(xsBFrameB)){

								xsT <- as.matrix(xsT)
								Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, xsT)
								total <- c(N_A, N_B, X)
							}
							else{

								xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
								XFrameB <- rbind(xsAFrameB, xsBFrameB)

								if (is.null(xsT)){

									Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, (delta_b + delta_ab + delta_ba) * XFrameB)
									total <- c(N_A, N_B, XB)
								}
								else {

									xsT <- as.matrix(xsT)
									Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, (delta_b + delta_ab + delta_ba) * XFrameB, xsT)
									total <- c(N_A, N_B, XB, X)
								}
							}
						}
						else {

							xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
							XFrameA <- rbind(xsAFrameA, xsBFrameA)

							if (is.null(xsBFrameB)){
						
								if (is.null(xsT)){

									Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA)
									total <- c(N_A, N_B, XA) 
								}
								else{

									xsT <- as.matrix(xsT)
									Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA, xsT)
									total <- c(N_A, N_B, XA, X) 
								}
							}
							else{

								xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
								XFrameB <- rbind(xsAFrameB, xsBFrameB)

								if (is.null(xsT)){

									Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA, (delta_b + delta_ab + delta_ba) * XFrameB)
									total <- c(N_A, N_B, XA, XB) 
								}
								else{

									xsT <- as.matrix(xsT)
									Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA, (delta_b + delta_ab + delta_ba) * XFrameB, xsT)
									total <- c(N_A, N_B, XA, XB, X) 
								}
							}

						}
				
					}
					else {

						if (is.null(N_A) & is.null(N_B)) {

							if (is.null(xsAFrameA)){

								if (is.null(xsBFrameB)){
									
									xsT <- as.matrix(xsT)
									Xs <- cbind(delta_ab + delta_ba, xsT)
									total <- c(N_ab, X)
								}
								else{

									xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
									XFrameB <- rbind(xsAFrameB, xsBFrameB)

									if (is.null(xsT)){

										Xs <- cbind(delta_ab + delta_ba, (delta_b + delta_ab + delta_ba) * XFrameB)
										total <- c(N_ab, XB)
									}
									else {

										Xs <- cbind(delta_ab + delta_ba, (delta_b + delta_ab + delta_ba) * XFrameB, xsT)
										total <- c(N_ab, XB, X)
									}
								}
							}
							else {

								xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
								XFrameA <- rbind(xsAFrameA, xsBFrameA)

								if (is.null(xsBFrameB)){
						
									if (is.null(xsT)){

										Xs <- cbind(delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA)
										total <- c(N_ab, XA) 
									}

									else{

										Xs <- cbind(delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA, xsT)
										total <- c(N_ab, XA, X) 
									}
								}
								else{

									xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
									XFrameB <- rbind(xsAFrameB, xsBFrameB)

									if (is.null(xsT)){

										Xs <- cbind(delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA, (delta_b + delta_ab + delta_ba) * XFrameB)                
										total <- c(N_ab, XA, XB)
									}
									else{

										Xs <- cbind(delta_ab + delta_ba, (delta_a + delta_ab + delta_ba) * XFrameA, (delta_b + delta_ab + delta_ba) * XFrameB, xsT)
										total <- c(N_ab, XA, XB, X) 
									}
								}

							}
						}
						else {

							if (is.null(xsAFrameA)){

								if (is.null(xsBFrameB)){

									xsT <- as.matrix(xsT)
									Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, xsT)
									total <- c(N_A - N_ab, N_ab, N_B - N_ab, X)
								}
								else {

									xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
									XFrameB <- rbind(xsAFrameB, xsBFrameB)

									if (is.null(xsT)){

										Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, (delta_b + delta_ab + delta_ba) * XFrameB)
										total <- c(N_A - N_ab, N_ab, N_B - N_ab, XB)
									}
									else {

										Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, (delta_b + delta_ab + delta_ba) * XFrameB, xsT)
										total <- c(N_A - N_ab, N_ab, N_B - N_ab, XB, X)
									}
								}
							}
							else {

								xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
								XFrameA <- rbind(xsAFrameA, xsBFrameA)

								if (is.null(xsBFrameB)){
						
									if (is.null(xsT)){

										Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, (delta_a + delta_ab + delta_ba) * XFrameA)
										total <- c(N_A - N_ab, N_ab, N_B - N_ab, XA) 
									}
									else {

										Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, (delta_a + delta_ab + delta_ba) * XFrameA, xsT)
										total <- c(N_A - N_ab, N_ab, N_B - N_ab, XA, X) 
									}
								}
								else {

									xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
									XFrameB <- rbind(xsAFrameB, xsBFrameB)

									if (is.null(xsT)){

										Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, (delta_a + delta_ab + delta_ba) * XFrameA, (delta_b + delta_ab + delta_ba) * XFrameB)                
										total <- c(N_A - N_ab, N_ab, N_B - N_ab, XA, XB)
									}
									else {

										Xs <- cbind(delta_a, delta_ab + delta_ba, delta_b, (delta_a + delta_ab + delta_ba) * XFrameA, (delta_b + delta_ab + delta_ba) * XFrameB, xsT)
										total <- c(N_A - N_ab, N_ab, N_B - N_ab, XA, XB, X) 
									}
								}
							}
						}
					}
				}
				g <- calib (Xs, d, total, method = met)

				if (k == 1)
					size_estimation <- sum (g * d * sample[,k])
				else
					total_estimation <- sum (g * d * sample[,k])

				if (k > 1) {

					mean_estimation <- total_estimation / size_estimation
					est[,k-1] <- c(total_estimation, mean_estimation)
					Vhat_Yhat_CalSF <- varest(sample[,k], Xs, 1/d, g)
					Vhat_Ymeanhat_CalSF <- 1/size_estimation^2 * Vhat_Yhat_CalSF
					varest[,k-1] <- c(Vhat_Yhat_CalSF, Vhat_Ymeanhat_CalSF)

					if (!is.null(conf_level)) {

						total_upper <- total_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_CalSF)
						total_lower <- total_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_CalSF)
						mean_upper <- mean_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_CalSF)
						mean_lower <- mean_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Ymeanhat_CalSF)
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