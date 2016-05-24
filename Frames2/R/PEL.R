#' @name PEL
#' @aliases PEL
#' @title Pseudo empirical likelihood estimator
#' 
#' @description Produces estimates for population totals using the pseudo empirical likelihood estimator from survey data obtained
#'  from a dual frame sampling design. Confidence intervals for the population total are also computed, if required.
#' 
#' @usage PEL(ysA, ysB, pi_A, pi_B, domains_A, domains_B, N_A = NULL, N_B = NULL, 
#' N_ab = NULL, xsAFrameA = NULL, xsBFrameA = NULL, xsAFrameB = NULL, xsBFrameB = NULL, 
#' XA = NULL, XB = NULL, conf_level = NULL)
#' @param ysA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{c} containing information about variable(s) of interest from \eqn{s_A}.
#' @param ysB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{c} containing information about variable(s) of interest from \eqn{s_B}.
#' @param pi_A A numeric vector of length \eqn{n_A} or a square numeric matrix of dimension \eqn{n_A} containing first order or first and second order inclusion probabilities for units included in \eqn{s_A}.
#' @param pi_B A numeric vector of length \eqn{n_B} or a square numeric matrix of dimension \eqn{n_B} containing first order or first and second order inclusion probabilities for units included in \eqn{s_B}.
#' @param domains_A A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domains_B A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".
#' @param N_A (Optional) A numeric value indicating the size of frame A. 
#' @param N_B (Optional) A numeric value indicating the size of frame B.
#' @param N_ab (Optional) A numeric value indicating the size of the overlap domain.
#' @param xsAFrameA (Optional) A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, containing auxiliary information in frame A for units included in \eqn{s_A}.
#' @param xsBFrameA (Optional) A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, containing auxiliary information in frame A for units included in \eqn{s_B}. For units in domain \eqn{b}, these values are 0.
#' @param xsAFrameB (Optional) A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, containing auxiliary information in frame B for units included in \eqn{s_A}. For units in domain \eqn{a}, these values are 0.
#' @param xsBFrameB (Optional) A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, containing auxiliary information in frame B for units included in \eqn{s_B}.
#' @param XA (Optional) A numeric value or vector of length \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, indicating the population totals for the auxiliary variables considered in frame A.
#' @param XB (Optional) A numeric value or vector of length \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, indicating the population totals for the auxiliary variables considered in frame B.
#' @param conf_level (Optional) A numeric value indicating the confidence level for the confidence intervals, if desired.
#' @return \code{PEL} returns an object of class "EstimatorDF" which is a list with, at least, the following components:
#'  \item{Call}{the matched call.}
#'  \item{Est}{total and mean estimation for main variable(s).}
#'  \item{VarEst}{variance estimation for main variable(s).}
#'  If parameter \code{conf_level} is different from \code{NULL}, object includes component
#'  \item{ConfInt}{total and mean estimation and confidence intervals for main variables(s).}
#'  In addition, components \code{TotDomEst} and \code{MeanDomEst} are available when estimator is based on estimators of the domains. Component \code{Param} shows value of parameters involded in calculation of the estimator (if any).
#'  By default, only \code{Est} component (or \code{ConfInt} component, if parameter \code{conf_level} is different from \code{NULL}) is shown. It is possible to access to all the components of the objects by using function \code{summary}.

#' @details Pseudo empirical likelihood estimator for the population mean is computed as
#'  \deqn{\hat{\bar{Y}}_{PEL} = \frac{N_a}{N}\hat{\bar{Y}}_a + \frac{\eta N_{ab}}{N}\hat{\bar{Y}}_{ab}^A + \frac{(1 - \eta) N_{ab}}{N}\hat{\bar{Y}}_{ab}^B + \frac{N_b}{N}\hat{\bar{Y}}_b}
#'  where \eqn{\hat{\bar{Y}}_a = \sum_{k \in s_a}\hat{p}_{ak}y_k, \hat{\bar{Y}}_{ab} = \sum_{k \in s_{ab}^A}\hat{p}_{abk}^Ay_k, \hat{\bar{Y}}_{ab}^B = \sum_{k \in s_{ab}^B}\hat{p}_{abk}^By_k}
#'  and \eqn{\hat{\bar{Y}}_b = \sum_{k \in s_b}\hat{p}_{bk}y_k} with \eqn{\hat{p}_{ak}, \hat{p}_{abk}^A, \hat{p}_{abk}^B} and \eqn{\hat{p}_{bk}} the weights resulting of applying the pseudo empirical likelihood procedure to a determined function under a determined set of constraints, depending on the case. 
#'  Furthermore, \eqn{\eta \in (0,1)}. In this case, \eqn{N_A, N_B} and \eqn{N_{ab}} have been supposed known and no additional auxiliary variables have been considered. This is not happening in some cases.
#'  Function covers following scenarios:
#'  \itemize{
#'	\item There is not any additional auxiliary variable 
#'      \itemize{
#'      	\item \eqn{N_A, N_B} and \eqn{N_{ab}} unknown
#'              \item \eqn{N_A} and \eqn{N_B} known and \eqn{N_{ab}} unknown
#'              \item \eqn{N_A, N_B} and \eqn{N_{ab}} known
#'      }
#'      \item At least, one additional auxiliary variable is available 
#'      \itemize{
#'              \item \eqn{N_A} and \eqn{N_B} known and \eqn{N_{ab}} unknown
#'              \item \eqn{N_A, N_B} and \eqn{N_{ab}} known
#'      }
#'  }
#'  Explicit variance of this estimator is not easy to obtain. Instead, confidence intervals can be computed through the bi-section method. This method constructs intervals in the form \eqn{\{\theta|r_{ns}(\theta) < \chi_1^2(\alpha)\}}, 
#'  where \eqn{\chi_1^2(\alpha)} is the \eqn{1 - \alpha} quantile from a \eqn{\chi^2} distribution with one degree of freedom and \eqn{r_{ns}(\theta)} represents the so called pseudo empirical log likelihood ratio statistic, 
#'  which can be obtained as a difference of two pseudo empirical likelihood functions.
#' @references Rao, J. N. K. and Wu, C. (2010)
#' \emph{Pseudo Empirical Likelihood Inference for Multiple Frame Surveys}.
#' Journal of the American Statistical Association, 105, 1494 - 1503.
#' @references Wu, C. (2005)
#' \emph{Algorithms and R codes for the pseudo empirical likelihood methods in survey sampling}.
#' Survey Methodology, Vol. 31, 2, pp. 239 - 243.
#' @seealso \code{\link{JackPEL}}
#' @examples
#' data(DatA)
#' data(DatB)
#' data(PiklA)
#' data(PiklB)
#' 
#' #Let calculate pseudo empirical likelihood estimator for variable Feeding, without
#' #considering any auxiliary information
#' PEL(DatA$Feed, DatB$Feed, PiklA, PiklB, DatA$Domain, DatB$Domain)
#' 
#' #Now, let calculate pseudo empirical estimator for variable Clothing when the frame
#' #sizes and the overlap domain size are known
#' PEL(DatA$Clo, DatB$Clo, PiklA, PiklB, DatA$Domain, DatB$Domain, 
#' N_A = 1735, N_B = 1191, N_ab = 601)
#' 
#' #Finally, let calculate pseudo empirical likelihood estimator and a 90% confidence interval
#' #for population total for variable Feeding, considering Income and Metres2 as auxiliary 
#' #variables and with frame sizes and overlap domain size known.
#' PEL(DatA$Feed, DatB$Feed, PiklA, PiklB, DatA$Domain, DatB$Domain, 
#' N_A = 1735, N_B =  1191, N_ab = 601, xsAFrameA = DatA$Inc, xsBFrameA = DatB$Inc, 
#' xsAFrameB = DatA$M2, xsBFrameB = DatB$M2, XA = 4300260, XB = 176553, 
#' conf_level = 0.90)
#' @export
PEL = function (ysA, ysB, pi_A, pi_B, domains_A, domains_B, N_A = NULL, N_B = NULL, N_ab = NULL, xsAFrameA = NULL, xsBFrameA = NULL, xsAFrameB = NULL, xsBFrameB = NULL, XA = NULL, XB = NULL, conf_level = NULL)
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
		stop("There are missing values in inclusion probabilities from frame A.")
	if (any(is.na(pi_B)))
		stop("There are missing values in inclusion probabilities from frame B.")
	if (any(is.na(domains_A)))
		stop("There are missing values in domains from frame A.")
	if (any(is.na(domains_B)))
		stop("There are missing values in domains from frame B.")
	if (nrow(ysA) != nrow(pi_A) | nrow(ysA) != length(domains_A) | length(domains_A) != nrow(pi_A))
		stop("Arguments from frame A have different sizes.")
	if (nrow(ysB) != nrow(pi_B) | nrow(ysB) != length(domains_B) | length(domains_B) != nrow(pi_B))
		stop("Arguments from frame B have different sizes.")
	if (length(which(domains_A == "a")) + length(which(domains_A == "ab")) != length(domains_A))
		stop("Domains from frame A are not correct.")
	if (length(which(domains_B == "b")) + length(which(domains_B == "ba")) != length(domains_B))
		stop("Domains from frame B are not correct.")
	if ((is.null (N_A) & !is.null (N_B)) | (!is.null (N_A) & is.null (N_B)))
		stop("Only one value has been indicated for N_A and N_B. This is not valid. Both or none should be indicated.")
	if (!is.null (N_ab) & (is.null (N_A) | is.null (N_B)))
		stop("A value for N_ab has been provided, but any value for N_A or N_B is missing. This is not a possible option.")

	sample <- rbind(ysA, ysB)
	domains <- factor(c(as.character(domains_A), as.character(domains_B)))
	
	cl <- match.call()
	
	n_A <- nrow (ysA)
	n_B <- nrow (ysB)
	n <-  n_A + n_B
	c <- ncol(ysA)
	
	delta_a <- Domains (rep (1, n), domains, "a")
	delta_ab <- Domains (rep (1, n), domains, "ab")
	delta_ba <- Domains (rep (1, n), domains, "ba")
	delta_b <- Domains (rep (1, n), domains, "b")

	ones_a_A <- Domains (rep (1, n_A), domains_A, "a")
	ones_ab_A <- Domains (rep (1, n_A), domains_A, "ab")
	ones_b_B <- Domains (rep (1, n_B), domains_B, "b")
	ones_ab_B <- Domains (rep (1, n_B), domains_B, "ba")

	est <- matrix(, 2, c, dimnames = list(c("Total", "Mean"), cnames))
	totdom <- matrix(, 4, c, dimnames = list(c("Total dom. a", "Total dom. ab", "Total dom. b", "Total dom. ba"), cnames))
	meandom <- matrix(, 4, c, dimnames = list(c("Mean dom. a", "Mean dom. ab", "Mean dom. b", "Mean dom. ba"), cnames))
	par <- 	matrix(, 1, c, dimnames = list("eta", cnames))
	if (is.null(conf_level))
		interv <- NULL
	else
		interv <- matrix(, 6, c, dimnames = list(c("Total", "Lower Bound", "Upper Bound", "Mean", "Lower Bound", "Upper Bound"), cnames))


	if (!is.null(dim(drop(pi_A))) & !is.null(dim(drop(pi_B)))) {

		pik <- c(diag(pi_A), diag(pi_B))	
		di <- 1 / pik

		Nhat_a_A <- HT (ones_a_A, diag(pi_A))
		Nhat_ab_A <- HT (ones_ab_A, diag(pi_A))
		Nhat_b_B <- HT (ones_b_B, diag(pi_B))
		Nhat_ab_B <- HT (ones_ab_B, diag(pi_B))

		di_star <- di / tapply(di, domains, sum)[domains]

		for (k in 1:c) {

			data_ab_A <- Domains (ysA[,k], domains_A, "ab")
			data_ab_B <- Domains (ysB[,k], domains_B, "ba")
			Vhat_Yhat_ab_A <- VarHT (data_ab_A, pi_A)
			Vhat_Yhat_ab_B <- VarHT (data_ab_B, pi_B)

			eta_0 <- Vhat_Yhat_ab_B / (Vhat_Yhat_ab_A + Vhat_Yhat_ab_B)
			par[,k] <- eta_0

			z1 <- delta_a
			z2 <- delta_ab
			z3 <- delta_ba
			z4 <- sample[,k] / eta_0 * delta_ab - sample[,k] / (1 - eta_0) * delta_ba

			if (is.null(xsAFrameA) & is.null(xsBFrameB)) {
	
				z <- as.matrix (cbind(z1, z2, z3, z4))

				if (is.null(N_ab)) {

					Vhat_Nhat_ab_B <- VarHT (ones_ab_B, pi_B)
					Vhat_Nhat_ab_A <- VarHT (ones_ab_A, pi_A)
					phihat <- Vhat_Nhat_ab_B / (Vhat_Nhat_ab_A + Vhat_Nhat_ab_B)
					Nhat_abP <- phihat * Nhat_ab_A + (1 - phihat) * Nhat_ab_B

					if (is.null(N_A) & is.null(N_B)) {

						Nhat_A <- HT (ones_a_A + ones_ab_A, diag(pi_A))
						Nhat_B <- HT (ones_b_B + ones_ab_B, diag(pi_B))
						Nhat_P <- Nhat_A + Nhat_B - Nhat_abP

						Wa <- (Nhat_A - Nhat_abP) / Nhat_P
						Wab_eta_0 <- eta_0 * Nhat_abP / Nhat_P
						Wba_eta_0 <- (1 - eta_0) * Nhat_abP / Nhat_P
						Wb <- (Nhat_B - Nhat_abP) / Nhat_P
						domain_size_estimation <- c(Nhat_A - Nhat_abP, Nhat_abP, Nhat_B - Nhat_abP, Nhat_abP)

						zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0))
					}
					else {
	
						Nhat_P <- N_A + N_B - Nhat_abP

						Wa <- (N_A - Nhat_abP) / Nhat_P
						Wab_eta_0 <- eta_0 * Nhat_abP / Nhat_P
						Wba_eta_0 <- (1 - eta_0) * Nhat_abP / Nhat_P
						Wb <- (N_B - Nhat_abP) / Nhat_P
						domain_size_estimation <- c(N_A - Nhat_abP, Nhat_abP, N_B - Nhat_abP, Nhat_abP)

						zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0))
					}

					Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)
					lambda <- Lag2 (z, di, zmean)

					phat <- as.vector(di_star / as.vector((rep(1,n) + crossprod(t(z - Zmean), lambda))))
					phat_star <- as.vector(phat / tapply(phat, domains, sum)[domains])

					Yhatmean_a <- sum (phat_star * sample[,k] * delta_a)
					Yhatmean_ab <- sum (phat_star * sample[,k] * delta_ab)
					Yhatmean_ba <- sum (phat_star * sample[,k] * delta_ba)
					Yhatmean_b <- sum (phat_star * sample[,k] * delta_b)

					meandom[,k] <- c(Yhatmean_a, Yhatmean_ab, Yhatmean_b, Yhatmean_ba)
					totdom[,k] <- meandom[,k] * domain_size_estimation

					mean_estimation <- Wa * Yhatmean_a + Wab_eta_0 * Yhatmean_ab + Wba_eta_0 * Yhatmean_ba + Wb * Yhatmean_b
					total_estimation <- mean_estimation * Nhat_P
					est[,k] <- c(total_estimation, mean_estimation)
	
					if (!is.null(conf_level)) {
					
						Di_star <- diag(di_star)

						suma1_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), z[domains == "a",] - Zmean[domains == "a",])
						suma2_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), sample[domains == "a",k] - mean_estimation)
						suma1_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), z[domains == "ab",] - Zmean[domains == "ab",])
						suma2_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), sample[domains == "ab",k] - mean_estimation)
						suma1_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), z[domains == "b",] - Zmean[domains == "b",])
						suma2_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), sample[domains == "b",k] - mean_estimation)
						suma1_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), z[domains == "ba",] - Zmean[domains == "ba",])
						suma2_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), sample[domains == "ba",k] - mean_estimation)
						Bhat <- solve(Wa * suma1_a + Wab_eta_0 * suma1_ab + Wba_eta_0 * suma1_ba + Wb * suma1_b, Wa * suma2_a + Wab_eta_0 * suma2_ab + Wba_eta_0 * suma2_ba + Wb * suma2_b)

						rhat <- sample[,k] - crossprod(t(z), Bhat)
						Rhatmean_a <- sum (phat_star[domains == "a"] * rhat[domains == "a"])
						Rhatmean_ab <- sum (phat_star[domains == "ab"] * rhat[domains == "ab"])
						Rhatmean_b <- sum (phat_star[domains == "b"] * rhat[domains == "b"])
						Rhatmean_ba <- sum (phat_star[domains == "ba"] * rhat[domains == "ba"])
						r_tilde <- (rhat - Rhatmean_a) * delta_a + (eta_0 * (rhat - Rhatmean_ab)) * delta_ab + (rhat - Rhatmean_b) * delta_b + ((1 - eta_0) * (rhat - Rhatmean_ba)) * delta_ba
						Vhat_A_Ar <- 1 / Nhat_P^2 * VarHT (r_tilde[domains == "a" | domains == "ab"], pi_A)
						Vhat_B_Br <- 1 / Nhat_P^2 * VarHT (r_tilde[domains == "b" | domains == "ba"], pi_B)

						rhatmean <- mean_estimation - crossprod(Bhat, zmean)
						Rhatmean <- rep (rhatmean, n)
						suma3_a <- crossprod(di_star[domains == "a"], (rhat[domains == "a"] - Rhatmean[domains == "a"])^2)
						suma3_ab <- crossprod(di_star[domains == "ab"], (rhat[domains == "ab"] - Rhatmean[domains == "ab"])^2)
						suma3_b <- crossprod(di_star[domains == "b"], (rhat[domains == "b"] - Rhatmean[domains == "b"])^2)
						suma3_ba <- crossprod(di_star[domains == "ba"], (rhat[domains == "ba"] - Rhatmean[domains == "ba"])^2)
						deff_P <- (Vhat_A_Ar + Vhat_B_Br) / (1 / n * (Wa * suma3_a + Wab_eta_0 * suma3_ab + Wba_eta_0 * suma3_ba + Wb * suma3_b))

						conf_interv <- PELConfInt (conf_level, sample[,k], di_star, mean_estimation, n/deff_P, di_star)
						upper <- conf_interv[1]
						lower <- conf_interv[2]
						interv[,k] <- c(total_estimation, lower * Nhat_P, upper * Nhat_P, mean_estimation, lower, upper)	
					}
				}
				else {
	
					N <- N_A + N_B - N_ab

					Wa <- (N_A - N_ab) / N
					Wab_eta_0 <- eta_0 * N_ab / N
					Wba_eta_0 <- (1 - eta_0) * N_ab / N
					Wb <- (N_B - N_ab) / N
					domain_size_estimation <- c(N_A - N_ab, N_ab, N_B - N_ab, N_ab)

					zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0))
					Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)
					lambda <- Lag2 (z, di, zmean)
		
					phat <- as.vector(di_star / as.vector((rep(1,n) + crossprod(t(z - Zmean), lambda))))
					phat_star <- as.vector(phat / tapply(phat, domains, sum)[domains])

					Yhatmean_a <- sum (phat_star * Domains(sample[,k], domains, "a"))
					Yhatmean_ab <- sum (phat_star * Domains(sample[,k], domains, "ab"))
					Yhatmean_b <- sum (phat_star * Domains(sample[,k], domains, "b"))
					Yhatmean_ba <- sum (phat_star * Domains(sample[,k], domains, "ba"))

					meandom[,k] <- c(Yhatmean_a, Yhatmean_ab, Yhatmean_b, Yhatmean_ba)
					totdom[,k] <- meandom[,k] * domain_size_estimation

					mean_estimation <- Wa * Yhatmean_a + Wab_eta_0 * Yhatmean_ab +  Wba_eta_0 * Yhatmean_ba + Wb * Yhatmean_b
					total_estimation <- N * mean_estimation
					est[,k] <- c(total_estimation, mean_estimation)

					if (!is.null(conf_level)) {

						Di_star <- diag(di_star)

						suma1_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), z[domains == "a",] - Zmean[domains == "a",])
						suma2_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), sample[domains == "a",k] - mean_estimation)
						suma1_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), z[domains == "ab",] - Zmean[domains == "ab",])
						suma2_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), sample[domains == "ab",k] - mean_estimation)
						suma1_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), z[domains == "b",] - Zmean[domains == "b",])
						suma2_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), sample[domains == "b",k] - mean_estimation)
						suma1_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), z[domains == "ba",] - Zmean[domains == "ba",])
						suma2_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), sample[domains == "ba",k] - mean_estimation)
						Bhat <- solve(Wa * suma1_a + Wab_eta_0 * suma1_ab + Wba_eta_0 * suma1_ba + Wb * suma1_b, Wa * suma2_a + Wab_eta_0 * suma2_ab + Wba_eta_0 * suma2_ba + Wb * suma2_b)

						rhat <- sample[,k] - crossprod(t(z), Bhat)
						Rhatmean_a <- sum (phat_star[domains == "a"] * rhat[domains == "a"])
						Rhatmean_ab <- sum (phat_star[domains == "ab"] * rhat[domains == "ab"])
						Rhatmean_b <- sum (phat_star[domains == "b"] * rhat[domains == "b"])
						Rhatmean_ba <- sum (phat_star[domains == "ba"] * rhat[domains == "ba"])
						r_tilde <- (rhat - Rhatmean_a) * delta_a + (eta_0 * (rhat - Rhatmean_ab)) * delta_ab + (rhat - Rhatmean_b) * delta_b + ((1 - eta_0) * (rhat - Rhatmean_ba)) * delta_ba
						Vhat_A_Ar <- 1 / N^2 * VarHT (r_tilde[domains == "a" | domains == "ab"], pi_A)
						Vhat_B_Br <- 1 / N^2 * VarHT (r_tilde[domains == "b" | domains == "ba"], pi_B)

						rhatmean <- mean_estimation - crossprod(Bhat, zmean)
						Rhatmean <- rep (rhatmean, n)
						suma3_a <- crossprod(di_star[domains == "a"], (rhat[domains == "a"] - Rhatmean[domains == "a"])^2)
						suma3_ab <- crossprod(di_star[domains == "ab"], (rhat[domains == "ab"] - Rhatmean[domains == "ab"])^2)
						suma3_b <- crossprod(di_star[domains == "b"], (rhat[domains == "b"] - Rhatmean[domains == "b"])^2)
						suma3_ba <- crossprod(di_star[domains == "ba"], (rhat[domains == "ba"] - Rhatmean[domains == "ba"])^2)
						deff_P <- (Vhat_A_Ar + Vhat_B_Br) / (1 / n * (Wa * suma3_a + Wab_eta_0 * suma3_ab + Wba_eta_0 * suma3_ba + Wb * suma3_b))

						conf_interv <- PELConfInt (conf_level, sample[,k], di_star, mean_estimation, n/deff_P, di_star)
						upper <- conf_interv[1]
						lower <- conf_interv[2]
						interv[,k] <- c(total_estimation, N * lower, N * upper, mean_estimation, lower, upper)	
					}
				}
			}
			else {

				if (is.null(N_ab)) {

					Vhat_Nhat_ab_B <- VarHT (ones_ab_B, pi_B)
					Vhat_Nhat_ab_A <- VarHT (ones_ab_A, pi_A)
					phihat <- Vhat_Nhat_ab_B / (Vhat_Nhat_ab_A + Vhat_Nhat_ab_B)
					Nhat_abP <- phihat * Nhat_ab_A + (1 - phihat) * Nhat_ab_B
					Nhat_P <- N_A + N_B - Nhat_abP

					Wa <- (N_A - Nhat_abP) / Nhat_P
					Wab_eta_0 <- eta_0 * Nhat_abP / Nhat_P
					Wba_eta_0 <- (1 - eta_0) * Nhat_abP / Nhat_P
					Wb <- (N_B - Nhat_abP) / Nhat_P
					domain_size_estimation <- c(N_A - Nhat_abP, Nhat_abP, N_B - Nhat_abP, Nhat_abP)

					if (is.null(xsAFrameA)){
					
						xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
						XFrameB <- rbind(xsAFrameB, xsBFrameB)
						z5 <- matrix(0, n, ncol(XFrameB))
						z5[domains == "b",] <- XFrameB[domains == "b",]
						z5[domains == "ba",] <- XFrameB[domains == "ba",]/(1 - eta_0)

						z <- as.matrix (cbind(z1, z2, z3, z4, z5))
						zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0, XB / Nhat_P))
						Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)
					}
					else {
				
						if (is.null(xsBFrameB)){
						
							xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
							XFrameA <- rbind(xsAFrameA, xsBFrameA)
							z5 <- matrix(0, n, ncol(XFrameA))
							z5[domains == "a",] <- XFrameA[domains == "a",]
							z5[domains == "ab",] <- XFrameA[domains == "ab",]/eta_0

							z <- as.matrix (cbind(z1, z2, z3, z4, z5))
							zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0, XA / Nhat_P))
							Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)	
						}
						else{
						
							xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
							XFrameA <- rbind(xsAFrameA, xsBFrameA)
							xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
							XFrameB <- rbind(xsAFrameB, xsBFrameB)
							z5 <- matrix(0, n, ncol(XFrameA))
							z5[domains == "a",] <- XFrameA[domains == "a",]
							z5[domains == "ab",] <- XFrameA[domains == "ab",]/eta_0
							z6 <- matrix(0, n, ncol(XFrameB))
							z6[domains == "b",] <- XFrameB[domains == "b",]
							z6[domains == "ba",] <- XFrameB[domains == "ba",]/(1 - eta_0)

							z <- as.matrix (cbind(z1, z2, z3, z4, z5, z6))
							zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0, XA / Nhat_P, XB / Nhat_P))
							Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)
						}
					}

					lambda <- Lag2 (z, di, zmean)
					phat <- as.vector(di_star / as.vector((rep(1,n) + crossprod(t(z - Zmean), lambda))))
					phat_star <- as.vector(phat / tapply(phat, domains, sum)[domains])

					Yhatmean_a <- sum (phat_star * Domains(sample[,k], domains, "a"))
					Yhatmean_ab <- sum (phat_star * Domains(sample[,k], domains, "ab"))
					Yhatmean_b <- sum (phat_star * Domains(sample[,k], domains, "b"))
					Yhatmean_ba <- sum (phat_star * Domains(sample[,k], domains, "ba"))

					meandom[,k] <- c(Yhatmean_a, Yhatmean_ab, Yhatmean_b, Yhatmean_ba)
					totdom[,k] <- meandom[,k] * domain_size_estimation

					mean_estimation <- Wa * Yhatmean_a + Wab_eta_0 * Yhatmean_ab + Wba_eta_0 * Yhatmean_ba + Wb * Yhatmean_b
					total_estimation <- Nhat_P * mean_estimation
					est[,k] <- c(total_estimation, mean_estimation)

					if (!is.null(conf_level)) {

						Di_star <- diag(di_star)

						suma1_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), z[domains == "a",] - Zmean[domains == "a",])
						suma2_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), sample[domains == "a",k] - mean_estimation)
						suma1_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), z[domains == "ab",] - Zmean[domains == "ab",])
						suma2_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), sample[domains == "ab",k] - mean_estimation)
						suma1_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), z[domains == "b",] - Zmean[domains == "b",])
						suma2_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), sample[domains == "b",k] - mean_estimation)
						suma1_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), z[domains == "ba",] - Zmean[domains == "ba",])
						suma2_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), sample[domains == "ba",k] - mean_estimation)
						Bhat <- solve(Wa * suma1_a + Wab_eta_0 * suma1_ab + Wba_eta_0 * suma1_ba + Wb * suma1_b, Wa * suma2_a + Wab_eta_0 * suma2_ab + Wba_eta_0 * suma2_ba + Wb * suma2_b)

						rhat <- sample[,k] - crossprod(t(z), Bhat)
						Rhatmean_a <- sum (phat_star[domains == "a"] * rhat[domains == "a"])
						Rhatmean_ab <- sum (phat_star[domains == "ab"] * rhat[domains == "ab"])
						Rhatmean_b <- sum (phat_star[domains == "b"] * rhat[domains == "b"])
						Rhatmean_ba <- sum (phat_star[domains == "ba"] * rhat[domains == "ba"])
						r_tilde <- (rhat - Rhatmean_a) * delta_a + (eta_0 * (rhat - Rhatmean_ab)) * delta_ab + (rhat - Rhatmean_b) * delta_b + ((1 - eta_0) * (rhat - Rhatmean_ba)) * delta_ba
						Vhat_A_Ar <- 1 / Nhat_P^2 * VarHT (r_tilde[domains == "a" | domains == "ab"], pi_A)
						Vhat_B_Br <- 1 / Nhat_P^2 * VarHT (r_tilde[domains == "b" | domains == "ba"], pi_B)

						rhatmean <- mean_estimation - crossprod(Bhat, zmean)
						Rhatmean <- rep (rhatmean, n)
						suma3_a <- crossprod(di_star[domains == "a"], (rhat[domains == "a"] - Rhatmean[domains == "a"])^2)
						suma3_ab <- crossprod(di_star[domains == "ab"], (rhat[domains == "ab"] - Rhatmean[domains == "ab"])^2)
						suma3_b <- crossprod(di_star[domains == "b"], (rhat[domains == "b"] - Rhatmean[domains == "b"])^2)
						suma3_ba <- crossprod(di_star[domains == "ba"], (rhat[domains == "ba"] - Rhatmean[domains == "ba"])^2)
						deff_P <- (Vhat_A_Ar + Vhat_B_Br) / (1 / n * (Wa * suma3_a + Wab_eta_0 * suma3_ab + Wba_eta_0 * suma3_ba + Wb * suma3_b))

						conf_interv <- PELConfInt (conf_level, sample[,k], di_star, mean_estimation, n/deff_P, di_star)
						upper <- conf_interv[1]
						lower <- conf_interv[2]
						interv[,k] <- c(total_estimation, lower * Nhat_P, upper * Nhat_P, mean_estimation, lower, upper)	
					}
				}
				else {
			
					N <- N_A + N_B - N_ab

					Wa <- (N_A - N_ab) / N
					Wab_eta_0 <- eta_0 * N_ab / N
					Wba_eta_0 <- (1 - eta_0) * N_ab / N
					Wb <- (N_B - N_ab) / N
					domain_size_estimation <- c(N_A - N_ab, N_ab, N_B - N_ab, N_ab)

					if (is.null(xsAFrameA)){
					
						xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
						XFrameB <- rbind(xsAFrameB, xsBFrameB)
						z5 <- matrix(0, n, ncol(XFrameB))
						z5[domains == "b",] <- XFrameB[domains == "b",]
						z5[domains == "ba",] <- XFrameB[domains == "ba",]/(1 - eta_0)

						z <- as.matrix (cbind(z1, z2, z3, z4, z5))
						zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0, XB / N))
						Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)
					}
					else {
				
						if (is.null(xsBFrameB)){
					
							xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
							XFrameA <- rbind(xsAFrameA, xsBFrameA)
							z5 <- matrix(0, n, ncol(XFrameA))
							z5[domains == "a",] <- XFrameA[domains == "a",]
							z5[domains == "ab",] <- XFrameA[domains == "ab",]/eta_0

							z <- as.matrix (cbind(z1, z2, z3, z4, z5))
							zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0, XA / N))
							Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)	
						}
						else{
						
							xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
							XFrameA <- rbind(xsAFrameA, xsBFrameA)
							xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
							XFrameB <- rbind(xsAFrameB, xsBFrameB)
							z5 <- matrix(0, n, ncol(XFrameA))
							z5[domains == "a",] <- XFrameA[domains == "a",]
							z5[domains == "ab",] <- XFrameA[domains == "ab",]/eta_0
							z6 <- matrix(0, n, ncol(XFrameB))
							z6[domains == "b",] <- XFrameB[domains == "b",]
							z6[domains == "ba",] <- XFrameB[domains == "ba",]/(1 - eta_0)

							z <- as.matrix (cbind(z1, z2, z3, z4, z5, z6))
							zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0, XA / N, XB / N))
							Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)
						}

					}

					lambda <- Lag2 (z, di, zmean)
					phat <- as.vector(di_star / as.vector((rep(1,n) + crossprod(t(z - Zmean), lambda))))
					phat_star <- as.vector(phat / tapply(phat, domains, sum)[domains])

					Yhatmean_a <- sum (phat_star * Domains(sample[,k], domains, "a"))
					Yhatmean_ab <- sum (phat_star * Domains(sample[,k], domains, "ab"))
					Yhatmean_b <- sum (phat_star * Domains(sample[,k], domains, "b"))
					Yhatmean_ba <- sum (phat_star * Domains(sample[,k], domains, "ba"))

					meandom[,k] <- c(Yhatmean_a, Yhatmean_ab, Yhatmean_b, Yhatmean_ba)
					totdom[,k] <- meandom[,k] * domain_size_estimation

					mean_estimation <- Wa * Yhatmean_a + Wab_eta_0 * Yhatmean_ab + Wba_eta_0 * Yhatmean_ba + Wb * Yhatmean_b
					total_estimation <- N * mean_estimation
					est[,k] <- c(total_estimation, mean_estimation)

					if (!is.null(conf_level)) {

						Di_star <- diag(di_star)

						suma1_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), z[domains == "a",] - Zmean[domains == "a",])
						suma2_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), sample[domains == "a",k] - mean_estimation)
						suma1_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), z[domains == "ab",] - Zmean[domains == "ab",])
						suma2_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), sample[domains == "ab",k] - mean_estimation)
						suma1_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), z[domains == "b",] - Zmean[domains == "b",])
						suma2_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), sample[domains == "b",k] - mean_estimation)
						suma1_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), z[domains == "ba",] - Zmean[domains == "ba",])
						suma2_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), sample[domains == "ba",k] - mean_estimation)
						Bhat <- solve(Wa * suma1_a + Wab_eta_0 * suma1_ab + Wba_eta_0 * suma1_ba + Wb * suma1_b, Wa * suma2_a + Wab_eta_0 * suma2_ab + Wba_eta_0 * suma2_ba + Wb * suma2_b)

						rhat <- sample[,k] - z %*% Bhat
						Rhatmean_a <- sum (phat_star[domains == "a"] * rhat[domains == "a"])
						Rhatmean_ab <- sum (phat_star[domains == "ab"] * rhat[domains == "ab"])
						Rhatmean_b <- sum (phat_star[domains == "b"] * rhat[domains == "b"])
						Rhatmean_ba <- sum (phat_star[domains == "ba"] * rhat[domains == "ba"])
						r_tilde <- (rhat - Rhatmean_a) * delta_a + (eta_0 * (rhat - Rhatmean_ab)) * delta_ab + (rhat - Rhatmean_b) * delta_b + ((1 - eta_0) * (rhat - Rhatmean_ba)) * delta_ba
						Vhat_A_Ar <- 1 / N^2 * VarHT (r_tilde[domains == "a" | domains == "ab"], pi_A)
						Vhat_B_Br <- 1 / N^2 * VarHT (r_tilde[domains == "b" | domains == "ba"], pi_B)

						rhatmean <- mean_estimation - crossprod(Bhat, zmean)
						Rhatmean <- rep (rhatmean, n)
						suma3_a <- crossprod(di_star[domains == "a"], (rhat[domains == "a"] - Rhatmean[domains == "a"])^2)
						suma3_ab <- crossprod(di_star[domains == "ab"], (rhat[domains == "ab"] - Rhatmean[domains == "ab"])^2)
						suma3_b <- crossprod(di_star[domains == "b"], (rhat[domains == "b"] - Rhatmean[domains == "b"])^2)
						suma3_ba <- crossprod(di_star[domains == "ba"], (rhat[domains == "ba"] - Rhatmean[domains == "ba"])^2)
						deff_P <- (Vhat_A_Ar + Vhat_B_Br) / (1 / n * (Wa * suma3_a + Wab_eta_0 * suma3_ab + Wba_eta_0 * suma3_ba + Wb * suma3_b))

						conf_interv <- PELConfInt (conf_level, sample[,k], di_star, mean_estimation, n/deff_P, di_star)
						upper <- conf_interv[1]
						lower <- conf_interv[2]
						interv[,k] <- c(total_estimation, lower * N, upper * N, mean_estimation, lower, upper)	
					}		
				}
			}
		}
	}
	else {

		if(is.null(dim(drop(pi_A))) & is.null(dim(drop(pi_B)))) {
		
			pik <- c(pi_A, pi_B)	
			di <- 1 / pik

			Nhat_a_A <- HT (ones_a_A, pi_A)
			Nhat_ab_A <- HT (ones_ab_A, pi_A)
			Nhat_b_B <- HT (ones_b_B, pi_B)
			Nhat_ab_B <- HT (ones_ab_B, pi_B)

			di_star <- rep (0, n)
			for (i in 1:n) {
				if (domains[i] == "a")
					di_star[i] <- di[i] / sum(di[domains == "a"])
				if (domains[i] == "ab")
					di_star[i] <- di[i] / sum (di[domains == "ab"])
				if (domains[i] == "b")
					di_star[i] <- di[i] / sum(di[domains == "b"])
				if (domains[i] == "ba")
					di_star[i] <- di[i] / sum(di[domains == "ba"])
			}

			for (k in 1:c) {

				data_ab_A <- Domains (ysA[,k], domains_A, "ab")
				data_ab_B <- Domains (ysB[,k], domains_B, "ba")
				Vhat_Yhat_ab_A <- varest (data_ab_A, pik = pi_A)
				Vhat_Yhat_ab_B <- varest (data_ab_B, pik= pi_B)

				eta_0 <- Vhat_Yhat_ab_B / (Vhat_Yhat_ab_A + Vhat_Yhat_ab_B)
				par[,k] <- eta_0

				z1 <- delta_a
				z2 <- delta_ab
				z3 <- delta_ba
				z4 <- sample[,k] / eta_0 * delta_ab - sample[,k] / (1 - eta_0) * delta_ba

				if (is.null(xsAFrameA) & is.null(xsBFrameB)) {
	
					z <- as.matrix (cbind(z1, z2, z3, z4))

					if (is.null(N_ab)) {

						Vhat_Nhat_ab_B <- varest (ones_ab_B, pik = pi_B)
						Vhat_Nhat_ab_A <- varest (ones_ab_A, pik = pi_A)
						phihat <- Vhat_Nhat_ab_B / (Vhat_Nhat_ab_A + Vhat_Nhat_ab_B)
						Nhat_abP <- phihat * Nhat_ab_A + (1 - phihat) * Nhat_ab_B

						if (is.null(N_A) & is.null(N_B)) {

							Nhat_A <- HT (ones_a_A + ones_ab_A, pi_A)
							Nhat_B <- HT (ones_b_B + ones_ab_B, pi_B)
							Nhat_P <- Nhat_A + Nhat_B - Nhat_abP
	
							Wa <- (Nhat_A - Nhat_abP) / Nhat_P
							Wab_eta_0 <- eta_0 * Nhat_abP / Nhat_P
							Wba_eta_0 <- (1 - eta_0) * Nhat_abP / Nhat_P
							Wb <- (Nhat_B - Nhat_abP) / Nhat_P
							domain_size_estimation <- c(Nhat_A - Nhat_abP, Nhat_abP, Nhat_B - Nhat_abP, Nhat_abP)
	
							zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0))
						}
						else {
	
							Nhat_P <- N_A + N_B - Nhat_abP

							Wa <- (N_A - Nhat_abP) / Nhat_P
							Wab_eta_0 <- eta_0 * Nhat_abP / Nhat_P
							Wba_eta_0 <- (1 - eta_0) * Nhat_abP / Nhat_P
							Wb <- (N_B - Nhat_abP) / Nhat_P
							domain_size_estimation <- c(N_A - Nhat_abP, Nhat_abP, N_B - Nhat_abP, Nhat_abP)

							zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0))
						}

						Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)
						lambda <- Lag2 (z, di, zmean)

						phat <- as.vector(di_star / as.vector((rep(1,n) + crossprod(t(z - Zmean), lambda))))
						phat_star <- as.vector(phat / tapply(phat, domains, sum)[domains])

						Yhatmean_a <- sum (phat_star * sample[,k] * delta_a)
						Yhatmean_ab <- sum (phat_star * sample[,k] * delta_ab)
						Yhatmean_ba <- sum (phat_star * sample[,k] * delta_ba)
						Yhatmean_b <- sum (phat_star * sample[,k] * delta_b)

						meandom[,k] <- c(Yhatmean_a, Yhatmean_ab, Yhatmean_b, Yhatmean_ba)
						totdom[,k] <- meandom[,k] * domain_size_estimation

						mean_estimation <-  Wa * Yhatmean_a + Wab_eta_0 * Yhatmean_ab + Wba_eta_0 * Yhatmean_ba + Wb * Yhatmean_b
						total_estimation <- Nhat_P * mean_estimation
						est[,k] <- c(total_estimation, mean_estimation)
	
						if (!is.null(conf_level)) {
					
							Di_star <- diag(di_star)

							suma1_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), z[domains == "a",] - Zmean[domains == "a",])
							suma2_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), sample[domains == "a",k] - mean_estimation)
							suma1_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), z[domains == "ab",] - Zmean[domains == "ab",])
							suma2_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), sample[domains == "ab",k] - mean_estimation)
							suma1_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), z[domains == "b",] - Zmean[domains == "b",])
							suma2_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), sample[domains == "b",k] - mean_estimation)
							suma1_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), z[domains == "ba",] - Zmean[domains == "ba",])
							suma2_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), sample[domains == "ba",k] - mean_estimation)
							Bhat <- solve(Wa * suma1_a + Wab_eta_0 * suma1_ab + Wba_eta_0 * suma1_ba + Wb * suma1_b, Wa * suma2_a + Wab_eta_0 * suma2_ab + Wba_eta_0 * suma2_ba + Wb * suma2_b)

							rhat <- sample[,k] - crossprod(t(z), Bhat)
							Rhatmean_a <- sum (phat_star[domains == "a"] * rhat[domains == "a"])
							Rhatmean_ab <- sum (phat_star[domains == "ab"] * rhat[domains == "ab"])
							Rhatmean_b <- sum (phat_star[domains == "b"] * rhat[domains == "b"])
							Rhatmean_ba <- sum (phat_star[domains == "ba"] * rhat[domains == "ba"])
							r_tilde <- (rhat - Rhatmean_a) * delta_a + (eta_0 * (rhat - Rhatmean_ab)) * delta_ab + (rhat - Rhatmean_b) * delta_b + ((1 - eta_0) * (rhat - Rhatmean_ba)) * delta_ba
							Vhat_A_Ar <- 1 / Nhat_P^2 * varest (r_tilde[domains == "a" | domains == "ab"], pik = pi_A)
							Vhat_B_Br <- 1 / Nhat_P^2 * varest (r_tilde[domains == "b" | domains == "ba"], pik = pi_B)

							rhatmean <- mean_estimation - crossprod(Bhat, zmean)
							Rhatmean <- rep (rhatmean, n)
							suma3_a <- crossprod(di_star[domains == "a"], (rhat[domains == "a"] - Rhatmean[domains == "a"])^2)
							suma3_ab <- crossprod(di_star[domains == "ab"], (rhat[domains == "ab"] - Rhatmean[domains == "ab"])^2)
							suma3_b <- crossprod(di_star[domains == "b"], (rhat[domains == "b"] - Rhatmean[domains == "b"])^2)
							suma3_ba <- crossprod(di_star[domains == "ba"], (rhat[domains == "ba"] - Rhatmean[domains == "ba"])^2)
							deff_P <- (Vhat_A_Ar + Vhat_B_Br) / (1 / n * (Wa * suma3_a + Wab_eta_0 * suma3_ab + Wba_eta_0 * suma3_ba + Wb * suma3_b))
	
							conf_interv <- PELConfInt (conf_level, sample[,k], di_star, mean_estimation, n/deff_P, di_star)
							upper <- conf_interv[1]
							lower <- conf_interv[2]
							interv[,k] <- c(total_estimation, lower * Nhat_P, upper * Nhat_P, mean_estimation, lower, upper)	
						}
					}
					else {
	
						N <- N_A + N_B - N_ab

						Wa <- (N_A - N_ab) / N
						Wab_eta_0 <- eta_0 * N_ab / N
						Wba_eta_0 <- (1 - eta_0) * N_ab / N
						Wb <- (N_B - N_ab) / N
						domain_size_estimation <- c(N_A - N_ab, N_ab, N_B - N_ab, N_ab)

						zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0))
						Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)
						lambda <- Lag2 (z, di, zmean)
						phat <- as.vector(di_star / as.vector((rep(1,n) + crossprod(t(z - Zmean), lambda))))
						phat_star <- as.vector(phat / tapply(phat, domains, sum)[domains])

						Yhatmean_a <- sum (phat_star * Domains(sample[,k], domains, "a"))
						Yhatmean_ab <- sum (phat_star * Domains(sample[,k], domains, "ab"))
						Yhatmean_b <- sum (phat_star * Domains(sample[,k], domains, "b"))
						Yhatmean_ba <- sum (phat_star * Domains(sample[,k], domains, "ba"))

						meandom[,k] <- c(Yhatmean_a, Yhatmean_ab, Yhatmean_b, Yhatmean_ba)
						totdom[,k] <- meandom[,k] * domain_size_estimation

						mean_estimation <- Wa * Yhatmean_a + Wab_eta_0 * Yhatmean_ab +  Wba_eta_0 * Yhatmean_ba + Wb * Yhatmean_b
						total_estimation <- N * mean_estimation
						est[,k] <- c(total_estimation, mean_estimation)
	
						if (!is.null(conf_level)) {

							Di_star <- diag(di_star)

							suma1_a <- (t(z[domains == "a",] - Zmean[domains == "a",])) %*% Di_star[domains == "a", domains == "a"] %*% (z[domains == "a",] - Zmean[domains == "a",])
							suma2_a <- t(z[domains == "a",] - Zmean[domains == "a",]) %*% Di_star[domains == "a", domains == "a"] %*% (sample[domains == "a",k] - mean_estimation)
							suma1_ab <- (t(z[domains == "ab",] - Zmean[domains == "ab",])) %*% Di_star[domains == "ab", domains == "ab"] %*% (z[domains == "ab",] - Zmean[domains == "ab",])
							suma2_ab <- t(z[domains == "ab",] - Zmean[domains == "ab",]) %*% Di_star[domains == "ab", domains == "ab"] %*% (sample[domains == "ab",k] - mean_estimation)
							suma1_b <- (t(z[domains == "b",] - Zmean[domains == "b",])) %*% Di_star[domains == "b", domains == "b"] %*% (z[domains == "b",] - Zmean[domains == "b",])
							suma2_b <- t(z[domains == "b",] - Zmean[domains == "b",]) %*% Di_star[domains == "b", domains == "b"] %*% (sample[domains == "b",k] - mean_estimation)
							suma1_ba <- (t(z[domains == "ba",] - Zmean[domains == "ba",])) %*% Di_star[domains == "ba", domains == "ba"] %*% (z[domains == "ba",] - Zmean[domains == "ba",])
							suma2_ba <- t(z[domains == "ba",] - Zmean[domains == "ba",]) %*% Di_star[domains == "ba", domains == "ba"] %*% (sample[domains == "ba",k] - mean_estimation)
							Bhat <- solve(Wa * suma1_a + Wab_eta_0 * suma1_ab + Wba_eta_0 * suma1_ba + Wb * suma1_b) %*% (Wa * suma2_a + Wab_eta_0 * suma2_ab + Wba_eta_0 * suma2_ba + Wb * suma2_b)
							
							rhat <- sample[,k] - crossprod(t(z), Bhat)
							Rhatmean_a <- sum (phat_star[domains == "a"] * rhat[domains == "a"])
							Rhatmean_ab <- sum (phat_star[domains == "ab"] * rhat[domains == "ab"])
							Rhatmean_b <- sum (phat_star[domains == "b"] * rhat[domains == "b"])
							Rhatmean_ba <- sum (phat_star[domains == "ba"] * rhat[domains == "ba"])
							r_tilde <- (rhat - Rhatmean_a) * delta_a + (eta_0 * (rhat - Rhatmean_ab)) * delta_ab + (rhat - Rhatmean_b) * delta_b + ((1 - eta_0) * (rhat - Rhatmean_ba)) * delta_ba
							Vhat_A_Ar <- 1 / N^2 * varest (r_tilde[domains == "a" | domains == "ab"], pik = pi_A)
							Vhat_B_Br <- 1 / N^2 * varest (r_tilde[domains == "b" | domains == "ba"], pik = pi_B)

							rhatmean <- mean_estimation - crossprod(Bhat, zmean)
							Rhatmean <- rep (rhatmean, n)
							suma3_a <- crossprod(di_star[domains == "a"], (rhat[domains == "a"] - Rhatmean[domains == "a"])^2)
							suma3_ab <- crossprod(di_star[domains == "ab"], (rhat[domains == "ab"] - Rhatmean[domains == "ab"])^2)
							suma3_b <- crossprod(di_star[domains == "b"], (rhat[domains == "b"] - Rhatmean[domains == "b"])^2)
							suma3_ba <- crossprod(di_star[domains == "ba"], (rhat[domains == "ba"] - Rhatmean[domains == "ba"])^2)
							deff_P <- (Vhat_A_Ar + Vhat_B_Br) / (1 / n * (Wa * suma3_a + Wab_eta_0 * suma3_ab + Wba_eta_0 * suma3_ba + Wb * suma3_b))

							conf_interv <- PELConfInt (conf_level, sample[,k], di_star, mean_estimation, n/deff_P, di_star)
							upper <- conf_interv[1]
							lower <- conf_interv[2]
							interv[,k] <- c(total_estimation, lower * N, upper * N, mean_estimation, lower, upper)	
						}
					}
				}
				else {

					if (is.null(N_ab)) {
	
						Vhat_Nhat_ab_B <- varest (ones_ab_B, pik = pi_B)
						Vhat_Nhat_ab_A <- varest (ones_ab_A, pik = pi_A)
						phihat <- Vhat_Nhat_ab_B / (Vhat_Nhat_ab_A + Vhat_Nhat_ab_B)
						Nhat_abP <- phihat * Nhat_ab_A + (1 - phihat) * Nhat_ab_B
						Nhat_P <- N_A + N_B - Nhat_abP

						Wa <- (N_A - Nhat_abP) / Nhat_P
						Wab_eta_0 <- eta_0 * Nhat_abP / Nhat_P
						Wba_eta_0 <- (1 - eta_0) * Nhat_abP / Nhat_P
						Wb <- (N_B - Nhat_abP) / Nhat_P
						domain_size_estimation <- c(N_A - Nhat_abP, Nhat_abP, N_B - Nhat_abP, Nhat_abP)

						if (is.null(xsAFrameA)){
					
							xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
							XFrameB <- rbind(xsAFrameB, xsBFrameB)
							z5 <- matrix(0, n, ncol(XFrameB))
							z5[domains == "b",] <- XFrameB[domains == "b",]
							z5[domains == "ba",] <- XFrameB[domains == "ba",]/(1 - eta_0)

							z <- as.matrix (cbind(z1, z2, z3, z4, z5))
							zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0, XB / Nhat_P))
							Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)
						}
						else {
				
							if (is.null(xsBFrameB)){
						
								xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
								XFrameA <- rbind(xsAFrameA, xsBFrameA)
								z5 <- matrix(0, n, ncol(XFrameA))
								z5[domains == "a",] <- XFrameA[domains == "a",]
								z5[domains == "ab",] <- XFrameA[domains == "ab",]/eta_0

								z <- as.matrix (cbind(z1, z2, z3, z4, z5))
								zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0, XA / Nhat_P))
								Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)	
							}
							else{
						
								xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
								XFrameA <- rbind(xsAFrameA, xsBFrameA)
								xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
								XFrameB <- rbind(xsAFrameB, xsBFrameB)
								z5 <- matrix(0, n, ncol(XFrameA))
								z5[domains == "a",] <- XFrameA[domains == "a",]
								z5[domains == "ab",] <- XFrameA[domains == "ab",]/eta_0
								z6 <- matrix(0, n, ncol(XFrameB))
								z6[domains == "b",] <- XFrameB[domains == "b",]
								z6[domains == "ba",] <- XFrameB[domains == "ba",]/(1 - eta_0)

								z <- as.matrix (cbind(z1, z2, z3, z4, z5, z6))
								zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0, XA / Nhat_P, XB / Nhat_P))
								Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)
							}
						}

						lambda <- Lag2 (z, di, zmean)
						phat <- as.vector(di_star / as.vector((rep(1,n) + crossprod(t(z - Zmean), lambda))))
						phat_star <- as.vector(phat / tapply(phat, domains, sum)[domains])

						Yhatmean_a <- sum (phat_star * Domains(sample[,k], domains, "a"))
						Yhatmean_ab <- sum (phat_star * Domains(sample[,k], domains, "ab"))
						Yhatmean_b <- sum (phat_star * Domains(sample[,k], domains, "b"))
						Yhatmean_ba <- sum (phat_star * Domains(sample[,k], domains, "ba"))

						meandom[,k] <- c(Yhatmean_a, Yhatmean_ab, Yhatmean_b, Yhatmean_ba)
						totdom[,k] <- meandom[,k] * domain_size_estimation

						mean_estimation <- Wa * Yhatmean_a + Wab_eta_0 * Yhatmean_ab + Wba_eta_0 * Yhatmean_ba + Wb * Yhatmean_b
						total_estimation <- Nhat_P * mean_estimation
						est[,k] <- c(total_estimation, mean_estimation)

						if (!is.null(conf_level)) {

							Di_star <- diag(di_star)

							suma1_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), z[domains == "a",] - Zmean[domains == "a",])
							suma2_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), sample[domains == "a",k] - mean_estimation)
							suma1_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), z[domains == "ab",] - Zmean[domains == "ab",])
							suma2_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), sample[domains == "ab",k] - mean_estimation)
							suma1_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), z[domains == "b",] - Zmean[domains == "b",])
							suma2_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), sample[domains == "b",k] - mean_estimation)
							suma1_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), z[domains == "ba",] - Zmean[domains == "ba",])
							suma2_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), sample[domains == "ba",k] - mean_estimation)
							Bhat <- solve(Wa * suma1_a + Wab_eta_0 * suma1_ab + Wba_eta_0 * suma1_ba + Wb * suma1_b, Wa * suma2_a + Wab_eta_0 * suma2_ab + Wba_eta_0 * suma2_ba + Wb * suma2_b)

							rhat <- sample[,k] - crossprod(t(z), Bhat)
							Rhatmean_a <- sum (phat_star[domains == "a"] * rhat[domains == "a"])
							Rhatmean_ab <- sum (phat_star[domains == "ab"] * rhat[domains == "ab"])
							Rhatmean_b <- sum (phat_star[domains == "b"] * rhat[domains == "b"])
							Rhatmean_ba <- sum (phat_star[domains == "ba"] * rhat[domains == "ba"])
							r_tilde <- (rhat - Rhatmean_a) * delta_a + (eta_0 * (rhat - Rhatmean_ab)) * delta_ab + (rhat - Rhatmean_b) * delta_b + ((1 - eta_0) * (rhat - Rhatmean_ba)) * delta_ba
							Vhat_A_Ar <- 1 / Nhat_P^2 * varest (r_tilde[domains == "a" | domains == "ab"], pik = pi_A)
							Vhat_B_Br <- 1 / Nhat_P^2 * varest (r_tilde[domains == "b" | domains == "ba"], pik = pi_B)

							rhatmean <- mean_estimation - crossprod(Bhat, zmean)
							Rhatmean <- rep (rhatmean, n)
							suma3_a <- crossprod(di_star[domains == "a"], (rhat[domains == "a"] - Rhatmean[domains == "a"])^2)
							suma3_ab <- crossprod(di_star[domains == "ab"], (rhat[domains == "ab"] - Rhatmean[domains == "ab"])^2)
							suma3_b <- crossprod(di_star[domains == "b"], (rhat[domains == "b"] - Rhatmean[domains == "b"])^2)
							suma3_ba <- crossprod(di_star[domains == "ba"], (rhat[domains == "ba"] - Rhatmean[domains == "ba"])^2)
							deff_P <- (Vhat_A_Ar + Vhat_B_Br) / (1 / n * (Wa * suma3_a + Wab_eta_0 * suma3_ab + Wba_eta_0 * suma3_ba + Wb * suma3_b))
							conf_interv <- PELConfInt (conf_level, sample[,k], di_star, mean_estimation, n/deff_P, di_star)
							upper <- conf_interv[1]
							lower <- conf_interv[2]
							interv[,k] <- c(total_estimation, lower * Nhat_P, upper * Nhat_P, mean_estimation, lower, upper)	
						}
					}
					else {
			
						N <- N_A + N_B - N_ab

						Wa <- (N_A - N_ab) / N
						Wab_eta_0 <- eta_0 * N_ab / N
						Wba_eta_0 <- (1 - eta_0) * N_ab / N
						Wb <- (N_B - N_ab) / N
						domain_size_estimation <- c(N_A - N_ab, N_ab, N_B - N_ab,  N_ab)

						if (is.null(xsAFrameA)){
					
							xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
							XFrameB <- rbind(xsAFrameB, xsBFrameB)
							z5 <- matrix(0, n, ncol(XFrameB))
							z5[domains == "b",] <- XFrameB[domains == "b",]
							z5[domains == "ba",] <- XFrameB[domains == "ba",]/(1 - eta_0)

							z <- as.matrix (cbind(z1, z2, z3, z4, z5))
							zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0, XB / N))
							Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)
						}
						else {
				
							if (is.null(xsBFrameB)){
					
								xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
								XFrameA <- rbind(xsAFrameA, xsBFrameA)
								z5 <- matrix(0, n, ncol(XFrameA))
								z5[domains == "a",] <- XFrameA[domains == "a",]
								z5[domains == "ab",] <- XFrameA[domains == "ab",]/eta_0
	
								z <- as.matrix (cbind(z1, z2, z3, z4, z5))
								zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0, XA / N))
								Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)	
							}
							else{
						
								xsAFrameA <- as.matrix(xsAFrameA); xsBFrameA <- as.matrix(xsBFrameA)
								XFrameA <- rbind(xsAFrameA, xsBFrameA)
								xsAFrameB <- as.matrix(xsAFrameB); xsBFrameB <- as.matrix(xsBFrameB)
								XFrameB <- rbind(xsAFrameB, xsBFrameB)
								z5 <- matrix(0, n, ncol(XFrameA))
								z5[domains == "a",] <- XFrameA[domains == "a",]
								z5[domains == "ab",] <- XFrameA[domains == "ab",]/eta_0
								z6 <- matrix(0, n, ncol(XFrameB))
								z6[domains == "b",] <- XFrameB[domains == "b",]
								z6[domains == "ba",] <- XFrameB[domains == "ba",]/(1 - eta_0)

								z <- as.matrix (cbind(z1, z2, z3, z4, z5, z6))
								zmean <- as.matrix (c(Wa, Wab_eta_0, Wba_eta_0, 0, XA / N, XB / N))
								Zmean <- matrix (zmean, nrow = nrow(z), ncol = ncol(z), byrow = T)
							}
						}

						lambda <- Lag2 (z, di, zmean)
						phat <- as.vector(di_star / as.vector((rep(1,n) + crossprod(t(z - Zmean), lambda))))
						phat_star <- as.vector(phat / tapply(phat, domains, sum)[domains])

						Yhatmean_a <- sum (phat_star * Domains(sample[,k], domains, "a"))
						Yhatmean_ab <- sum (phat_star * Domains(sample[,k], domains, "ab"))
						Yhatmean_b <- sum (phat_star * Domains(sample[,k], domains, "b"))
						Yhatmean_ba <- sum (phat_star * Domains(sample[,k], domains, "ba"))

						meandom[,k] <- c(Yhatmean_a, Yhatmean_ab, Yhatmean_b, Yhatmean_ba)
						totdom[,k] <- meandom[,k] * domain_size_estimation

						mean_estimation <- Wa * Yhatmean_a + Wab_eta_0 * Yhatmean_ab + Wba_eta_0 * Yhatmean_ba + Wb * Yhatmean_b
						total_estimation <- N * mean_estimation
						est[,k] <- c(total_estimation, mean_estimation)

						if (!is.null(conf_level)) {

							Di_star <- diag(di_star)

							suma1_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), z[domains == "a",] - Zmean[domains == "a",])
							suma2_a <- crossprod(t(crossprod(z[domains == "a",] - Zmean[domains == "a",], Di_star[domains == "a", domains == "a"])), sample[domains == "a",k] - mean_estimation)
							suma1_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), z[domains == "ab",] - Zmean[domains == "ab",])
							suma2_ab <- crossprod(t(crossprod(z[domains == "ab",] - Zmean[domains == "ab",], Di_star[domains == "ab", domains == "ab"])), sample[domains == "ab",k] - mean_estimation)
							suma1_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), z[domains == "b",] - Zmean[domains == "b",])
							suma2_b <- crossprod(t(crossprod(z[domains == "b",] - Zmean[domains == "b",], Di_star[domains == "b", domains == "b"])), sample[domains == "b",k] - mean_estimation)
							suma1_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), z[domains == "ba",] - Zmean[domains == "ba",])
							suma2_ba <- crossprod(t(crossprod(z[domains == "ba",] - Zmean[domains == "ba",], Di_star[domains == "ba", domains == "ba"])), sample[domains == "ba",k] - mean_estimation)
							Bhat <- solve(Wa * suma1_a + Wab_eta_0 * suma1_ab + Wba_eta_0 * suma1_ba + Wb * suma1_b, Wa * suma2_a + Wab_eta_0 * suma2_ab + Wba_eta_0 * suma2_ba + Wb * suma2_b)

							rhat <- sample[,k] - crossprod(t(z), Bhat)
							Rhatmean_a <- sum (phat_star[domains == "a"] * rhat[domains == "a"])
							Rhatmean_ab <- sum (phat_star[domains == "ab"] * rhat[domains == "ab"])
							Rhatmean_b <- sum (phat_star[domains == "b"] * rhat[domains == "b"])
							Rhatmean_ba <- sum (phat_star[domains == "ba"] * rhat[domains == "ba"])
							r_tilde <- (rhat - Rhatmean_a) * delta_a + (eta_0 * (rhat - Rhatmean_ab)) * delta_ab + (rhat - Rhatmean_b) * delta_b + ((1 - eta_0) * (rhat - Rhatmean_ba)) * delta_ba
							Vhat_A_Ar <- 1 / N^2 * varest (r_tilde[domains == "a" | domains == "ab"], pik = pi_A)
							Vhat_B_Br <- 1 / N^2 * varest (r_tilde[domains == "b" | domains == "ba"], pik = pi_B)

							rhatmean <- mean_estimation - crossprod(Bhat, zmean)
							Rhatmean <- rep (rhatmean, n)
							suma3_a <- crossprod(di_star[domains == "a"], (rhat[domains == "a"] - Rhatmean[domains == "a"])^2)
							suma3_ab <- crossprod(di_star[domains == "ab"], (rhat[domains == "ab"] - Rhatmean[domains == "ab"])^2)
							suma3_b <- crossprod(di_star[domains == "b"], (rhat[domains == "b"] - Rhatmean[domains == "b"])^2)
							suma3_ba <- crossprod(di_star[domains == "ba"], (rhat[domains == "ba"] - Rhatmean[domains == "ba"])^2)
							deff_P <- (Vhat_A_Ar + Vhat_B_Br) / (1 / n * (Wa * suma3_a + Wab_eta_0 * suma3_ab + Wba_eta_0 * suma3_ba + Wb * suma3_b))

							conf_interv <- PELConfInt (conf_level, sample[,k], di_star, mean_estimation, n/deff_P, di_star)
							upper <- conf_interv[1]
							lower <- conf_interv[2]
							interv[,k] <- c(total_estimation, lower * N, upper * N, mean_estimation, lower, upper)	
						}		
					}
				}
			}	
		}
		else
			stop("Invalid option: Probability vector in one frame and probability matrix in the other frame. Type of probabilities structures must match.")
	}
   	results = list(Call = cl, Est = est, TotDomEst = totdom, MeanDomEst = meandom, Param = par, ConfInt = interv)
   	class(results) = "EstimatorDF"
   	attr(results, "attributesDF") = conf_level
   	return(results)
}