#' @name MLCDF
#' @aliases MLCDF
#' @title Multinomial logistic calibration estimator under dual frame approach with auxiliary information from each frame
#' 
#' @description Produces estimates for class totals and proportions using multinomial logistic regression from survey data obtained
#'  from a dual frame sampling design using a model calibrated dual frame approach with a possibly different set of auxiliary variables for each frame.
#'  Confidence intervals are also computed, if required.
#' 
#' @usage MLCDF (ysA, ysB, pik_A, pik_B, domains_A, domains_B, xsA, xsB, xA, xB, ind_samA, 
#'  ind_samB, ind_domA, ind_domB, N, N_ab = NULL, met = "linear", conf_level = NULL)
#' @param ysA A data frame containing information about one or more factors, each one of dimension \eqn{n_A}, collected from \eqn{s_A}.
#' @param ysB A data frame containing information about one or more factors, each one of dimension \eqn{n_B}, collected from \eqn{s_B}.
#' @param pik_A A numeric vector of length \eqn{n_A} containing first order inclusion probabilities for units included in \eqn{s_A}.
#' @param pik_B A numeric vector of length \eqn{n_B} containing first order inclusion probabilities for units included in \eqn{s_B}.
#' @param domains_A A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domains_B A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".
#' @param xsA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, containing auxiliary information in frame A for units included in \eqn{s_A}.
#' @param xsB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, containing auxiliary information in frame B for units included in \eqn{s_B}.
#' @param xA A numeric vector or length \eqn{N_A} or a numeric matrix or data frame of dimensions \eqn{N_A} x \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, containing auxiliary information for the units in frame A.
#' @param xB A numeric vector or length \eqn{N_B} or a numeric matrix or data frame of dimensions \eqn{N_B} x \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, containing auxiliary information for the units in frame B.
#' @param ind_samA A numeric vector of length \eqn{n_A} containing the identificators of units of the frame A (from 1 to \eqn{N_A}) that belongs to \eqn{s_A}.
#' @param ind_samB A numeric vector of length \eqn{n_B} containing the identificators of units of the frame B (from 1 to \eqn{N_B}) that belongs to \eqn{s_B}.
#' @param ind_domA A character vector of length \eqn{N_A} indicating the domain each unit from frame A belongs to. Possible values are "a" and "ab".
#' @param ind_domB A character vector of length \eqn{N_B} indicating the domain each unit from frame B belongs to. Possible values are "b" and "ba".
#' @param N A numeric value indicating the size of the population.
#' @param N_ab (Optional) A numeric value indicating the size of the overlap domain
#' @param met (Optional) A character vector indicating the distance that must be used in calibration process. Possible values are "linear", "raking" and "logit". Default is "linear".
#' @param conf_level (Optional) A numeric value indicating the confidence level for the confidence intervals, if desired.
#' @details Multinomial logistic calibration estimator in dual frame using auxiliary information from each frame for a proportion is given by
#'  \deqn{\hat{P}_{MLCi}^{DF} = \frac{1}{N} \left(\sum_{k \in s_A \cup s_B} w_k^{\circ} z_{ki}\right), \hspace{0.3cm} i = 1,...,m}
#'  with \eqn{m} the number of categories of the response variable, \eqn{z_i} the indicator variable for the i-th category of the response variable,
#'  and \eqn{w^{\circ}} calibration weights which are calculated having into account a different set of constraints, depending on the case. For instance, if \eqn{N_A, N_B} and \eqn{N_{ab}} are known, calibration constraints are
#'  \deqn{\sum_{k \in s_a}w_k^{\circ} = N_a, \sum_{k \in s_{ab}}w_k^{\circ} = \eta N_{ab}, \sum_{k \in s_{ba}}w_k^{\circ} = (1 - \eta) N_{ab}\sum_{k \in s_{b}}w_k^{\circ} = N_{b},} \deqn{\sum_{k \in s_A}w_k^\circ p_{ki}^A = \sum_{k \in U_a} p_{ki}^A + \eta \sum_{k \in U_{ab}} p_{ki}^A}
#'  and \deqn{\sum_{k \in s_B}w_k^\circ p_{ki}^B = \sum_{k \in U_b} p_{ki}^B + (1 - \eta) \sum_{k \in U_{ba}} p_{ki}^B}
#' with \eqn{\eta \in (0,1)} and \deqn{p_{ki}^A = \frac{exp(x_k^{'}\beta_i^A)}{\sum_{r=1}^m exp(x_k^{'}\beta_r^A)},}
#' being \eqn{\beta_i^A} the maximum likelihood parameters of the multinomial logistic model considering original design weights \eqn{d^A}. \eqn{p_{ki}^B} can be defined similarly.
#' @return \code{MLCDF} returns an object of class "MultEstimatorDF" which is a list with, at least, the following components:
#'  \item{Call}{the matched call.}
#'  \item{Est}{class frequencies and proportions estimations for main variable(s).}
#' @references Molina, D., Rueda, M., Arcos, A. and Ranalli, M. G. (2015)
#'  \emph{Multinomial logistic estimation in dual frame surveys}
#'  Statistics and Operations Research Transactions (SORT). To be printed.
#' @seealso \code{\link{JackMLCDF}}
#' @examples
#' data(DatMA)
#' data(DatMB)
#' data(DatPopM) 
#'
#' N <- nrow(DatPopM)
#' levels(DatPopM$Domain) <- c(levels(DatPopM$Domain), "ba")
#' DatPopMA <- subset(DatPopM, DatPopM$Domain == "a" | DatPopM$Domain == "ab", stringAsFactors = FALSE)
#' DatPopMB <- subset(DatPopM, DatPopM$Domain == "b" | DatPopM$Domain == "ab", stringAsFactors = FALSE)
#' DatPopMB[DatPopMB$Domain == "ab",]$Domain <- "ba"
#'
#' #Let calculate proportions of categories of variable Prog using MLCDF estimator
#' #using Read as auxiliary variable
#' MLCDF(DatMA$Prog, DatMB$Prog, DatMA$ProbA, DatMB$ProbB, DatMA$Domain, DatMB$Domain, 
#' DatMA$Read, DatMB$Read, DatPopMA$Read, DatPopMB$Read, DatMA$Id_Frame, DatMB$Id_Frame, 
#' DatPopMA$Domain, DatPopMB$Domain, N)
#'
#' #Let obtain 95% confidence intervals together with the estimations
#' MLCDF(DatMA$Prog, DatMB$Prog, DatMA$ProbA, DatMB$ProbB, DatMA$Domain, DatMB$Domain, 
#' DatMA$Read, DatMB$Read, DatPopMA$Read, DatPopMB$Read, DatMA$Id_Frame, DatMB$Id_Frame, 
#' DatPopMA$Domain, DatPopMB$Domain, N, conf_level = 0.95)
#' @export
MLCDF = function (ysA, ysB, pik_A, pik_B, domains_A, domains_B, xsA, xsB, xA, xB, ind_samA, ind_samB, ind_domA, ind_domB, N, N_ab = NULL, met = "linear", conf_level = NULL){

	ysA <- as.data.frame(ysA)
	ysB <- as.data.frame(ysB)
	xsA <- as.matrix(xsA)
	xsB <- as.matrix(xsB)
	xA <- as.matrix(xA)
	xB <- as.matrix(xB)

	if (any(is.na(ysA)))
		stop("There are missing values in sample from frame A.")
	if (any(is.na(ysB)))
		stop("There are missing values in sample from frame B.")
	if (any(is.na(pik_A)))
		stop("There are missing values in pikl from frame A.")
	if (any(is.na(pik_B)))
		stop("There are missing values in pikl from frame B.")
	if (any(is.na(domains_A)))
		stop("There are missing values in domains from frame A.")
	if (any(is.na(domains_B)))
		stop("There are missing values in domains from frame B.")
	if (nrow(ysA) != length(pik_A) | nrow(ysA) != length(domains_A) | length(domains_A) != length(pik_A))
		stop("Arguments from frame A have different sizes.")
	if (nrow(ysB) != length(pik_B) | nrow(ysB) != length(domains_B) | length(domains_B) != length(pik_B))
		stop("Arguments from frame B have different sizes.")
	if (ncol(ysA) != ncol(ysB))
		stop("Number of variables does not match.")
	if (length(which(domains_A == "a")) + length(which(domains_A == "ab")) != length(domains_A))
		stop("Domains from frame A are not correct.")
	if (length(which(domains_B == "b")) + length(which(domains_B == "ba")) != length(domains_B))
		stop("Domains from frame B are not correct.")

	cl <- match.call()

	estimations <- list()
	interv <- list()
	c <- ncol(ysA)
	R <- ncol(xA)
	n_A <- nrow(ysA)
	n_B <- nrow(ysB)
	n <- n_A + n_B
	N_A <- nrow(xA)
	N_B <- nrow(xB)

	ones_ab_A <- Domains (rep (1, n_A), domains_A, "ab")
	ones_ab_B <- Domains (rep (1, n_B), domains_B, "ba")
	Vhat_Nhat_ab_A <- varest(ones_ab_A, pik = pik_A)
	Vhat_Nhat_ab_B <- varest(ones_ab_B, pik = pik_B)

	eta_0 <- Vhat_Nhat_ab_B / (Vhat_Nhat_ab_A + Vhat_Nhat_ab_B) 
	
	domains <- factor(c(as.character(domains_A), as.character(domains_B)))
	delta_a <- Domains (rep (1, n), domains, "a")
	delta_ab <- Domains (rep (1, n), domains, "ab")
	delta_b <- Domains (rep (1, n), domains, "b")
	delta_ba <- Domains (rep (1, n), domains, "ba")

	for (k in 1:c){

		ys <- factor(c(as.character(ysA[,k]),as.character(ysB[,k])))
		y_sA <- as.factor(ysA[,k])
		y_sB <- as.factor(ysB[,k])

		lev <- sort(levels(ys))
		levA <- sort(levels(y_sA))
		levB <- sort(levels(y_sB))
		m <- length(lev)
		mA <- length(levA)
		mB <- length(levB)
	
		mat <- matrix (NA, 2, m)
		rownames(mat) <- c("Class Tot.", "Prop.")
		colnames(mat) <- lev

		z <- disjunctive(ys)
		pik <- c(pik_A, pik_B)
		d <- 1/pik

		modA <- multinom(formula = y_sA ~ 0 + xsA, weights = 1/pik_A, trace = FALSE)		
		beta_tilde_A <- rbind(rep(0, R), summary(modA)$coefficients)
		modB <- multinom(formula = y_sB ~ 0 + xsB, weights = 1/pik_B, trace = FALSE)		
		beta_tilde_B <- rbind(rep(0, R), summary(modB)$coefficients)

		denomA <- rowSums(exp(xA %*% t(beta_tilde_A)))
		denomB <- rowSums(exp(xB %*% t(beta_tilde_B)))

		pA <- exp (xA %*% t(beta_tilde_A)) / denomA
		pB <- exp (xB %*% t(beta_tilde_B)) / denomB

		if (mA < m){
			common <- which(lev %in% levA)
			pA2 <- matrix(0, N_A, m)
			pA2[,common] <- pA
			pA <- pA2
		}
		if (mB < m){
			common <- which(lev %in% levB)
			pB2 <- matrix(0, N_B, m)
			pB2[,common] <- pB 			
			pB <- pB2
		}
	
		psA <- pA[ind_samA,]
		psB <- pB[ind_samB,]
		ps <- rbind(psA, psB)

		pA[ind_domA == "ab",] <- eta_0 * pA[ind_domA == "ab",]
		pB[ind_domB == "ba",] <- (1 - eta_0) * pB[ind_domB == "ba",]

		if (is.null(N_ab)){

			Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, ps * (delta_a + delta_ab), ps * (delta_b + delta_ba))
			total <- c(N_A, N_B, colSums(pA), colSums(pB))
		}
		else {
			Xs <- cbind(delta_a, delta_ab, delta_ba, delta_b, ps * (delta_a + delta_ab), ps * (delta_b + delta_ba))
			total <- c(N_A - N_ab, eta_0 * N_ab, (1 - eta_0) * N_ab, N_B - N_ab, colSums(pA), colSums(pB))
		}

		g <- calib (Xs, d, total, method = met)

		mat[1,] <- colSums (g * d * z)
		mat[2,] <- 1/N * mat[1,]
		estimations[[k]] <- mat

		if (!is.null(conf_level)){

			alpha <- ginv(t(Xs) %*% diag(d) %*% Xs) %*% t(Xs) %*% diag(d) %*% z
			e <- z - Xs %*% alpha
			e <- e * d
			e[domains == "ab",] <- eta_0 * e[domains == "ab",]
			e[domains == "ba",] <- (1 - eta_0) * e[domains == "ba",]
			Vhat_AMLCDF <- apply(e, 2, var)
			Vhat_PMLCDF <- 1/N^2 * Vhat_AMLCDF

			interval <- matrix (NA, 6, m)
			rownames(interval) <- c("Class Tot.", "Lower Bound", "Upper Bound", "Prop.", "Lower Bound", "Upper Bound")
			colnames(interval) <- lev

			interval[1,] <- mat[1,]		
			interval[2,] <- mat[1,] + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_AMLCDF)
			interval[3,] <- mat[1,] - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_AMLCDF)
			interval[4,] <- mat[2,]
			interval[5,] <- mat[2,] + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_PMLCDF)
			interval[6,] <- mat[2,] - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_PMLCDF)
			interv[[k]] <- interval
		}
	}

	results = list(Call = cl, Est = estimations, ConfInt = interv)
	class(results) = "EstimatorMDF"
	attr(results, "attributesMDF") = conf_level
	return(results)
}