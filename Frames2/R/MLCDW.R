#' @name MLCDW
#' @aliases MLCDW
#' @title Multinomial logistic calibration estimator under dual frame approach with auxiliary information from the whole population
#' 
#' @description Produces estimates for class totals and proportions using multinomial logistic regression from survey data obtained
#'  from a dual frame sampling design using a model calibrated dual frame approach with auxiliary information from the whole population. Confidence intervals are also computed, if required.
#' 
#' @usage MLCDW (ysA, ysB, pik_A, pik_B, domains_A, domains_B, xsA, xsB, x, ind_sam, N_A, 
#'  N_B, N_ab = NULL, met = "linear", conf_level = NULL)
#' @param ysA A data frame containing information about one or more factors, each one of dimension \eqn{n_A}, collected from \eqn{s_A}.
#' @param ysB A data frame containing information about one or more factors, each one of dimension \eqn{n_B}, collected from \eqn{s_B}.
#' @param pik_A A numeric vector of length \eqn{n_A} containing first order inclusion probabilities for units included in \eqn{s_A}.
#' @param pik_B A numeric vector of length \eqn{n_B} containing first order inclusion probabilities for units included in \eqn{s_B}.
#' @param domains_A A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domains_B A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".
#' @param xsA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{m}, with \eqn{m} the number of auxiliary variables, containing auxiliary information in frame A for units included in \eqn{s_A}.
#' @param xsB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{m}, with \eqn{m} the number of auxiliary variables, containing auxiliary information in frame B for units included in \eqn{s_B}.
#' @param x A numeric vector or length \eqn{N} or a numeric matrix or data frame of dimensions \eqn{N} x \eqn{m}, with \eqn{m} the number of auxiliary variables, containing auxiliary information for every unit in the population.
#' @param ind_sam A numeric vector of length \eqn{n = n_A + n_B} containing the identificators of units of the population (from 1 to \eqn{N}) that belongs to \eqn{s_A} or \eqn{s_B}
#' @param N_A A numeric value indicating the size of frame A
#' @param N_B A numeric value indicating the size of frame B
#' @param N_ab (Optional) A numeric value indicating the size of the overlap domain
#' @param met (Optional) A character vector indicating the distance that must be used in calibration process. Possible values are "linear", "raking" and "logit". Default is "linear".
#' @param conf_level (Optional) A numeric value indicating the confidence level for the confidence intervals, if desired.
#' @details Multinomial logistic calibration estimator in dual frame using auxiliary information from the whole population for a proportion is given by
#'  \deqn{\hat{P}_{MLCi}^{DW} = \frac{1}{N} \left(\sum_{k \in s_A \cup s_B} w_k^{\circ} z_{ki}\right), \hspace{0.3cm} i = 1,...,m}
#'  with \eqn{m} the number of categories of the response variable, \eqn{z_i} the indicator variable for the i-th category of the response variable,
#'  and \eqn{w^{\circ}} calibration weights which are calculated having into account a different set of constraints, depending on the case. For instance, if \eqn{N_A, N_B} and \eqn{N_{ab}} are known, calibration constraints are
#'  \deqn{\sum_{k \in s_a}w_k^{\circ} = N_a, \sum_{k \in s_{ab}}w_k^{\circ} = \eta N_{ab}, \sum_{k \in s_{ba}}w_k^{\circ} = (1 - \eta) N_{ab}, \sum_{k \in s_{b}}w_k^{\circ} = N_{b}} and \deqn{\sum_{k \in s_A \cup s_B}w_k^\circ p_{ki}^{\circ} = \sum_{k \in U} p_{ki}^\circ}
#' with \eqn{\eta \in (0,1)} and \deqn{p_{ki}^{\circ} = \frac{exp(x_k^{'}\beta_i^{\circ})}{\sum_{r=1}^m exp(x_k^{'}\beta_r^{\circ})},}
#' being \eqn{\beta_i^\circ} the maximum likelihood parameters of the multinomial logistic model considering weights \eqn{d_k^{\circ} =\left\{\begin{array}{lcc}
#'  d_k^A & \textrm{if } k \in a\\
#'  \eta d_k^A & \textrm{if } k \in ab\\
#'  (1 - \eta) d_k^B & \textrm{if } k \in ba \\
#'  d_k^B & \textrm{if } k \in b
#'  \end{array}
#'  \right.}.
#' @return \code{MLCDW} returns an object of class "MultEstimatorDF" which is a list with, at least, the following components:
#'  \item{Call}{the matched call.}
#'  \item{Est}{class frequencies and proportions estimations for main variable(s).}
#' @references Molina, D., Rueda, M., Arcos, A. and Ranalli, M. G. (2015)
#'  \emph{Multinomial logistic estimation in dual frame surveys}
#'  Statistics and Operations Research Transactions (SORT). To be printed. 
#' @seealso \code{\link{JackMLCDW}}
#' @examples
#' data(DatMA)
#' data(DatMB)
#' data(DatPopM) 
#'
#' IndSample <- c(DatMA$Id_Pop, DatMB$Id_Pop)
#' N_FrameA <- nrow(DatPopM[DatPopM$Domain == "a" | DatPopM$Domain == "ab",])
#' N_FrameB <- nrow(DatPopM[DatPopM$Domain == "b" | DatPopM$Domain == "ab",])
#' N_Domainab <- nrow(DatPopM[DatPopM$Domain == "ab",])
#	
#' #Let calculate proportions of categories of variable Prog using MLCDW estimator
#' #using Read as auxiliary variable
#' MLCDW(DatMA$Prog, DatMB$Prog, DatMA$ProbA, DatMB$ProbB, DatMA$Domain, DatMB$Domain, 
#' DatMA$Read, DatMB$Read, DatPopM$Read, IndSample, N_FrameA, N_FrameB)
#'
#' #Now, let suppose that the overlap domian size is known
#' MLCDW(DatMA$Prog, DatMB$Prog, DatMA$ProbA, DatMB$ProbB, DatMA$Domain, DatMB$Domain, 
#' DatMA$Read, DatMB$Read, DatPopM$Read, IndSample, N_FrameA, N_FrameB, N_Domainab)
#'
#' #Let obtain 95% confidence intervals together with the estimations
#' MLCDW(DatMA$Prog, DatMB$Prog, DatMA$ProbA, DatMB$ProbB, DatMA$Domain, DatMB$Domain, 
#' DatMA$Read, DatMB$Read, DatPopM$Read, IndSample, N_FrameA, N_FrameB, N_Domainab,
#' conf_level = 0.95)
#' @export
MLCDW = function (ysA, ysB, pik_A, pik_B, domains_A, domains_B, xsA, xsB, x, ind_sam, N_A, N_B, N_ab = NULL, met = "linear", conf_level = NULL){

	ysA <- as.data.frame(ysA)
	ysB <- as.data.frame(ysB)
	xsA <- as.matrix(xsA)
	xsB <- as.matrix(xsB)
	x <- as.matrix(x)

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
	if (is.null (N_ab) == "FALSE" & (is.null (N_A) == "TRUE" | is.null (N_B) == "TRUE"))
		stop("A value for N_ab has been provided, but values for N_A or N_B are missing. This is not a possible option.")

	cl <- match.call()

	estimations <- list()
	interv <- list()
	c <- ncol(ysA)
	xs <- rbind(xsA, xsB)
	N <- nrow(x)
	R <- ncol(x)
	n_A <- nrow(ysA)
	n_B <- nrow(ysB)
	n <- n_A + n_B

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

	pik <- c(pik_A, pik_B)
	dd <- 1/pik
	d <- dd * delta_a + dd * eta_0 * delta_ab + dd * (1-eta_0) * delta_ba + dd * delta_b

	for (k in 1:c){

		ys <- factor(c(as.character(ysA[,k]),as.character(ysB[,k])))

		lev <- levels(ys)
		m <- length(lev)

		mat <- matrix (NA, 2, m)
		rownames(mat) <- c("Class Tot.", "Prop.")
		colnames(mat) <- lev

		z <- disjunctive(ys)

		mod <- multinom(formula = ys ~ 0 + xs, weights = d, trace = FALSE)		
		beta_tilde <- rbind(rep(0, R), summary(mod)$coefficients)

		denom <- rowSums(exp(x %*% t(beta_tilde)))

		p <- exp (x %*% t(beta_tilde)) / denom

		if (is.null(N_ab)){

			Xs <- cbind(delta_a + delta_ab + delta_ba, delta_b + delta_ab + delta_ba, p[ind_sam,])
			total <- c(N_A, N_B, colSums(p))
		}
		else {
			Xs <- cbind(delta_a, delta_ab, delta_ba, delta_b, p[ind_sam,])
			total <- c(N_A - N_ab, eta_0 * N_ab, (1 - eta_0) * N_ab, N_B - N_ab, colSums(p))
		}

		g <- calib (Xs, d, total, method = met)

		mat[1,] <- colSums (g * d * z)
		mat[2,] <- 1/N * mat[1,]
		estimations[[k]] <- mat

		if (!is.null(conf_level)){

			alpha <- ginv(t(Xs) %*% diag(d) %*% Xs) %*% t(Xs) %*% diag(d) %*% z
			e <- z - Xs %*% alpha
			e <- e * d
			Vhat_AMLCDW <- apply(e, 2, var)
			Vhat_PMLCDW <- 1/N^2 * Vhat_AMLCDW

			interval <- matrix (NA, 6, m)
			rownames(interval) <- c("Class Tot.", "Lower Bound", "Upper Bound", "Prop.", "Lower Bound", "Upper Bound")
			colnames(interval) <- lev

			interval[1,] <- mat[1,]		
			interval[2,] <- mat[1,] + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_AMLCDW)
			interval[3,] <- mat[1,] - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_AMLCDW)
			interval[4,] <- mat[2,]
			interval[5,] <- mat[2,] + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_PMLCDW)
			interval[6,] <- mat[2,] - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_PMLCDW)
			interv[[k]] <- interval
		}
	}

	results = list(Call = cl, Est = estimations, ConfInt = interv)
	class(results) = "EstimatorMDF"
	attr(results, "attributesMDF") = conf_level
	return(results)
}