#' @name WeightsCalSF
#' @aliases WeightsCalSF
#' @title g-weights for the SF calibration estimator
#' 
#' @description Computes the g-weights for the SF calibration estimator. 
#' 
#' @usage WeightsCalSF(ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, 
#' N_A = NULL, N_B = NULL, N_ab = NULL, xsAFrameA = NULL, xsBFrameA = NULL, 
#' xsAFrameB = NULL, xsBFrameB = NULL, xsT = NULL, XA = NULL, XB = NULL, X = NULL, 
#' met = "linear")
#' @param ysA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{c} containing information about variable(s) of interest from \eqn{s_A}.
#' @param ysB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{c} containing information about variable(s) of interest from \eqn{s_A}.
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
#' @details Function provides g-weights in following scenarios:
#'  \itemize{
#'	\item There is not any additional auxiliary variable 
#'      \itemize{
#'      	\item \eqn{N_A, N_B} and \eqn{N_{ab}} unknown
#'              \item \eqn{N_{ab}} known and \eqn{N_A} and \eqn{N_B} unknown
#'              \item \eqn{N_A} and \eqn{N_B} known and \eqn{N_{ab}} unknown
#'              \item \eqn{N_A, N_B} and \eqn{N_{ab}} known
#'      }
#'      \item At least, one additional auxiliary variable is available 
#'      \itemize{
#'              \item \eqn{N_{ab}} known and \eqn{N_A} and \eqn{N_B} unknown
#'              \item \eqn{N_A} and \eqn{N_B} known and \eqn{N_{ab}} unknown
#'              \item \eqn{N_A, N_B} and \eqn{N_{ab}} known
#'      }
#'  }
#' @return A numeric vector containing the g-weights for the SF calibration estimator.
#' @references Ranalli, M. G., Arcos, A., Rueda, M. and Teodoro, A. (2013)
#'  \emph{Calibration estimationn in dual frame surveys}. arXiv:1312.0761 [stat.ME]
#' @references Deville, J. C., S\"arndal, C. E. (1992)
#'  \emph{Calibration estimators in survey sampling.}
#'  Journal of the American Statistical Association, 87, 376 - 382
#' @examples
#' data(DatA)
#' data(DatB)
#' data(PiklA)
#' data(PiklB)
#' 
#' #Let calculate g-weights for the SF calibration estimator for variable Clothing,
#' #without considering any auxiliary information
#' WeightsCalSF(DatA$Clo, DatB$Clo, PiklA, PiklB, DatA$ProbB, DatB$ProbA, 
#' DatA$Domain, DatB$Domain)
#' 
#' #Now, let calculate g-weights for the SF calibration estimator for variable Leisure
#' #when the frame sizes and the overlap domain size are known
#' WeightsCalSF(DatA$Lei, DatB$Lei, PiklA, PiklB, DatA$ProbB, DatB$ProbA, 
#' DatA$Domain, DatB$Domain, N_A = 1735, N_B = 1191, N_ab = 601)
#' 
#' #Finally, let calculate g-weights for the SF calibration estimator
#' #for variable Feeding, considering Income and Metres2 as auxiliary 
#' #variables and with frame sizes and overlap domain size known.
#' WeightsCalSF(DatA$Feed, DatB$Feed, PiklA, PiklB, DatA$ProbB, DatB$ProbA, 
#' DatA$Domain, DatB$Domain, N_A = 1735, N_B =  1191, N_ab = 601, xsAFrameA = DatA$Inc, 
#' xsBFrameA = DatB$Inc, xsAFrameB = DatA$M2, xsBFrameB = DatB$M2, 
#' XA = 4300260, XB = 176553)
#' @export
WeightsCalSF = function (ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, N_A = NULL, N_B = NULL, N_ab = NULL, xsAFrameA = NULL, xsBFrameA = NULL, xsAFrameB = NULL, xsBFrameB = NULL, xsT = NULL, XA = NULL, XB = NULL, X = NULL, met = "linear") {

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
		stop("A value for N_ab has been provided, but values for N_A or N_B are missing. This is not a possible option.")
	if ((is.null (xsAFrameA) & !is.null (xsBFrameA)) | (!is.null (xsAFrameA) & is.null (xsBFrameA)))
		stop("Auxiliary information from frame A is available only in one frame. This is not a possible option.")
	if ((is.null (xsAFrameB) & !is.null (xsBFrameB)) | (!is.null (xsAFrameB) & is.null (xsBFrameB)))
		stop("Auxiliary information from frame B is available only in one frame. This is not a possible option.")

	sample <- rbind(ysA, ysB)
	domains <- factor(c(as.character(domains_A), as.character(domains_B)))
	n_A <- nrow(ysA)
	n_B <- nrow(ysB)
	n <- n_A + n_B
	c <- ncol(sample)

	results <- matrix(NA, nrow = n, ncol = c)
	colnames (results) <- cnames

	delta_a <- Domains (rep (1, n), domains, "a")
	delta_ab <- Domains (rep (1, n), domains, "ab")
	delta_b <- Domains (rep (1, n), domains, "b")
	delta_ba <- Domains (rep (1, n), domains, "ba")

	ones_a_A <- Domains(rep (1, n_A), domains_A, "a")
	ones_ab_A <- Domains(rep (1, n_A), domains_A, "ab")
	ones_b_B <- Domains(rep (1, n_B), domains_B, "b")
	ones_ab_B <- Domains(rep (1, n_B), domains_B, "ba")

	if (!is.null(dim(drop(pi_A))) & !is.null(dim(drop(pi_B)))) {

		pik_A <- diag(pi_A)
		pik_B <- diag(pi_B)
		pik <- c(pik_A, pik_B)
		
		w_tilde_iS_A <- (1 / pik_A) * (domains_A == "a") + (1 / (pik_A + pik_ab_B)) * (domains_A == "ab")
		w_tilde_iS_B <- (1 / pik_B) * (domains_B == "b") + (1 / (pik_B + pik_ba_A)) * (domains_B == "ba")
		d <- c(w_tilde_iS_A, w_tilde_iS_B)

		Nhat_a_A <- HT (ones_a_A, pik_A)
		Nhat_ab_A <- HT (ones_ab_A, pik_A)
		Nhat_b_B <- HT (ones_b_B, pik_B)
		Nhat_ab_B <- HT (ones_ab_B, pik_B)

		for (k in 1:c) {

			if (is.null(xsAFrameA) & is.null(xsBFrameB) & is.null(xsT)) {

				if (is.null(N_ab)) {
				
					if (is.null(N_A) & is.null(N_B)) {

						Nhat_abS <- 0
						for (i in 1:n_A)
							if (domains_A[i] == "ab")
								Nhat_abS <- Nhat_abS + 1 / (1 / pik_A[i] + 1 / pik_ab_B[i])
				
						for (i in 1:n_B)
							if (domains_B[i] == "ba")
								Nhat_abS <- Nhat_abS + 1 / (1 / pik_B[i] + 1 / pik_ba_A[i])

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
			results[,k] <- calib (Xs, d, total, method = met)
		}
	}
	else {

		if (is.null(dim(drop(pi_A))) & is.null(dim(drop(pi_B)))){

			pik <- c(pi_A, pi_B)
			w_tilde_iS_A <- (1 / pi_A) * (domains_A == "a") + (1 / (pi_A + pik_ab_B)) * (domains_A == "ab")
			w_tilde_iS_B <- (1 / pi_B) * (domains_B == "b") + (1 / (pi_B + pik_ba_A)) * (domains_B == "ba")
			d <- c(w_tilde_iS_A, w_tilde_iS_B)

			Nhat_a_A <- HT (ones_a_A, pi_A)
			Nhat_ab_A <- HT (ones_ab_A, pi_A)
			Nhat_b_B <- HT (ones_b_B, pi_B)
			Nhat_ab_B <- HT (ones_ab_B, pi_B)

			for (k in 1:ncol(sample)) {

				if (is.null(xsAFrameA) & is.null(xsBFrameB) & is.null(xsT)) {

					if (is.null(N_ab)) {
				
						if (is.null(N_A) & is.null(N_B)) {

							Nhat_abS <- 0
							for (i in 1:n_A)
								if (domains_A[i] == "ab")
									Nhat_abS <- Nhat_abS + 1 / (1 / pi_A[i] + 1 / pik_ab_B[i])
				
							for (i in 1:n_B)
								if (domains_B[i] == "ba")
									Nhat_abS <- Nhat_abS + 1 / (1 / pi_B[i] + 1 / pik_ba_A[i])

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
				results[,k] <- calib (Xs, d, total, method = met)
			}
		}
		else
			stop("Invalid option: Probability vector in one frame and probability matrix in the other frame. Type of both structures must match.")
	}
	return (results)
}