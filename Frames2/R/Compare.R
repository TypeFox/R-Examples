#' @name Compare
#' @aliases Compare
#' @title Summary of estimators
#' 
#' @description Returns all possible estimators that can be computed according to the information provided
#' 
#' @usage Compare(ysA, ysB, pi_A, pi_B, domains_A, domains_B, pik_ab_B = NULL, pik_ba_A = NULL, 
#' N_A = NULL, N_B = NULL, N_ab = NULL, xsAFrameA = NULL, xsBFrameA = NULL,  
#' xsAFrameB = NULL, xsBFrameB = NULL, XA = NULL, XB = NULL, met = "linear", 
#' conf_level = NULL)
#' @param ysA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{c} containing information about variable(s) of interest from \eqn{s_A}.
#' @param ysB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{c} containing information about variable(s) of interest from \eqn{s_B}.
#' @param pi_A A numeric vector of length \eqn{n_A} or a square numeric matrix of dimension \eqn{n_A} containing first order or first and second order inclusion probabilities for units included in \eqn{s_A}.
#' @param pi_B A numeric vector of length \eqn{n_B} or a square numeric matrix of dimension \eqn{n_B} containing first order or first and second order inclusion probabilities for units included in \eqn{s_B}.
#' @param domains_A A character vector of length \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domains_B A character vector of length \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".
#' @param pik_ab_B (Optional) A numeric vector of size \eqn{n_A} containing first order inclusion probabilities according to sampling desing in frame B for units belonging 
#'  to overlap domain that have been selected in \eqn{s_A}.
#' @param pik_ba_A (Optional) A numeric vector of size \eqn{n_B} containing first order inclusion probabilities according to sampling desing in frame A for units belonging 
#'  to overlap domain that have been selected in \eqn{s_B}.
#' @param N_A (Optional) A numeric value indicating the size of frame A. 
#' @param N_B (Optional) A numeric value indicating the size of frame B.
#' @param N_ab (Optional) A numeric value indicating the size of the overlap domain.
#' @param xsAFrameA (Optional) A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, containing auxiliary information in frame A for units included in \eqn{s_A}.
#' @param xsBFrameA (Optional) A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, containing auxiliary information in frame A for units included in \eqn{s_B}. For units in domain \eqn{b}, these values are 0.
#' @param xsAFrameB (Optional) A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, containing auxiliary information in frame B for units included in \eqn{s_A}. For units in domain \eqn{a}, these values are 0.
#' @param xsBFrameB (Optional) A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, containing auxiliary information in frame B for units included in \eqn{s_B}.
#' @param XA (Optional) A numeric value or vector of length \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, indicating the population totals for the auxiliary variables considered in frame A.
#' @param XB (Optional) A numeric value or vector of length \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, indicating the population totals for the auxiliary variables considered in frame B.
#' @param met (Optional) A character vector indicating the distance that must be used in calibration process. Possible values are "linear", "raking" and "logit". Default is "linear".
#' @param conf_level (Optional) A numeric value indicating the confidence level for the confidence intervals, if desired.
#' @examples
#' data(DatA)
#' data(DatB)
#' data(PiklA)
#' data(PiklB)
#' 
#' Compare(DatA$Feed, DatB$Feed, PiklA, PiklB, DatA$Domain, DatB$Domain)
#' @export
Compare = function (ysA, ysB, pi_A, pi_B, domains_A, domains_B, pik_ab_B = NULL, pik_ba_A = NULL, N_A = NULL, N_B = NULL, N_ab = NULL, xsAFrameA = NULL, xsBFrameA = NULL, xsAFrameB = NULL, xsBFrameB = NULL, XA = NULL, XB = NULL, met = "linear", conf_level = NULL)
{

	EstPEL <- PEL(ysA, ysB, pi_A, pi_B, domains_A, domains_B, N_A = N_A, N_B = N_B, N_ab = N_ab, xsAFrameA = xsAFrameA, xsBFrameA = xsBFrameA, xsAFrameB = xsAFrameB, xsBFrameB = xsBFrameB, XA = XA, XB = XB, conf_level = conf_level)
	EstCalDF <- CalDF(ysA, ysB, pi_A, pi_B, domains_A, domains_B, N_A = N_A, N_B = N_B, N_ab = N_ab, xsAFrameA = xsAFrameA, xsBFrameA = xsBFrameA, xsAFrameB = xsAFrameB, xsBFrameB = xsBFrameB, XA = XA, XB = XB, met = met, conf_level = conf_level)
	estimators <- list(PEL = EstPEL, Calibration_DF = EstCalDF)	
	
	if (!is.null(pik_ab_B) & !is.null(pik_ba_A)){

		EstCalSF <- CalSF(ysA, ysB, pi_A, pi_B, domains_A, domains_B, pik_ab_B = pik_ab_B, pik_ba_A = pik_ba_A, N_A = N_A, N_B = N_B, N_ab = N_ab, xsAFrameA = xsAFrameA, xsBFrameA = xsBFrameA, xsAFrameB = xsAFrameB, xsBFrameB = xsBFrameB, XA = XA, XB = XB, met = met, conf_level = conf_level)
		estimators <- list(PEL = EstPEL, Calibration_DF = EstCalDF, Calibration_SF = EstCalSF)
	}

	if (!is.null(N_A) & !is.null(N_B)){

		EstPML <- PML(ysA, ysB, pi_A, pi_B, domains_A, domains_B, N_A = N_A, N_B = N_B, conf_level = conf_level)
		estimators <- list(PML = EstPML, PEL = EstPEL, Calibration_DF = EstCalDF)

		if (!is.null(pik_ab_B) & !is.null(pik_ba_A)){
			EstSFRR <- SFRR(ysA, ysB, pi_A, pi_B, domains_A, domains_B, pik_ab_B = pik_ab_B, pik_ba_A = pik_ba_A, N_A = N_A, N_B = N_B, conf_level = conf_level)
			estimators <- list(PML = EstPML, PEL = EstPEL, Calibration_DF = EstCalDF, SFRR = EstSFRR, Calibration_SF = EstCalSF)
		 }
	}
	else{

		EstHartley <- Hartley(ysA, ysB, pi_A, pi_B, domains_A, domains_B, conf_level = conf_level)
		EstFB <- FB(ysA, ysB, pi_A, pi_B, domains_A, domains_B, conf_level = conf_level)
		estimators <- list(Hartley = EstHartley, FullerBurmeister = EstFB, PEL = EstPEL, Calibration_DF = EstCalDF)

		if (!is.null(pik_ab_B) & !is.null(pik_ba_A)){
			EstBKA <- BKA(ysA, ysB, pi_A, pi_B, domains_A, domains_B, pik_ab_B = pik_ab_B, pik_ba_A = pik_ba_A, conf_level = conf_level)
			estimators <- list(Hartley = EstHartley, FullerBurmeister = EstFB, PEL = EstPEL, Calibration_DF = EstCalDF, BKA = EstBKA, Calibration_SF = EstCalSF)
		}

	}
   	return(estimators)
}