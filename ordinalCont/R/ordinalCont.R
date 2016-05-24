#' ordinalCont-package
#'
#' @name ordinalCont-package
#' @docType package
#' @details Ordinal regression analysis is a convenient tool for analyzing ordinal response variables 
#' in the presence of covariates. We extend this methodology to the case of continuous self-rating 
#' scales such as the Visual Analog Scale (VAS) used in pain assessment, or the Linear Analog 
#' Self-Assessment (LASA) scales in quality of life studies. Subjects are
#' typically given a linear scale of 100 mm and asked to put a mark where they perceive
#' themselves. These scales  measure subjects' 
#' perception of an intangible quantity, and cannot be handled as ratio variables because of their 
#' inherent nonlinearity.  Instead we treat them as ordinal variables, measured on a continuous scale. We 
#' express  the likelihood in terms of a function (the ``g function'')
#'  connecting the  
#' scale with an underlying continuous latent  variable. In the current version the g function 
#' is taken as 
#' the generalized logistic function (Richards 1959). This has 3 parameters: 
#' \code{M}, the offset, \code{B}, the slope, and \code{T}, the symmetry of the curve.
#' The link function is the inverse of the CDF of the assumed underlying distribution of the 
#' latent variable. Currently 
#' the logit link, which corresponds to a standard logistic distribution, is implemented. 
#' (This implies a proportional odds model.)  The likelihood is 
#' maximized using \code{optim {stats}} with a quasi-Newton method (\code{"BFGS"}). Fixed-effects models are implemented in the function \code{\link{ocm}}, and mixed models in  \code{\link{ocmm}}. 
#' @references   Manuguerra M, Heller GZ (2010). Ordinal Regression Models for Continuous 
#'  Scales, \emph{The International Journal of Biostatistics}: 6(1), Article 14.
#'@references Richards, F. (1959). A flexible growth function for empirical use, 
#' \emph{Journal of Experimental Botany}, 10, 290-301.
#' @author Maurizio Manuguerra, Gillian Heller
NULL


#' @title ANZ0001 trial
#' 
#' @details  The ANZ0001 trial, conducted by the ANZ Breast Cancer Trials Group, is an unblinded, multi-centre, randomized trial with three chemotherapy treatment arms, concluded in 2005 (Stockler et al 2007). 
#' Health-related quality of life measures (Overall quality of life, Physical Well-Being, Mood, Pain, Nausea and Vomiting, Appetite) are assessed at each chemotherapy treatment cycle, from randomization until disease progression, when treatment is interrupted. 
#' The treatments Intermittent Capecitabine (IC) and Continuous Capecitabine (CC) are compared with the standard combination treatment CMF, each with its own protocol. 
#' There is no maximum duration of treatment, but it is interrupted on disease progression, or when patient intolerance or unacceptable toxicity are recorded.
#' The data set is extracted from the ANZ0001 trial and contains information from 292 patients with complete quality of life measurements.
#'
#' The variables are as follows:
#'
#' \tabular{ll}{\code{randno}\tab{patient ID number}\cr
#'\code{cycleno}\tab{chemotherapy cycle number}\cr
#'\code{age}\tab{age of patient at entry to study}\cr
#'\code{bsa}\tab{Body Surface Area (m\eqn{^2})}\cr
#'\code{treatment}\tab  treatment received by  patient (1,2,3)\cr
#'\code{overall}\tab Overall quality of life as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#'\code{phys}\tab Physical Well-Being as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#'\code{mood}\tab Mood as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#'\code{pain}\tab Pain as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#'\code{nausvom}\tab Nausea and Vomiting as recorded by the patient on a LASA scale,  normalized to  (0, 1)\cr
#'\code{appetite}\tab Appetite as recorded by the patient on a LASA scale,  normalized to  (0, 1)
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ANZ0001
#' @usage data(ANZ0001)
#' @references Stockler, M., T. Sourjina, P. Grimison, V. Gebski, M. Byrne, V. Harvey, P. Francis et al. ``A randomized trial of capecitabine (C) given intermittently (IC) rather than continuously (CC) compared to classical CMF as first-line chemotherapy for advanced breast cancer (ABC).'' In \emph{ASCO Annual Meeting Proceedings}, vol. 25, no. 18_suppl, p. 1031. 2007.
#' @format A data frame with 2473 rows and 11 variables
NULL

