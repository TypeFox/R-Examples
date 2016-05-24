#' Design and trial settings used in the Lai, Lavori, Liao paper simulations
#'
#' A list of design and trial design settings used for analysis and simulations in
#' the Lai, Lavori, Liao paper displayed in Tables 1 and 2. The
#' elements of the list are the following
#' \describe{
#'   \item{trialParameters}{
#'     \describe{
#'      \item{N}{the sample size at each of three interim looks, the last being the final one;
#'               The length of this also determines the number of interim looks}
#'      \item{type1Error}{the overall type I error}
#'      \item{eps}{the fraction of type I error spent at each interim look}
#'      \item{type2Error}{the type II error desired}
#'     }
#'   }
#'   \item{scenarios}{
#'     A list of the 10 settings used in the simulations named \code{S0}, \code{S1}, ...,
#'     \code{S10} as in the paper, each with three elements
#'     \describe{
#'      \item{mean}{a \eqn{2\times J} matrix of means, the first row for the null setting,
#'                  the second for the alternative}
#'      \item{sd}{a \eqn{2\times J} matrix of standard deviations, the first row for the
#'                null setting, the second for the alternative}
#'     }
#'   }
#'   \item{prevalences}{
#'     A list of two elements with prevalence vectors used in the paper; the lengths of these
#'     vectors implicitly define the number of groups.
#'     \describe{
#'      \item{table1}{a vector of equal prevalences for six groups used in table 1}
#'      \item{table2}{a vector of prevalences used in table 2 of the paper}
#'     }
#'   }
#' }
#' @name LLL.SETTINGS
#' @docType data
#' @references Adaptive Choice of Patient Subgroup for Comparing Two Treatments
#' by Tze Leung Lai and Philip W. Lavori and Olivia Yueh-Wen Liao. Contemporary Clinical Trials,
#' Vol. 39, No. 2, pp 191-200 (2014). \url{http://www.sciencedirect.com/science/article/pii/S1551714414001311}
#' @keywords data
#'
NULL

