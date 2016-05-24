#' @name DatA
#' @aliases DatA
#' @docType data
#' @title Database of household expenses for frame A
#' 
#' @description This dataset contains some variables regarding household expenses for a sample of 105 households selected from a list of landline phones (let say, frame A) in a particular city in a specific month.
#' @usage DatA
#' @details The sample, of size \eqn{n_A = 105}, has been drawn from a population of \eqn{N_A = 1735} households with landline phone according to a stratified random sampling. Population units were divided in 6 different strata. 
#'  Population sizes of these strata are \eqn{N_A^h = (727, 375, 113, 186, 115, 219)}. \eqn{N_{ab} = 601} of the households composing the population have, also, mobile phone. On the other hand, frame totals for auxiliary variables in this frame are \eqn{X_{Income}^A = 4300260} and \eqn{X_{Taxes}^A = 215577}.
#' @format
#' \describe{
#'    \item{Domain}{A string indicating the domain each household belongs to. Possible values are "a" if household belongs to domain a or "ab" if household belongs to overlap domain.}
#'    \item{Feed}{Feeding expenses (in euros) at the househould}
#'    \item{Clo}{Clothing expenses (in euros) at the household}
#'    \item{Lei}{Leisure expenses (in euros) at the household}
#'    \item{Inc}{Household income (in euros). Values for this variable are only available for households included in frame A. For households included in domain b, value of this variable is set to 0.}
#'    \item{Tax}{Household municipal taxes (in euros) paid. Values for this variable are only available for households included in frame A. For households included in domain b, value of this variable is set to 0.}
#'    \item{M2}{Square meters of the house. Values for this variable are only available for households included in frame B. For households included in domain a, value of this variable is set to 0.}
#'    \item{Size}{Household size. Values for this variable are only available for households included in frame B. For households included in domain a, value of this variable is set to 0.}
#'    \item{ProbA}{First order inclusion probability in frame A. This probability is 0 for households included in domain b.}
#'    \item{ProbB}{First order inclusion probability in frame B. This probability is 0 for households included in domain a.}
#'    \item{Stratum}{A numeric value indicating the stratum each household belongs to.}
#' }
#' @seealso \code{\link{PiklA}}
#' @examples
#' data(DatA)
#' attach(DatA)
#' #Let perform a brief descriptive analysis for the three main variables
#' param <- data.frame(Feed, Clo, Lei)
#' summary (param)
#' hist (Feed)
#' hist (Clo)
#' hist (Lei)
#' @keywords datasets
NULL