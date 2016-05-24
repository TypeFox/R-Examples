#' @name DatB
#' @aliases DatB
#' @docType data
#' @title Database of household expenses for frame B
#' 
#' @description This dataset contains some variables regarding household expenses for a sample of 135 households selected from a list of mobile phones (let say, frame B) in a particular city in a specific month.
#' @usage DatB
#' @details The sample, of size \eqn{n_B = 135}, has been drawn from a population of \eqn{N_B = 1191} households with mobile phone according to a simple random sampling without replacement design. 
#'  \eqn{N_{ab} = 601} of these households have, also, landline phone. On the other hand, frame totals for auxiliary variables in this frame are \eqn{X_{Metres2}^B = 176553} and \eqn{X_{Size}^B = 3529}
#' @format
#' \describe{
#'    \item{Domain}{A string indicating the domain each household belongs to. Possible values are "b" if household belongs to domain b or "ba" if household belongs to overlap domain.}
#'    \item{Feed}{Feeding expenses (in euros) at the househould}
#'    \item{Clo}{Clothing expenses (in euros) at the household}
#'    \item{Lei}{Leisure expenses (in euros) at the household}
#'    \item{Inc}{Household income (in euros). Values for this variable are only available for households included in frame A. For households included in domain b, value of this variable is set to 0.}
#'    \item{Tax}{Household municipal taxes (in euros) paid. Values for this variable are only available for households included in frame A. For households included in domain b, value of this variable is set to 0.}
#'    \item{M2}{Square meters of the house. Values for this variable are only available for households included in frame B. For households included in domain a, value of this variable is set to 0.}
#'    \item{Size}{Household size. Values for this variable are only available for households included in frame B. For households included in domain a, value of this variable is set to 0.}
#'    \item{ProbA}{First order inclusion probability in frame A. This probability is 0 for households included in domain b.}
#'    \item{ProbB}{First order inclusion probability in frame B. This probability is 0 for households included in domain a.}
#' }
#' @seealso \code{\link{PiklB}}
#' @examples
#' data(DatB)
#' attach(DatB)
#' #Let perform a brief descriptive analysis for the three main variables
#' param <- data.frame(Feed, Clo, Lei)
#' summary (param)
#' hist (Feed)
#' hist (Clo)
#' hist (Lei)
#' @keywords datasets
NULL