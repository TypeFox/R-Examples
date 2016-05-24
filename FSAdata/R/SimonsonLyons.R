#'Catches in removal events of trout at various locations.
#'
#'Catches of various species in consecutive removal events at various locations.
#'
#'@name SimonsonLyons
#'@docType data
#'@format A data frame of 58 observations on the following 7 variables:
#'\describe{
#' \item{species}{Species of fish.}
#' \item{stream}{Stream name.}
#' \item{first}{Catch on the first removal pass.}
#' \item{second}{Catch on the second removal pass.}
#' \item{third}{Catch on the third removal pass.}
#' \item{fourth}{Catch on the fourth removal pass.}
#' \item{pop.cs}{Population estimate by Carle-Strub method.}
#'}
#'@section Topic(s): \itemize{
#' \item Population size
#' \item Abundance
#' \item Removal
#'}
#'@concept Abundance 'Population Size' Removal
#'@source From Appendix in Simonson, T.D. and J. Lyons.  1995.  Comparison of
#'catch per effort and removal procedures for sampling stream fish assemblages.
#'North American Journal of Fisheries Management, 15:419-427.
#'@keywords datasets
#'@examples
#'data(SimonsonLyons)
#'str(SimonsonLyons)
#'head(SimonsonLyons)
#'
#'## extract data for one species and stream (e.g., 3rd row)
#'SimonsonLyons[3,]
#'
NULL