#' pubdata1
#'
#' @name pubdata1
#'
#' @aliases pd1.count.matrix pd1.similarity pd1.uncat.refs
#'
#' @description Small example dataset with 5 publications that have most of their references categorized into disciplines.
#' The dataset contains the following information: A matrix of counts of referenced disciplines for each publication, a vector of counts of uncategorized references in each publication, and a matrix that contains a measure of similarity between disciplines.
#'
#' @format
#' \describe{
#'  \item{\code{pd1.count.matrix}}{the count of referenced disciplines for each publication}
#'  \item{\code{pd1.uncat.refs}}{the count of referenced disciplines for each publication}
#'  \item{\code{pd1.similarity}}{between disciplines as given in Porter and Rafols, 2009.}
#' }
#' @docType data
#' @usage data("pubdata1")
#' @references
#' Porter, A. and Rafols, I. (2009) Is science becoming more interdisciplinary? Measuring and mapping six research fields over time. Scientometrics, Vol. 81, No. 3 (719-745). DOI:10.1007/s11192-008-2197-2
#' @keywords datasets
NULL

#' pubdata2
#'
#' @name pubdata2
#'
#' @aliases pd2.count.matrix pd2.similarity pd2.uncat.refs
#'
#' @description Small example dataset with 2 publications.
#' The first publication references rather dissimilar disciplines and has uncategorized references.
#' Therefore, the computation of the interval of uncertainty of the Rao-Stirling index requires longer computation time.
#' The dataset contains the following information: A matrix of counts of referenced disciplines for each publication, a vector of counts of uncategorized references in each publication, and a matrix that contains a measure of similarity between disciplines.
#'
#' @format
#' \describe{
#'  \item{\code{pd2.count.matrix}}{the count of referenced disciplines for each publication}
#'  \item{\code{pd2.uncat.refs}}{the count of referenced disciplines for each publication}
#'  \item{\code{pd2.similarity}}{between disciplines as given in Porter and Rafols, 2009.}
#' }
#' @docType data
#' @usage data("pubdata2")
#' @references
#' Porter, A. and Rafols, I. (2009) Is science becoming more interdisciplinary? Measuring and mapping six research fields over time. Scientometrics, Vol. 81, No. 3 (719-745). DOI:10.1007/s11192-008-2197-2
#' @keywords datasets
NULL
