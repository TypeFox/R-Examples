#' Elementary School Teacher Workload Data
#' 
#' Selected variables from a study on elementary school teacher workload
#' in Maricopa County, Arizona.
#' @name teachers
#' @docType data
#' @format data frame with the following 6 variables:
#' \describe{
#'   \item{dist}{school district size; factor with levels \code{large} and
#'     \code{me/sm} (medium/small)}
#'   \item{school}{school identifier}
#'   \item{hrwork}{number of hours required to work at school per week}
#'   \item{size}{class size}
#'   \item{preprmin}{minutes spent per week in school on preparation}
#'   \item{assist}{minutes per week that a teacher's aide works with the
#'      teacher in the classroom}   
#' }
#' @note The study is described in Exercise 16 of Chapter 15. The psu sizes
#' are given in \code{\link{teachmi}}. The large stratum had 245 schools; the 
#' small/medium stratum had 66 schools.
#' @source Data courtesy of Rita Gnap (1995).
#' @references Gnap, R. (1995). Teacher load in Arizona elementary school 
#' districts in Maricopa County. Ph.D. diss., Arizona State University 
#' 
#' Lohr (1999). Sampling: Design and Analysis, Duxbury, p. TODO and
#'   446.
#' @export
NULL
