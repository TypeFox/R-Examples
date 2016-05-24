#' @name DatPopM
#' @aliases DatPopM
#' @docType data
#' @title Database of auxiliary information for the whole population of students
#' 
#' @description This dataset contains population information about the auxiliary variables of the population of students 
#' @usage DatPopM
#' @details The population size is \eqn{N = 10000}.
#' @format
#' \describe{
#'    \item{Ses}{An ordinal factor with three categories (low, middle and high) indicating the socio-economical status of the student.}
#'    \item{Read}{A number indicating the mark of the student in a reading test.}
#'    \item{Write}{A number indicating the mark of the student in a writing test.}
#'    \item{Domain}{A string indicating the domain each student belongs to. Possible values are "a" if student belongs to domain a, "b" if student belongs to domain b or "ab" if student belongs to overlap domain.}
#' }
#' @seealso \code{\link{DatMA}} \code{\link{DatMB}}
#' @examples
#' data(DatPopM)
#' attach(DatPopM)
#' #Let perform a brief descriptive analysis for the three auxiliary variables
#' summary (Ses)
#' summary(Read)
#' summary(Write)
#' @keywords datasets
NULL