#' @name DatMB
#' @aliases DatMB
#' @docType data
#' @title Database of students' program choice for frame B
#' 
#' @description This dataset contains some variables regarding the program choice for a sample of 232 students included in the sampling frame B.
#' @usage DatMB
#' @details The sample, of size \eqn{n_B = 232}, has been drawn from a population of \eqn{N_B = 6500} students according to a simple random sampling design.
#'  \eqn{N_{ab} = 2000} of the students composing the population belongs also to frame A.
#' @format
#' \describe{
#'    \item{Id_Pop}{An integer from 1 to \eqn{N}, with \eqn{N} the number of students in the whole population, identifying the student within the population.}
#'    \item{Id_Frame}{An integer from 1 to \eqn{N_B}, with \eqn{N_B} the number of students in the frame, identifying the student within the frame.}
#'    \item{Prog}{A factor with three categories (academic, general and vocation) indicating the program choice of the student.}
#'    \item{Ses}{An ordinal factor with three categories (low, middle and high) indicating the socio-economical status of the student.}
#'    \item{Read}{A number indicating the mark of the student in a reading test.}
#'    \item{Write}{A number indicating the mark of the student in a writing test.}
#'    \item{Sch_Size}{A number indicating the size of the school the students belongs to.}
#'    \item{Domain}{A string indicating the domain each student belongs to. Possible values are "b" if student belongs to domain b or "ba" if student belongs to overlap domain.}
#'    \item{ProbA}{First order inclusion probability in frame A. This probability is 0 for students included in domain b.}
#'    \item{ProbB}{First order inclusion probability in frame B.}
#' }
#' @seealso \code{\link{DatPopM}}
#' @examples
#' data(DatMB)
#' attach(DatMB)
#' #Let perform a brief descriptive analysis for the main variable
#' summary (Prog)
#' #And let do the same for the numerical auxiliary variables Read and Write
#' summary(Read)
#' summary(Write)
#' @keywords datasets
NULL