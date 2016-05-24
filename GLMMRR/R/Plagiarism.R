#' An Experimental Survey Measuring Plagiarism Using the Crosswise Model
#'
#' A dataset containing the reponses to sensitive questions about plagiarism
#' and other attributes of 812 students. The crosswise model (CM) and direct questioning (DQ)
#' were utilized to gather the data. Each row holds the response to one question for one student.
#' The variables are as follows:
#'
#' \itemize{
#'   \item id. identification code of the student
#'   \item question. which question was asked (1 and 3: Partial Plagiarism, 2 and 4: Severe Plagiarism)
#'   \item gender. gender of the student (0: male, 1: female)
#'   \item age. age in years
#'   \item nationality. nationality of the student (0: German or Swiss, 1: other)
#'   \item no_papers. number of papers
#'   \item uni. location of data collection (1: ETH Zurich, 2: LMU Munich, 3: University Leipzig)
#'   \item course. course in which the data was collected
#'   \item Aspired_Degree. aspired degree of the student
#'   \item Semester. semesters enrolled
#'   \item ur_none. used resources: none
#'   \item ur_books. used resources: books
#'   \item ur_art. used resources: articles
#'   \item ur_int. used resources: internet
#'   \item ur_fsp. used resources: fellow students' papers
#'   \item ur_other. used resources: other
#'   \item preading. proofreading
#'   \item gradesf. satisfaction with grades
#'   \item pp. Plagiarism indicator (0: Severe Plagiarism, 1: Partial Plagiarism)
#'   \item RR. Randomized Response indicator (0: DQ, 1: Crosswise)
#'   \item RRp1. Randomized Response parameter p1
#'   \item RRp2. Randomized Response parameter p2
#'   \item RRmodel. Randomized Response Model
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Plagiarism
#' @author Ben Jann and Laurcence Brandenberger
#' @references \url{http://dx.doi.org/10.7892/boris.51190}
#' @usage data(Plagiarism)
#' @format A data frame with 812 rows and 24 variables
NULL
