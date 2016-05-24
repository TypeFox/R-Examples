#' This packages provides functions to estimate and visualize propensity
#' score analyses including matching for non-binary treatments.
#' 
#' @name TriMatch-package
#' @aliases TriMatch
#' @docType package
#' @title Propensity Score Analysis for Non-Binary Treatments
#' @author Jason Bryer \email{jason@@bryer.org}
#' @keywords propensity score analysis psa matching
#' @seealso \code{PSAgraphics} \code{multilevelPSA}
#' @import PSAgraphics
#' @importFrom psych describe
#' @importFrom psych describeBy
#' @import compiler
#' @import stats
#' @import gridExtra
#' @import ggplot2
#' @import ez
#' @import reshape2
#' @importFrom graphics plot
#' @importFrom scales percent
#' @importFrom grid grob editGrob vpPath viewport vpTree grid.layout getGrob gTree
#'             grobWidth grobHeight pushViewport grid.draw upViewport grid.newpage
NA

#' Results from a study examining the effects of tutoring services on course grades.
#' 
#' \itemize{
#' \item \code{treat} Treatment indicator.
#' \item \code{Course} The course id the student was enrolled in.
#' \item \code{Grade} The course grade the student earned (4=A, 3=B, 2=C, 1=D, 0=F or W).
#' \item \code{Gender} Gender of the student.
#' \item \code{Ethnicity} Ethnicity of the student, either White, Black, or Other.
#' \item \code{Military} Is the student an active military student.
#' \item \code{ESL} English second language student.
#' \item \code{EdMother} Education level of the mother (1 = did not finish high school;
#' 2 = high school grad; 3 = some college; 4 = earned associate degree;
#' 5 = earned baccalaureate degree; 6 = Earned Master's degree; 7 = earned doctorate).
#' \item \code{EdFather} Education level of the father (levels same as EdMother).
#' \item \code{Age} Age at the start of the course.
#' \item \code{Employment} Employment level at college enrollment (1 = No; 2 = part-time;
#' 3 = full-time).
#' \item \code{Income} Household income level at college enrollment (1 = <25K; 2 = <35K; 3 = <45K;
#' 4 = <55K; 5 = <70K; 6 = <85K; 7 = <100K; 8 = <120K; 9 = >120K).
#' \item \code{Transfer} Number of transfer credits at the start of the course.
#' \item \code{GPA} GPA as of the start of the course.
#' \item \code{GradeCode} Letter grade.
#' \item \code{Level} Level of the course, either Lower or Upper.
#' \item \code{ID} Randomly assigned student ID.
#' }
#'
#' @name tutoring
#' @docType data
#' @format a data frame with 17 variables.
#' @keywords datasets
NULL

#' Results from the 1987 National Medical Expenditure Study
#' 
#' This file was originally prepared by Anders Corr (corr@@fas.harvard.edu) who
#' reports on December 8, 2007 that the resulting numbers closely
#' match with those reported in the published article. It was later modified by
#' Jason Bryer (jason@@bryer.org) to an R data object to be included in this 
#' package. See \url{http://imai.princeton.edu/research/pscore.html} for more
#' information
#'
#' @name nmes
#' @docType data
#' @format a data frame with 9,708 observations of 12 variables.
#' @keywords datasets
#' @source http://imai.princeton.edu/research/pscore.html
#' @author United States Department of Health and Human Services. 
#'         Agency for Health Care Policy and Research
#' @references National Center For Health Services Research, 1987. National 
#'         Medical Expenditure Survey. Methods II. Questionnaires and data 
#'         collection methods for the household survey and the Survey of 
#'         American Indians and Alaska Natives. National Center for Health 
#'         Services Research and Health Technology Assessment.
#'         
#'         Imai, K., & van Dyk, D.A. (2004). Causal Inference With General 
#'         Treatment Regimes: Generalizing the Propensity Score, Journal of the 
#'         American Statistical Association, 99(467), pp. 854-866.
#'         
#'         Elizabeth Johnson, E., Dominici, F., Griswold, M., & Zeger, S.L. (2003).
#'         Disease cases and their medical costs attributable to smoking: An 
#'         analysis of the national medical expenditure survey. Journal of 
#'         Econometrics, 112.
NULL
