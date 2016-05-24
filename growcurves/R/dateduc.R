#' Student test scores and associated teachers for a single school in a large urban school district
#' 
#' Response captures vertically linked mathematics and reading scores on a norm-referenced 
#' standardized test administered during the spring of the years 1998 to 2002 obtained from
#' a large urban school district.  Included student are in grade 1 during the 1997 - 1998 school
#' year and followed successively until grade 5 for the 2001-2002 school year.  These data focus
#' on a single school of 227 students and 34 teachers.   Teacher links to students did not vary
#' by subject.  The data are configured to support model runs using engine functions \code{dpgrowmm} and \code{ddpgrow}.
#' 
#' \itemize{
#'   \item y. Numeric vector of \code{N = 562} student-year test scores.			  
#'   \item subject. Numeric vector, \code{1, ..., n = 227}, student identifiers
#'   \item trt. A numeric vector of \code{N} 0's to indicate there are not separate treatment arms.  	 
#'   \item time. An {N x 1} numeric vector of associated meaurement times in \code{1-5}.
#'   \item subj.aff. Same as \code{subject} as all students receive the "school" treatement.
#'   \item W.subj.aff. An \code{n = 227 x S = 34} matrix object that links the \code{n = 227} students to the \code{S = 34} teachers.  The rows sum to 1.
#'			This object is used in engine function \code{dpgrowmm}.
#'   \item dosemat. An \code{n = 227 x S = 34} numeric matrix where the first column is 1's for the intercept and the first teacher is left out for identifiability.
#'			This matrix is the same as \code{W.subj.aff} except that the first column is replaced with an intercept.  This object is use for engine function \code{ddpgrow}.
#'   \item tchr.aff. A numeric vector with values \code{1,...,T=34} encoding the id's for the participating teachers.
#'   \item labt. A character input providing the label \code{"school"} for the treatment delivered to students
#' }
#' 
#' @references
#' 	J.R. Lockwood, Daniel F. McCaffrey, Louis T. Mariano and Claude Setodji (2007) Bayesian Methods for Scalable Multivariate
#'	Value-Added Assessment, Journal of Educational and Behavioral Statistics, 32(2), 125 - 150.
#' @references
#'	S. M. Paddock and T. D. Savitsky (2012) Bayesian Hierarchical Semiparametric Modeling of Longitudinal Post-treatment Outcomes from Open-enrollment Therapy Groups, invited re-submission to: 
#'      JRSS Series A (Statistics in Society).
#' @docType data
#' @keywords datasets
#' @name dateduc
#' @usage dateduc
#' @format A list object for 562 total observations on 227 subjects
NULL