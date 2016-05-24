#' Repeated measures for two groups of subjects drawn from mmcar model with 2 nuisance covariates
#' 
#' A simulation dataset containing repeated subject measures for 2 treatment groups,
#' (control = 0, treatment = 1), constructed from an 'mmcar' model with correlation
#' between adjacent sessions equal to 0.25.  Subject effects were randomly drawn from 10 clusters
#' with weights/probabilities drawn from a Dirichlet distribution.  Cluster location
#' values were generated from a Gaussian base distribution.
#' Two nuisance demographic variables, age and income, are included.
#' 
#' \itemize{
#'   \item subject. subject identifier \code{(1,2,...,264}
#'   \item trt. treatment group identifier of length \code{N} (e.g. \code{(0,0,0,...,1,1,1,...)} , either \code{{0,1}} for control and treatment. 
#'   \item time. times in months for each repeated subject measure of length \code{N}.  There are 3 distinct time points.  e.g. \code{(0,3,6,0,3,6,0,0,3,,,,)}
#'   \item n.random. number of random effects per subject.  Set = 3.
#'   \item n.fix_degree. order of fixed effects.  Set = 2, for quadratic, meaning 3 effects (intercept, slope, quadratic) each, for treatment and control groups. 
#'   \item coefs. true fixed effect coefficient values used to generate data.
#'   \item subj.aff. indexes subjects receiving treatment.
#'   \item W.subj.aff.  multiple membership weight matrix that maps the \code{P_aff = 132} affected subjects (in \code{subj.aff}) to any of \code{S = 245} treatment sessions.
#'   \item group.  treatment group membership for each of the \code{S} sessions.
#'   \item Omega. the \code{S x S} CAR adjacency matrix used to model prior dependence among sessions
#'   \item gamma. true session effect values (of length \code{S}) used to generate model response.
#'   \item s. true cluster memberships for each of the \code{P} subjects.
#'   \item b.star. a list object of true cluster location values for each of \code{M = 10} clusters.  Each entry contains the \code{n.random = 3} location values for that cluster.
#'   \item b. a list object true random effect coefficient values for each of \code{P} subjects.  Each entry contains the \code{n.random = 3} effect values for that subject.
#'   \item tau.b. true values for the prior precisions of the base Gaussian distribution for each of \code{n.random = 3} subject effects.
#'   \item tau.e. true value for overall model error.
#'   \item coefs. true coefficient values for the time-based quadratic fixed effects generated from the \code{trt}, and \code{time} inputs, as well as the 2 nuisance covariates.
#'	e.g. X = c(1, time, time^2, trt_1,trt_1*time, trt_1*time^2, age, income).
#'   \item formula.  the additive formula containing the response and nuisance fixed effects.  In this case, \code{y ~ age + income}.
#'   \item data. \code{(N = 792) x 3} \code{data.frame} associated to \code{formula}.  Includes model response, \code{y}, a  \code{N = 792 x 1} numeric vector
#'				capturing repeated measures for \code{P = 264} subjects.  Also contained in \code{data} are two nuisance fixed effects, \code{age} and \code{income}.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name datsimcov
#' @usage datsimcov
#' @format A list object of 19 variables for 792 total observations on 264 subjects
NULL


