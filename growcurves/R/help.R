#' Bayesian Semi and Nonparametric Growth Curve Models with Employment 
#' of Multiple Membership Random Effects for Longitudinal Data
#'
#' \tabular{ll}{
#' Package: \tab growcurves\cr
#' Type: \tab Package\cr
#' Version: \tab 0.2.4.0\cr
#' Date: \tab 2015-11-13\cr
#' License: \tab GPL (>= 3) \cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' Specifies a non-parametric prior for subject random effects and adds additional sets of dose or exposure random effects
#' that are linked to subjects through a multiple membership construction for application to repeated measures data on subjects.
#' One class of models employs a set of by-subject random effect parameters under a Dirichlet process (DP) prior
#' with the addition of one or more sets of multiple membership (MM) effects that map by-dose effects to subject
#' for data characterized by repeated measures on subject.  There may be specified \code{q} random effect parameters per
#' subject, possibly equal to the number of measurement waves, \code{T}, as the DP prior borrows estimation strength
#' across subjects.  Another class of models employs a single set of by-subject effects under a dependent DP (DDP) prior for a 
#' collection of random distributions that are indexed by MM dose or sequence.  Each set of subject random effects 
#' under the DDP formulation are now indexed by a full set of MM sequences (for all subjects).
#'
#' CORE "ENGINE" FUNCTIONS
#'
#' \code{\link{dpgrow}} performs Bayesian mixed effects estimation for data characterized by repeated subject measures
#' (typically in time).  The user inputs a \code{subject} identifier vector, a vector of \code{time} measurements,
#' and a \code{trt} vector for treatment/group assignments.  Fixed and random effects are then automatically
#' generated and both subject and treatment level growth curves are constructed.
#'
#' \code{\link{dpgrowmm}} is very similar to \code{dpgrow}, but it adds an additional set of exposure random effects
#' which aren't grouped by subject that may be used to inject treatment dosage or attendance patterns that
#' is mapped back to clients via a multiple membership weight matrix.  The option \code{multi = TRUE} specifies
#' each exposure random effect as polynomial in time as is done for the by-subject effects.
#'
#' \code{\link{dpgrowmult}} is very similar to \code{dpgrowmm}, but it allows more than one set of multiple
#' membership effects. Each multiple membership effects term (block) may apply to any sub-set of subjects
#' through specification of the weight matrix and identification of affected subjects for that term.
#'
#' \code{\link{ddpgrow}} generalizes \code{dpgrowmm} and \code{dpgrowmult} by indexing the 
#' subject random effects with a set of exposures linked to subjects from an MM weight matrix.  This model 
#' brings the MM term inside the DP by specifying a set of dose-based (MM) random effects for each subject.   
#' The model formally employs a dependent DP (DDP) prior on a set of subject effects with the a single unknown
#' prior distribution now replaced by a collection of unknown distributions indexed by dosage pattern.  For example,
#' a dosage or exposure may be characterized by a sequence of cognitive behavior therapy sessions attended.
#' 
#' CORE "ACCESSOR" (PLOT) FUNCTIONS
#'
#' \code{\link{growplot}} uses model outputs from \code{dpgrow}, \code{dpgrowmm}, \code{dpgrowmult} and \code{ddpgrow}
#' to provide by-subject growth curves in two forms: 1. Growth curves aggregated under specified groupings;
#' 2. Indvidual growth curves plotted along with data for selected (or random subset of) subjects.
#'
#' \code{\link{trtplot}} uses model outputs from \code{dpgrow}, \code{dpgrowmm}, \code{dpgrowmult} and \code{ddpgrow}
#' to plot a distribution for fixed mean effects difference between two selected treatments 
#' across a range of chosen models for a one or more chosen time points.
#' Outputs include a set of boxplots for each time point that span 95% credible interval.
#'
#' \code{\link{effectsplot}} uses model outputs from \code{dpgrow}, \code{dpgrowmm} and \code{dpgrowmult}
#' to overlay plots of multiple membership effects for a given term under use of different prior 
#' formulations and/or from distinctly formulated models (e.g. with varying numbers of 
#' multiple membership terms).
#'
#' \code{\link{ddpEffectsplot}} uses model outputs from \code{ddpgrow}
#' to produce by subject and by clusters of subjects summaries for the \code{q x (S+1)} multivariate random
#' effects for each subject, where \code{S} denotes the number of unique dosages across all subjects,
#' and \code{q} denotes the polynomial order for each of the \code{S+1} effects.   There is also
#' a \code{q x 1} set of subject intercept effects.  This function is analogous to \code{effectsplot},
#' only each client now has its own set of MM random effects.  
#'
#' SIMULATED DATA SETS
#'
#' There are 3 simulated data sets available in order to allow exploration of the engine 
#' and associated accessor functions.
#'
#' \code{\link{datsim}} Simulated dataset with two treatment arms (treatment and control) composed from a 
#' model with a Dirichlet process (DP) prior on the set of client effects and a single MM term under a 
#' \code{"mmcar"} formulation.  Structured to express similar properties as the case example in both 
#' Savitsky and Paddock (2012) references, below.
#'
#' \code{\link{datsimcov}} Of similar structre to \code{simdat}, only the data generating model now
#' additionally employs 2 nuisance fixed effects.
#'
#' \code{\link{datsimmult}} Simulated data under 2 treatment arms generated from a model with now
#' 2 multiple membership terms.   The terms are generated under \code{c("mmi","mmcar")} prior formulations.  
#'
#' BENCHMARK DATA SETS
#'
#' \code{\link{datbrghtterms}} Data derived from BRIGHT study reviewed in reference and includes
#' BDI-II depressive symptom scores for client experimental units.  Associated data objects
#' are included to facilitate runs under engine functions \code{dpgrow}, \code{dpgrowmm} and \code{dpgrowmult}.
#'
#' \code{\link{datbrghtmodterms}} Data derived from BRIGHT study reviewed in reference and includes
#' BDI-II depressive symptom scores for client experimental units.  CBT treatment sessions are collected into
#' higher level modules.  All objects are then specified by module, rather than session.  These data
#' may be used with any engine function, though were created to facilitate use of \code{ddpgrow}.  Data
#' objects are included, therefore, to enable employment of \code{ddpgrow}, as well as the other engine
#' functions.
#'
#' \code{\link{dateduc}} Tests for students tracked for 5 years from grades 1 - 5 for a single school in 
#' a large urban school district.   Associated data objects are included to facilitate runs under engine
#' functions \code{dpgrow}, \code{dpgrowmm}, \code{ddpgrow}.
#'
#' @examples 
#' \dontrun{
#' ## extract simulated dataset
#' library(growcurves)
#' data(datsim)
#' ## attach(datsim)
#' ## run dpgrow mixed effects model, returning object of class "dpgrow"
#' shape.dp	= 4
#' res		= dpgrow(y = datsim$y, subject = datsim$subject, 
#'				trt = datsim$trt, time = datsim$time,
#'				n.random = 3, n.fix_degree = 2, 
#'				n.iter = 10000, n.burn = 2000, 
#'				n.thin = 10, shape.dp = shape.dp, 
#'				option = "dp")
#' ## Each plot is a "ggplot2" object saved in 
#' ## a list to plot.results
#' plot.results	= plot(res) ## includes subject and 
#' ##                    treatment growth curves
#' ## Extract credible intervals (2.5%, mean, 97.5%).
#' ## Includes fit statistics:  Dbar, DIC, pD, lpml.  
#' ## Note: DIC is the DIC3 of Celeaux et. al. (2006) 
#' ## for option = "dp".  Finally, the constructed fixed
#' ## and random effects matrices, X and Z, are returned 
#' ## with growth curve covariates appended 
#' ## to user submitted nuisance covariates. 
#' summary.results = summary(res)
#' ## View the summary results in the console
#' print(summary.results)
#' ## Collect posterior sampled values over 
#' ## the (n.iter - n.burn) retained iterations 
#' ## for each sampled parameter.  
#' samples.posterior	= samples(res)
#' ## model residuals (y - fit)
#' residuals		= resid(res) 
#' ## Model with DP on clients effects, but 
#' ## now INCLUDE session random effects
#' ## in a multiple membership construction 
#' ## communicated with the N x S matrix, W.subj.aff.
#' ## Returns object, res.mm, of class "dpgrowmm".
#' shape.dp	= 4
#' strength.mm	= 0.1
#' res.mm	= dpgrowmm(y = datsim$y, subject = datsim$subject, 
#'                      trt = datsim$trt, time = datsim$time, 
#'			n.random = 3, 
#'			Omega = datsim$Omega, group = datsim$group, 
#'			subj.aff = datsim$subj.aff,
#'			W.subj.aff = datsim$W.subj.aff, 
#'			n.iter = 10000, n.burn = 2000, n.thin = 10,
#'			strength.mm = strength.mm, 
#'			shape.dp = shape.dp, 
#'			option = "mmcar")
#' plot.results		= plot(res.mm)
#' }
#'
#' @name growcurves-package
#' @aliases growcurves package-growcurves
#' @docType package
#' @author Terrance Savitsky \email{tds151@@gmail.com} Susan Paddock \email{paddock@@rand.org}
#' @references 
#'	S. M. Paddock and T. D. Savitsky (2013) Bayesian Hierarchical Semiparametric Modeling of Longitudinal
#'             Post-treatment Outcomes from Open-enrollment Therapy Groups., 
#'             JRSS Series A (Statistics in Society), 2013, Volume 176, Part 2, pp. 797 - 808.
#' @references
#'	T. D. Savitsky and S. M. Paddock (2013) Bayesian Non-Parametric Hierarchical Modeling for Multiple
#'              Membership Data, Annals of Applied Statistics, Volume 7, Number 2, pp. 1074 - 1094.
#' @references
#' 	T. D. Savitsky and S. M. Paddock (2014) {B}ayesian Semi- and Non-Parametric Models for Longitudinal Data with Multiple Membership Effects in {R}, Journal of Statistical Software,
#'      Volume 57, Number 3, pp. 1 -- 35, http://www.jstatsoft.org/v57/i03/
#' @import reshape2 scales ggplot2 testthat Rcpp RcppArmadillo Formula
#' @useDynLib growcurves
#' @keywords package
NULL