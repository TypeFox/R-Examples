#' Bayesian Non-Parametric Models for Estimating a Set of Denoised, Latent Functions 
#' From an Observed Collection of Domain-Indexed Time-Series
#'
#' \tabular{ll}{
#' Package: \tab growfunctions\cr
#' Type: \tab Package\cr
#' Version: \tab 0.12\cr
#' Date: \tab 2015-11-13\cr
#' License: \tab GPL (>= 3) \cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' Parameterizes a model for a collection of noisy time-series indexed by 
#' domain or observation unit as an additive function process plus 
#' noise process.  The latent functions are modeled under both
#' Gaussian process (GP) and intrinsic Gaussian Markov random (iGMRF) 
#' field priors.  Dependence is estimated among the set of functions by 
#' selecting a prior, G, on the covariance or precision
#' parameters of the GP and iGMRF, respectively, where G receives a 
#' Dirichlet process (DP) prior. The resulting marginal prior on the 
#' functions is a nonparametric scale mixture over the 
#' covariance or precision parameters, with G the unknown mixing measure.  
#' Draws from a DP are almost surely discrete, allowing for ties 
#' (interpreted as clusters) among the covariance or
#' precision parameters indexed by domain or observation unit.  
#' Functions are included that permit additional inference on the 
#' clustering properties over the domains. The mixture models 
#' allow specification of multiple additive latent functions for 
#' each of the GP and iGMRF DP mixture formulations. Each 
#' function (under a GP prior) may be specified with a covariance 
#' function selected from a set of available options that allow 
#' for combinations of trend and seasonal components.  
#' The same is true for precision matrix 
#' constructions under the iGMRF DP mixture.
#'
#' ESTIMATION FUNCTIONS
#'
#' \code{\link{gpdpgrow}} performs Bayesian nonparametric estimation of a set of dependent, denoised
#' latent functions under a DP mixture of GP's in an unsupervised fashion based on user input of
#' an \code{N x T} data matrix, \code{y}, where \code{N} denotes the number of domains or units
#' and \code{T} denotes the number of time points.  The DP prior estimates the dependent structure 
#' among the covariance parameters generating each \code{T x 1} function across the \code{N} domains.
#' The user may specify multiple latent functions in an additive 
#' regression formulation, where
#' each covariance kernel may be selected from the squared exponential (\code{"se"}), the
#' rational quadratic (\code{"rq"}), or a quasi-periodic (product of a period and 
#' squared exponential
#' to allow the periodic kernel to evolve) (\code{"sn"}).
#'
#' \code{\link{gmrfdpgrow}} also inputs an \code{N x T} data matrix, \code{y}, but replaces the 
#' GP prior used in \code{gpdpgrow()} with an iGMRF prior.  The DP mixture is over the precision 
#' parameter, \code{kappa}, that multiplies a fixed matrix, \code{Q}, which specifies 
#' the length-scale
#' of dependence.  The user may specify multiple functions in an additive formulation, 
#' each with its own precision matrix.  The precision matrix is specified 
#' via two options.  Input \code{q_type} 
#' indicates whether the precision matrix is a \emph{trend}  (\code{"tr"}) 
#' precision matrix or a 
#' \emph{seasonal} (\code{"sn"}) precision matrix. Integer innput \code{q_order} 
#' controls the length scale of the estimated functions and specifies the 
#' difference order of the precision matrix in the case of \code{q_type = "tr"}, 
#' or the periodicity, in the case \code{q_type = "sn"}.  Typically, \code{q_order} is set to
#' \code{1} or \code{2} when \code{q_type = "tr"}, though other values are possible.
#' 
#' \code{\link{MSPE}} Inputs a \code{gpdpgrow()} or \code{gmrfdpgrow()} object estimated where
#' some data values deliberately set to missing (\code{NA}) and produces an out-of-sample 
#' mean square prediction error (MSPE), and a normalized MSPE (normalized by the variance of
#' the missing test set). Both \code{gpdpgrow()} and \code{gmrfdpgrow()} returned objects also
#' include a leave-one-out log-pseudo marginal likelihood \code{LPML} fit statistic estimated
#' on the training set (so that no data are required to be left out).  One would expect the
#' \code{nMSPE} statistic to provide a greater penalty for model complexity since it is computed
#' on data not used for estimation.
#' 
#' \code{\link{predict_functions}} Uses the model-estimated GP covariance parameters from 
#' \code{gpdpgrow()} or iGMRF precision parameters from \code{gmrfdpgrow()}
#' to predict functions at \emph{future} time points beyond the range of the data 
#' from their posterior predictive distributions. Both \code{gpdpgrow()} and \code{gmrfdpgrow()}
#' will predict missing function values in the range of the data.  
#' 
#' PLOT FUNCTIONS
#'
#' \code{\link{cluster_plot}} inputs a returned object from either of \code{gpdpgrow()} or
#' \code{gmrfdpgrow()} and produces two plots.  The first plot creates panels for the clusters
#' and aggregates line plots of posterior mean estimates for member denoised functions.  
#' An optional 
#' smoother may be drawn through the set of functions in each panel to differentiate the
#' patterns expressed across clusters.  A second plot renders the estimated posterior mean
#' values (with an option for credible intervals) for a single or group of randomly-selected
#' latent functions against the actual data values, \code{y}.
#' 
#' \code{\link{informative_plot}} inputs a list of returned objects, all from either
#' of \code{gpdpgrow()} or \code{gmrfdpgrow()} (where the model type, GP or 
#' iGMRF, is communicated with input \code{model}) that compares 
#' credible intervals for covariance or precision
#' parameters where the data are drawn from an informative 
#' sampling design (rather than as \emph{iid}), so that the distribution 
#' for the population is not the same as that for
#' the sample.  One model conducts estimation in a fashion that ignores the informative
#' design and the other incorporates sampling weights to account for the informativeness.
#' Comparing the resulting estimated credible intervals provides a means to assess the 
#' sensitivity of estimated parameters to the sampling design.  The set of objects are
#' labeled with \code{objects_labels} with allowable inputs \code{c("ignore","weight","iid")}.
#' One of the objects must have label \code{"ignore"}, and the other must have label,
#' \code{weight}. An additional object may also be input that is estimated from an iid sample
#' drawn from the same population as the informative sample used for estimation under
#' objects \code{"ignore"} and \code{"weight"}.  Both informative and iid samples are generated
#' from the synthetic data function, \code{gen_informative_sample()}.
#'
#' \code{\link{fit_compare}} inputs a list of returned objects, each from either
#' \code{gpdpgrow()} or \code{gmrfdpgrow()}, and plots the posterior mean estimate
#' for a randomly-selected latent function as compared to the actual data in a set
#' panels indexed by cluster and object.  Allows comparison of fit performance of
#' functions to data across varied model specifications.
#'
#' \code{\link{predict_plot}} uses a returned object from \code{predict_functions()}
#' to plot both estimated and predicted function values (with the option for 
#' credible intervals) where the prediction interval is outside the range of data. 
#' The plot function, \code{\link{cluster_plot}}, should be used where 
#' the predicted values are in the range
#' of the data (and may, hence, be treated as missing values.  The estimated functions are
#' plotted in \code{cluster_plot} for the missing, as well as observed, data points).
#'
#' DATA SETS (and functions to generate synthetic data sets)
#'
#' \code{\link{cps}} Data derived from the Current Population Survey (CPS) administered 
#' to households in
#' local areas within each of 51 states.  These data capture a rectangular matrix of 
#' monthly-indexed
#' all-other employment direct estimates computed from the set of household responses.  
#' The data capture
#' 156 month observations (from 2000 - 2013) for each of the 51 states.  Two objects are included:
#' \code{y_raw} contains the 51 x 156 matrix of all-other employment counts.  
#' \code{y} contains the 51 x 156 matrix of all-other 
#' employment counts are standardization to (0,1), by state, to all
#' them to be modeled together.
#'
#' \code{\link{gen_informative_sample}} Generates an \code{N x T} population data matrix, 
#' \code{y}, and an associated \code{n x T} sample data matrix, \code{y_obs}, where the 
#' sample is drawn using a 1 or 2-stage informative process.  The 1-stage sample uses
#' unequal probability stratified sampling, while the 2-stage process samples the first stage
#' in blocks, while the second stage samples with unequal probability from strata for selected
#' blocks.  Both block and strata memberships for population units are generated based on the
#' variance of their time-series, \code{y}.  The resulting sample is informative because the
#' block and cluster memberships and selection probabilities are determined based on \code{y}.
#'
#'
#' @examples 
#' {
#' library(growfunctions)
#' 
#' ## load the monthly employment count 
#' ## data for a collection of 
#' ## U.S. states from the Current 
#' ## Population Survey (cps)
#' data(cps)
#' ## subselect the columns of N x T, y, 
#' ## associated with 
#' ## the years 2011 - 2013
#' ## to examine the state level 
#' ## employment levels 
#' ## during the "great recession"
#' y_short <- cps$y[,(cps$yr_label %in% 
#'                   c(2011:2013))]
#'
#' ## run DP mixture of GP's to estimate 
#' ## posterior distributions 
#' ## for model parameters
#' ## uses default setting of a single 
#' ## "rational quadratic" 
#' ## covariance formula
#' ## A short number of iterations is used 
#' ## for illustration
#' ## Run for 500 iterations with half 
#' ## as burn-in to 
#' ## get a more useful result
#' res_gp     <- gpdpgrow(
#'                    y = y_short, 
#'                    n.iter = 3, 
#'                    n.burn = 1, 
#'                    n.thin = 1, 
#'                    n.tune = 0)  
#' ## 2 plots of estimated functions: 
#' ## 1. faceted by cluster and fit;
#' ## 2.  data for experimental units.
#' ## for a group of randomly-selected 
#' ## functions
#' fit_plots_gp   <- cluster_plot( 
#'   object = res_gp, units_name = "state", 
#'   units_label = cps$st, single_unit = FALSE, 
#'   credible = TRUE )
#'                                    
#' ## Run the DP mixture of iGMRF's 
#' ## to estimate posterior 
#' ## distributions for model parameters
#' ## Under default RW2(kappa) = order 2 trend 
#' ## precision term
#' ## A short number of iterations 
#' ## is used for illustration
#' ## Run for 2000 iterations with 
#' ## half as burn-in to 
#' ## get a more useful result
#' res_gmrf <- gmrfdpgrow(
#'                  y = y_short, 
#'                  n.iter = 11, 
#'                  n.burn = 4, 
#'                  n.thin = 1) 
#'                                      
#' ## 2 plots of estimated functions: 
#' ## 1. faceted by cluster and fit;
#' ## 2.  data for experimental units.
#' ## for a group of 
#' ## randomly-selected functions
#' fit_plots_gmrf <- cluster_plot( 
#'   object = res_gmrf, units_name = "state", 
#'   units_label = cps$st, single_unit = FALSE, 
#'   credible = TRUE )                                    
#'                                      
#' ## visual comparison of fit performance 
#' ## between gpdpgrow() and gmrfdpgrow()
#' ## or any two objects returned from any
#' ## combination of these estimation
#' ## functions
#' objects        <- vector("list",2)
#' objects[[1]]   <- res_gmrf
#' objects[[2]]   <- res_gp
#' label.object   <- c("gmrf_tr2","gp_rq")
#' ## the map data.frame 
#' ## object from fit_plots gp 
#' ## includes a field that 
#' ## identifies cluster assignments
#' ## for each unit (or domain)
#' H        <- fit_plots_gp$map$cluster
#' fit_plot_compare_facet <- 
#' fit_compare( objects = objects, 
#'  H = H, label.object = label.object,
#'  y.axis.label = "normalized y",
#'  units_name = "state", units_label = cps$st)                                
#' }
#'
#' @name growfunctions-package
#' @aliases growfunctions package-growfunctions
#' @docType package
#' @author Terrance Savitsky \email{tds151@@gmail.com} 
#' @references
#'     T. D. Savitsky (2014) Bayesian Non-parametric Functional Mixture
#'     Estimation for Time-indexed data. submitted to: Survey Methodology.
#' @references
#'	T. D. Savitsky, D. Toth and M. Sverchkov (2014) Bayesian Estimation 
#'	Under Informative Sampling, 
#'     Submitted to: Electronic Journal of Statistics.
#' @references
#'     T. D. Savitsky (2014) Bayesian Non-Parametric Mixture Estimation 
#'     for Time-Indexed Functional
#'     Data for \code{R}. To Appear in: Journal of Statistical Software.     
#' @import reshape2 scales ggplot2 Rcpp RcppArmadillo Matrix spam mvtnorm
#' @useDynLib growfunctions
#' @keywords package
NULL