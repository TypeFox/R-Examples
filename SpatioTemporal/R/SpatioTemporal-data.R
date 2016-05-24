### Raw data ###
##' The raw data that was used to create the \code{\link{mesa.model}} structures.
##' \cr
##' The data structure contains raw data from the \strong{MESA Air} project. The
##' example below describes how to create the \code{\link{mesa.model}} structure
##' from raw data.
##' 
##' @title Data used in the examples
##' @name mesa.data.raw
##' @docType data
##' @format The structure contains observations, temporal trends, locations,
##'   geographic covariates, and spatio-temporal covariates. The data is
##'   stored as a list with elements:
##'   \describe{
##'     \item{X}{A data.frame containing names, locations, and (geographic)
##'              covariates for all the (observation) locations.}
##'     \item{obs}{A time-by-location matrix for the observed data, missing data
##'                marked as \code{NA}}
##'     \item{lax.conc.1500}{A time-by-location matrix of a spatio-temporal
##'                          covariate based on output from Caline3QHC.}
##'   }
##' 
##' @source Contains monitoring data from the \strong{MESA Air} project, see
##'   Cohen et.al. (2009) for details.
##' @keywords datasets
##'
##' @example Rd_examples/Ex_mesa_data_raw.R
##' 
##' @references M. A. Cohen, S. D. Adar, R. W. Allen, E. Avol, C. L. Curl, T.
##'   Gould, D. Hardie, A. Ho, P. Kinney, T. V. Larson, P. D. Sampson, L.
##'   Sheppard, K. D. Stukovsky, S. S. Swan, L. S. Liu, J. D. Kaufman. (2009)
##'   Approach to Estimating Participant Pollutant Exposures in the Multi-Ethnic
##'   Study of Atherosclerosis and Air Pollution (MESA Air). Environmental Science
##'   & Technology: 43(13), 4687-4693.
##' 
##' @seealso \code{\link{createSTdata}} for creation of \code{STdata} objects.
##' @family data matrix
##' @family example data
NULL

### STmodel ###
##' Example of a model structure holding observations, geographic covariates,
##' observation locations, smooth temporal trends, spatio-temporal covariates,
##' and covariance specifications for the model.
##'
##' A \code{STmodel} object consists of a list with, some or all of, the following
##' elements:
##' \describe{
##'   \item{obs}{A data.frame with columns:
##'     \describe{
##'       \item{obs}{The value of each observation.}
##'       \item{date}{The observations time, preferably of class
##'                   \code{\link[base:Date]{Date}}.}
##'       \item{ID}{A \code{character}-class giving observation locations;
##'                 should match elements in \code{locations$ID}.}
##'       \item{idx}{match between \code{obs$ID} and \code{locations$ID} for
##'                  faster computations.}
##'     }
##'     The data.frame is sorted by \code{date} and \code{idx}.
##'   }
##'   \item{locations.list,locations}{Specification of locations and data.frame
##'     with locations for observations (and predictions), see
##'     \code{\link{processLocation}}.
##'   }
##'   \item{D.nu,D.beta}{Distance matrices for the locations in the, possibly
##'     different coordinate systems for beta- and nu-fields. See
##'     \code{\link{processLocation}}.
##'   }
##'   \item{cov.beta,cov.nu}{Covariance structure for beta- and nu-fields, see
##'     \code{\link{updateCovf}}.
##'   }
##'   \item{LUR.list,LUR}{Specification of covariates for the beta-fields and
##'     a list with covariates for each of the beta-fields, see
##'     \code{\link{processLUR}} and \code{\link{createLUR}}.
##'   }
##'   \item{trend,trend.fnc}{The temporal trends with \emph{one of the} columns
##'     being named \code{date}, preferably of class \code{\link[base:Date]{Date}}
##'     providing the time alignment for the temporal trends.}
##'   \item{F}{A matrix contaning  smooth temporal trends for each observation;
##'     elements taken from \code{trend}.}
##' 
##'   \item{ST.list,ST,ST.all}{Spatio-termporal covariates, \code{NULL} if no
##'     covariates. For the observations and all space-time locations respectively,
##'     see \code{\link{processST}} and \code{\link{createST}}.
##'   }
##'   \item{old.trend,fit.trend}{Additional components added if the observations
##'                              have been detrended, see
##'                              \code{\link{detrendSTdata}}.
##'   }
##' }
##' 
##' @title Example of a \code{STmodel} structure
##' @name mesa.model
##' @docType data
##' @format A list with elements, a detailed description of each elements is
##' given in details below
##' 
##' @source Contains monitoring data from the \strong{MESA Air} project, see
##'   Cohen et.al. (2009) and \code{\link{mesa.data.raw}} for details.
##' @keywords datasets
##'
##' @example Rd_examples/Ex_mesa_model.R
##' 
##' @references M. A. Cohen, S. D. Adar, R. W. Allen, E. Avol, C. L. Curl, T.
##'   Gould, D. Hardie, A. Ho, P. Kinney, T. V. Larson, P. D. Sampson, L.
##'   Sheppard, K. D. Stukovsky, S. S. Swan, L. S. Liu, J. D. Kaufman. (2009)
##'   Approach to Estimating Participant Pollutant Exposures in the Multi-Ethnic
##'   Study of Atherosclerosis and Air Pollution (MESA Air). Environmental Science
##'   & Technology: 43(13), 4687-4693.
##' 
##' @seealso \code{\link{createSTmodel}} for creation of \code{STmodel}
##'   objects. \cr \code{\link{createSTdata}} for creation of the originating
##'   \code{STdata} object. 
##' @family example data
NULL

### estimateSTmodel ###
##' Example of a model structure holding parameter estimates for the model in
##' \code{\link{mesa.model}} using \code{\link{estimate.STmodel}}. Estimation
##' results are also provided for models including spatio-temporal covariates.
##' 
##' @title Examples of \code{estimateSTmodel} structure
##' @name est.mesa.model
##' @docType data
##' @format A list with elements, see the return description in
##'   \code{\link{estimate.STmodel}}.
##' 
##' @source Contains parametere estimates for the Spatio-Temporal model applied
##'   to monitoring data from the \strong{MESA Air} project, see
##'   Cohen et.al. (2009) and \code{\link{mesa.data.raw}} for details.
##' @keywords datasets
##'
##' @example Rd_examples/Ex_estimate_STmodel.R
##' 
##' @references M. A. Cohen, S. D. Adar, R. W. Allen, E. Avol, C. L. Curl, T.
##'   Gould, D. Hardie, A. Ho, P. Kinney, T. V. Larson, P. D. Sampson, L.
##'   Sheppard, K. D. Stukovsky, S. S. Swan, L. S. Liu, J. D. Kaufman. (2009)
##'   Approach to Estimating Participant Pollutant Exposures in the Multi-Ethnic
##'   Study of Atherosclerosis and Air Pollution (MESA Air). Environmental Science
##'   & Technology: 43(13), 4687-4693.
##' 
##' @seealso \code{\link{estimate.STmodel}} for parameter estimation. \cr
##'   \code{\link{createSTmodel}} for creation of the originating \code{STmodel}
##'   object.
##' @family example data
NULL

### predictSTmodel ###
##' Example of a predictions for the model in \code{\link{mesa.model}} using
##' \code{\link{predict.STmodel}}. Two sets of predictions are presented,
##' \code{pred.mesa.model} and \code{pred.mesa.model.obs}.
##' 
##' @title Example of a \code{predictSTmodel} structure
##' @name pred.mesa.model
##' @docType data
##' @format A list with elements, see the return description in
##'   \code{\link{predict.STmodel}}.
##' 
##' @source Contains parametere estimates for the Spatio-Temporal model applied
##'   to monitoring data from the \strong{MESA Air} project, see
##'   Cohen et.al. (2009) and \code{\link{mesa.data.raw}} for details.
##' @keywords datasets
##'
##' @example Rd_examples/Ex_predict_STmodel.R
##' 
##' @references M. A. Cohen, S. D. Adar, R. W. Allen, E. Avol, C. L. Curl, T.
##'   Gould, D. Hardie, A. Ho, P. Kinney, T. V. Larson, P. D. Sampson, L.
##'   Sheppard, K. D. Stukovsky, S. S. Swan, L. S. Liu, J. D. Kaufman. (2009)
##'   Approach to Estimating Participant Pollutant Exposures in the Multi-Ethnic
##'   Study of Atherosclerosis and Air Pollution (MESA Air). Environmental Science
##'   & Technology: 43(13), 4687-4693.
##' 
##' @seealso \code{\link{predict.STmodel}} for prediction. \cr
##'   \code{\link{createSTmodel}} for creation of the originating \code{STmodel}
##'   object.
##' @family example data
NULL

### estCVSTmodel & predCVSTmodel ###
##' Example of 10-fold cross-validated for the model in \code{\link{mesa.model}}
##' using \code{\link{estimateCV.STmodel}} and \code{\link{predictCV.STmodel}}.
##' 
##' @title Example of \code{estCVSTmodel} and \code{predCVSTmodel} structures
##' @name est.cv.mesa
##' @aliases pred.cv.mesa
##' @docType data
##' @format A list with elements, see the return description in
##'   \code{\link{estimateCV.STmodel}} and \code{\link{predictCV.STmodel}}.
##' 
##' @source Contains parametere estimates for the Spatio-Temporal model applied
##'   to monitoring data from the \strong{MESA Air} project, see
##'   Cohen et.al. (2009) and \code{\link{mesa.data.raw}} for details.
##' @keywords datasets
##'
##' @example Rd_examples/Ex_estimateCV_STmodel.R
##' 
##' @references M. A. Cohen, S. D. Adar, R. W. Allen, E. Avol, C. L. Curl, T.
##'   Gould, D. Hardie, A. Ho, P. Kinney, T. V. Larson, P. D. Sampson, L.
##'   Sheppard, K. D. Stukovsky, S. S. Swan, L. S. Liu, J. D. Kaufman. (2009)
##'   Approach to Estimating Participant Pollutant Exposures in the Multi-Ethnic
##'   Study of Atherosclerosis and Air Pollution (MESA Air). Environmental Science
##'   & Technology: 43(13), 4687-4693.
##' 
##' @seealso \code{\link{estimateCV.STmodel}} and
##'   \code{\link{predictCV.STmodel}} for cross-validation.  \cr
##'   \code{\link{createSTmodel}} for creation of the originating \code{STmodel}
##'   object.
##' @family example data
NULL

### MCMC.mesa.model ###
##' The output from a Metropolis-Hastings algorithm, implemented in
##' \code{\link{MCMC.STmodel}}), run for the model in \code{\link{mesa.model}}
##' 
##' @title Example of a \code{mcmcSTmodel} structure
##' @name MCMC.mesa.model
##' @docType data
##' @format A list with elements, see the return description in
##'   \code{\link{MCMC.STmodel}}.
##' 
##' @source Contains parametere estimates for the Spatio-Temporal model applied
##'   to monitoring data from the \strong{MESA Air} project, see
##'   Cohen et.al. (2009) and \code{\link{mesa.data.raw}} for details.
##' @keywords datasets
##'
##' @example Rd_examples/Ex_MCMC_mesa_model.R
##' 
##' @references M. A. Cohen, S. D. Adar, R. W. Allen, E. Avol, C. L. Curl, T.
##'   Gould, D. Hardie, A. Ho, P. Kinney, T. V. Larson, P. D. Sampson, L.
##'   Sheppard, K. D. Stukovsky, S. S. Swan, L. S. Liu, J. D. Kaufman. (2009)
##'   Approach to Estimating Participant Pollutant Exposures in the Multi-Ethnic
##'   Study of Atherosclerosis and Air Pollution (MESA Air). Environmental Science
##'   & Technology: 43(13), 4687-4693.
##' 
##' @seealso \code{\link{createSTmodel}} for creation of the originating
##'   \code{STmodel} object.
##' @family example data
NULL
