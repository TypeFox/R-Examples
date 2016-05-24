#' Bayesian variable selection for the Pogit model
#' 
#' This function performs Bayesian variable selection for a Poisson-Logistic (Pogit) 
#' model via spike and slab priors. For posterior inference, a MCMC sampling scheme 
#' is used that relies on augmenting the observed data by the unobserved counts and 
#' involves only Gibbs sampling steps. 
#' 
#' The method provides Bayesian variable selection for regression
#' models of count data subject to under-reporting using mixture priors with a spike 
#' and a slab component. 
#' By augmenting the observed count data with the unobserved counts, the resulting 
#' model can be factorized into a Poisson and a binomial logit model part. Hence, 
#' for this two-part model, sampling algorithms for a Poisson and a binomial 
#' logit model can be used which are described in \code{\link{poissonBvs}} and 
#' \code{\link{logitBvs}}. 
#' Bayesian variable selection is incorporated in both parts of the joint model 
#' using mixture priors with a Dirac spike and (by default) a Student-t slab. 
#' The implementation relies on the representation of the respective model as a 
#' Gaussian regression model in auxiliary variables (see again the help for the 
#' respective function). Though variable selection is primarily used to identify 
#' regressors with a non-zero effect, it can also be useful for identification of
#' the Pogit model. 
#' 
#' By default, identification of the Pogit model is achieved by additional
#' information on the reporting process through validation data and 
#' incorporation of variable selection. If the information on the parameters
#' of the reporting process is not provided by validation data, the
#' identification of the model parameters has to be guaranteed by specifying an 
#' informative prior distribution (see arguments).  
#' 
#' To model under-reported clustered data, a cluster-specific random intercept can
#' be included in both model parts of the Pogit model to account for dependence
#' within clusters. Bayesian variance selection is applied to determine whether
#' there is within-cluster dependence in either part of the model.  
#' Note that an observation-specific random intercept in the
#' Poisson sub-model yields an overdispersed Pogit model for unobserved
#' heterogeneity.   
#' 
#' For details concerning the sampling algorithm see Dvorzak and Wagner (2016).
#' 
#' Details for the model specification (see arguments):  
#' \describe{ 
#'  \item{\code{model}}{\describe{\item{}{A list:}
#'    \item{\code{deltaBetafix, deltaAlphafix}}{indicator vectors of length 
#'    \code{ncol(X)-1} and \code{ncol(W)-1}, respectively, for the Poisson and the 
#'    logit sub-model, to specify which regression effects are subject to selection 
#'    (i.e., 0 = subject to selection, 1 = fix in the model); defaults to vectors 
#'    of zeros.} 
#'    \item{\code{gammaBetafix, gammaAlphafix}}{indicators for variance selection 
#'    of the random intercept term in the Poisson and the logit sub-model 
#'    (i.e., 0 = with variance selection (default), 1 = no variance selection); 
#'    only used if a random intercept is included in either part of the joint 
#'    model (see \code{riBeta} and \code{riAlpha}, respectively).} 
#'    \item{\code{riBeta, riAlpha}}{logical. If \code{TRUE}, a cluster-specific
#'    random intercept is included in the respective part of the joint model;
#'    defaults to \code{FALSE}.} 
#'    \item{\code{clBetaID, clAlphaID}}{numeric vectors of length equal to the 
#'    number of observations containing the cluster ID c = 1,...,C for each unit
#'    (or sub-category) in the respective sub-model (required if 
#'    \code{riBeta=TRUE} or \code{riAlpha=TRUE}, respectively).} 
#'    \item{\code{subcat}}{a factor variable of length equal to the number of
#'    units that specifies for which sub-category validation data are available 
#'    (is required if \code{W} is a subset of \code{X}). 
#'    If \code{NULL} (default), it is presumed that validation data are available 
#'    for each unit (see also examples).}
#' }}
#' 
#'  \item{\code{prior}}{\describe{\item{}{A list:}
#'    \item{\code{slabP, slabL}}{distribution of the slab component in the
#'    Poisson and logit sub-model, i.e. "\code{Student}" (default) or "\code{Normal}".}
#'    \item{\code{psi.nuP, psi.nuL}}{hyper-parameter of the Student-t slab in
#'    the respective sub-model (used for a Student-t slab); defaults to 5.}
#'    \item{\code{m0b, m0a}}{prior mean for the intercept parameter in the 
#'    Poisson and the logit model; defaults to 0. If the argument \code{method} = 
#'    "\code{infprior}", the specification of \code{m0a} is required.}
#'    \item{\code{M0b, M0a}}{prior variance for the intercept parameter in the
#'    Poisson and the logit model; defaults to 100.}
#'    \item{\code{bj0, aj0}}{a vector of prior means for the regression effects 
#'    in the Poisson and the logit sub-model (which is encoded in a normal distribution, 
#'    see notes); defaults to a vector of zeros. If the argument \code{method} = 
#'    "\code{infprior}", the specification of \code{aj0} is mandatory.}
#'    \item{\code{VP, VL}}{variance of the slab in the respective sub-model; 
#'    defaults to 5.}
#'    \item{\code{wBeta, wAlpha}}{hyper-parameters of the Beta-prior for the mixture 
#'    weights \eqn{\omega_\beta} and \eqn{\omega_\alpha} in the respective sub-model; 
#'    defaults to \code{c(wa0=1, wb0=1)}, i.e. a uniform distribution.}
#'    \item{\code{piBeta, piAlpha}}{hyper-parameters of the Beta-prior for the mixture 
#'    weights \eqn{\pi_\beta} and \eqn{\pi_\alpha} in the respective sub-model; 
#'    defaults to \code{c(pa0=1, pb0=1)}, i.e. a uniform distribution.}
#' }}
#' 
#'  \item{\code{mcmc}}{\describe{\item{}{A list:}
#'    \item{\code{M}}{number of MCMC iterations after the burn-in phase; 
#'    defaults to 8000.} 
#'    \item{\code{burnin}}{number of MCMC iterations discarded as burn-in; 
#'    defaults to 2000.} 
#'    \item{\code{thin}}{thinning parameter; defaults to 1.} 
#'    \item{\code{startsel}}{number of MCMC iterations drawn from the unrestricted
#'    model (e.g., \code{burnin/2}); defaults to 1000.} 
#'    \item{\code{verbose}}{MCMC progress report in each \code{verbose}-th 
#'    iteration step; defaults to 500. If \code{verbose=0}, no output is 
#'    generated.} 
#'    \item{\code{msave}}{returns additional output with variable 
#'    selection details (i.e. posterior samples for \eqn{\omega_\beta}, 
#'    \eqn{\omega_\alpha}, \eqn{\delta_\beta}, \eqn{\delta_\alpha}, 
#'    \eqn{\pi_\beta}, \eqn{\pi_\alpha}, \eqn{\gamma_\beta},
#'    \eqn{\gamma_\alpha}); defaults to \code{FALSE}.} 
#' }}
#' 
#'  \item{\code{start}}{\describe{\item{}{A list:}
#'    \item{\code{beta}}{a vector of length \code{ncol(X)} containing starting 
#'    values for the regression parameters \eqn{\beta} in the Poisson model part.
#'    By default, a Poisson glm is fitted to the observed counts.}
#'    \item{\code{alpha}}{a vector of length \code{ncol(W)} containing starting
#'    values for the regression parameters \eqn{\alpha} in the logit model part. 
#'    By default, a binomial glm is fitted to the validation data for 
#'    \code{method} = "\code{val}". If \code{method} = "\code{infprior}", 
#'    starting values for \eqn{\alpha} are sampled from the (informative) prior 
#'    distribution.}
#'    \item{\code{firth}}{logical. If \code{TRUE}, a logistic regression model 
#'    applying Firth's correction to the likelihood using 
#'    \code{\link[logistf]{logistf}} is fitted to the validation data
#'    (only used if \code{method} = "\code{val}").}
#' }}}
#' 
#' 
#' @param y an integer vector of observed counts for units i = 1,...,I
#' @param E an (optional) vector containing total exposure times (offset); 
#'  should be \code{NULL} or an integer vector of length equal to the number of 
#'  counts.
#' @param X a design matrix in the Poisson part of the joint model 
#' @param W a design matrix in the logit part of the joint model (can be a subset
#'  of \code{X}) or \code{NULL}, if the same design matrix is used in both 
#'  sub-models, i.e. \code{W} = \code{X}. 
#' @param validation a two-column data frame or list with the number of 
#'  reported cases (= \code{v}) in the validation sample and the number of 
#'  total cases (= \code{m}) subject to the fallible reporting process 
#'  (i.e. validation sample size) for each unit (or sub-category); 
#'  required if \code{method =} "\code{val}", otherwise \code{NULL}. 
#'  The number of rows must conform with the number of rows in \code{W} or with 
#'  the number of units I (if \code{X = W}), respectively. 
#' @param method the method to be used to obtain parameter identification: 
#'  The default method "\code{val}" requires a small sample of validation data 
#'  (see \code{validation}). If the information on all or some parameters
#'  of the reporting process is not provided by validation data, an informative
#'  prior distribution for the regression effects in the logit sub-model 
#'  can be used (\code{method} = "\code{infprior}"). This prior information is encoded 
#'  in a normal distribution instead of the spike and slab prior (see the details 
#'  for \code{prior}). 
#' @param model a list specifying the structure of the model (see details)
#' @param prior an (optional) list of prior settings and hyper-parameters 
#'   controlling the priors (see details)
#' @param mcmc an (optional) list of MCMC sampling options (see details)
#' @param start an (optional) list containing starting values for the 
#'   regression effects in both sub-models (see details)
#' @param BVS if \code{TRUE} (default), Bayesian variable selection 
#'  (in at least one part of the joint model) is performed to identify 
#'  regressors with a non-zero effect; otherwise, an unrestricted model is 
#'  estimated (without variable selection).
#'  
#' @return The function returns an object of class "\code{pogit}" with methods
#'  \code{\link{print.pogit}}, \code{\link{summary.pogit}} and 
#'  \code{\link{plot.pogit}}. 
#' 
#'  An object of class "\code{pogit}" is a list containing the following elements:
#'
#'  \item{\code{samplesL}}{a named list containing the samples from the posterior
#'    distribution of the parameters in the logit part of the joint model 
#'    (see also \code{msave}):
#'    \describe{
#'    \item{\code{alpha, thetaAlpha}}{regression coefficients \eqn{\alpha} and
#'    \eqn{\theta_\alpha}}
#'    \item{\code{pdeltaAlpha}}{P(\eqn{\delta_\alpha}=1)}
#'    \item{\code{psiAlpha}}{scale parameter \eqn{\psi_\alpha} of the slab component}
#'    \item{\code{pgammaAlpha}}{P(\eqn{\gamma_\alpha}=1)}
#'    \item{\code{ai}}{cluster-specific random intercept}
#'    }}
#'  \item{\code{samplesP}}{a named list containing the samples from the posterior
#'    distribution of the parameters in the Poisson part of the joint model 
#'    (see also \code{msave}):
#'    \describe{
#'    \item{\code{beta, thetaBeta}}{regression coefficients \eqn{\beta} and
#'    \eqn{\theta_\beta}}
#'    \item{\code{pdeltaBeta}}{P(\eqn{\delta_\beta}=1)}
#'    \item{\code{psiBeta}}{scale parameter \eqn{\psi_\beta} of the slab component}
#'    \item{\code{pgammaBeta}}{P(\eqn{\gamma_\beta}=1)}
#'    \item{\code{bi}}{cluster-specific random intercept}
#'    }}
#'  \item{\code{data}}{a list containing the data \code{y}, \code{E}, 
#'    \code{X}, \code{W}, \code{val} and \code{subcat}}
#'  \item{\code{model.logit}}{a list containing details on the model specification in
#'  the logit sub-model, see details for \code{model}}
#'  \item{\code{model.pois}}{a list containing details on the model specification in
#'  the Poisson sub-model, see details for \code{model}}
#'  \item{\code{mcmc}}{see details for \code{mcmc}}
#'  \item{\code{prior.logit}}{see details for \code{prior}}
#'  \item{\code{prior.pois}}{see details for \code{prior}}
#'  \item{\code{dur}}{a list containing the total runtime (\code{total}) 
#'    and the runtime after burn-in (\code{durM}), in seconds}
#'  \item{\code{BVS}}{see arguments}
#'  \item{\code{method}}{see arguments}
#'  \item{\code{start}}{a list containing starting values, see arguments}
#'  \item{\code{family}}{"pogit"}
#'  \item{\code{call}}{function call}
#' 
#' @note If \code{method} = "\code{infprior}", an 
#' informative prior for the regression parameters in the logit model is 
#' required to guarantee identification of the model parameters. 
#' Otherwise, identification of the Pogit model may be weak and inference
#' will be biased.
#' 
#' @seealso \code{\link{logitBvs}}, \code{\link{poissonBvs}} 
#' @keywords models
#'  
#' @references Dvorzak, M. and Wagner, H. (2016). Sparse Bayesian modelling
#'  of underreported count data. \emph{Statistical Modelling}, \strong{16}(1),
#'  24 - 46, \url{http://dx.doi.org/10.1177/1471082x15588398}.
#'
#' @author Michaela Dvorzak <m.dvorzak@@gmx.at>, Helga Wagner
#' @importFrom BayesLogit compute.mixture
#' @importFrom plyr rbind.fill
#' @import stats
#' @import utils
#' @export
#' 
#' @examples
#' \dontrun{
#' ## Examples below (except for m3) should take 3-4 minutes. 
#' 
#' ## ------ (use simul1) ------
#' # load simulated data set 'simul1'
#' data(simul1)
#' y <- simul1$y
#' E <- simul1$E
#' X <- as.matrix(simul1[, -c(1,2,8,9)]) # W = X
#' validation <- simul1[, c("m", "v"), drop = FALSE]
#' 
#' # function call (with specific MCMC settings)
#' m1 <- pogitBvs(y = y, E = E, X = X, validation = validation, 
#'                mcmc = list(M = 4000, thin = 5, verbose = 1000))
#'
#' # print, summarize and plot results
#' print(m1)
#' summary(m1)
#' plot(m1)
#' 
#' # show traceplots disregarding burn-in and thinning
#' plot(m1, burnin = FALSE, thin = FALSE)
#' # show density plot of MCMC draws
#' plot(m1, type = "density")
#' 
#' # informative prior instead of validation data (change prior settings)
#' # e.g. available prior information on reporting probabilities 
#' p.a0 <- 0.9 
#' p.a  <- c(0.125, 0.5, 0.5, 0.5)
#' m0a_inf <- log(p.a0/(1 - p.a0))  # prior information for alpha_0
#' aj0_inf <- log(p.a/(1 - p.a))    # prior information for alpha
#' 
#' prior.set <- list(m0a = m0a_inf, aj0 = aj0_inf, VL = 0.005, slabL = "Normal")
#' m2 <- pogitBvs(y = y, E = E, X = X, method = "infprior", prior = prior.set, 
#'                mcmc = list(M = 4000, burnin = 2000, thin = 2), BVS = FALSE)
#' print(m2)
#' summary(m2)
#' plot(m2)
#' plot(m2, type = "acf", lag.max = 30)
#'
#' ## ------ (use simul2) ------
#' # complex model (with a long (!) runtime)
#' 
#' # load simulated data set 'simul2'
#' data(simul2)
#' y <- simul2$y
#' E <- simul2$E
#' cID <- simul2$cID
#' X <- as.matrix(simul2[, -c(1:3,9,10)])
#' validation <- simul2[, c("v", "m"), drop = FALSE]
#'  
#' # function call (with random intercept in both sub-models)
#' model <- list(riBeta = 1, riAlpha = 1, clBetaID = cID, clAlphaID = cID)
#' m3 <- pogitBvs(y = y, E = E, X = X, validation = validation, model = model, 
#'                mcmc = list(M = 6000, burnin = 200, thin = 10), BVS = TRUE)
#' print(m3)
#' summary(m3)
#' plot(m3)
#'
#' ## ------ (use cervical cancer data) ------
#' # load cervical cancer data
#' data(cervical)
#' data(cervical_validation)
#' y <- cervical$y
#' E <- cervical$E
#' X <- as.matrix(cervical[, -c(1:4)])
#' validation <- cervical_validation[, c(1, 2), drop = FALSE]
#' W          <- as.matrix(cervical_validation[, -c(1:3)])
#' subcat     <- factor(as.numeric(cervical$country))
#' 
#' # function call 
#' m4 <- pogitBvs(y = y, E = E, X = X, W = W, validation = validation, 
#'                model = list(subcat = subcat), mcmc = list(M = 10000, 
#'                burnin = 2000, thin = 10), start = list(firth = TRUE), 
#'                BVS = TRUE)             
#' print(m4)
#' # additionally compute estimated risks and reporting probabilities
#' summary(m4, printRes = TRUE) 
#' plot(m4, burnin = FALSE, thin = FALSE)
#' plot(m4, type = "acf", lag.max = 50)
#' 
#' # informative prior instead of validation data (change prior settings)
#' # e.g. prior information on country-specific reporting probabilities 
#' p.a0 <- 0.85
#' p.a  <- c(0.99, 0.70, 0.85)
#' m0a_inf <- log(p.a0/(1 - p.a0))  # prior information for alpha_0
#' aj0_inf <- log(p.a/(1 - p.a))    # prior information for alpha
#' 
#' prior.set <- list(m0a = m0a_inf, aj0 = aj0_inf, VL = 0.005, slabL = "Normal")
#' m5 <- pogitBvs(y = y, X = X, W = W, E = E, method = "infprior", 
#'                model = list(subcat = subcat), prior = prior.set, 
#'                mcmc = list(M = 10000, burnin = 2000, thin = 10))
#' print(m5)
#' summary(m5, printRes = TRUE)
#' plot(m5, burnin = FALSE, thin = FALSE)
#' plot(m5, type = "acf", lag.max = 50)
#' }

pogitBvs <- function(y, E = NULL, X, W = NULL, validation = NULL, 
                     method = "val", model = list(), prior = list(),
                     mcmc = list(), start = list(), BVS = TRUE){
  
  starttime <- proc.time()[3]
  cl <- match.call()
  
  #### --- Check data and get dimensions
  if (is.null(y) || is.null(X)) stop("need 'y' and 'X' argument")
  if (!is.vector(y)) stop("number of observed counts 'y' must be a vector")
  if (!is.matrix(X)) stop("'X' must be a matrix")
  if (any(is.na(y))) stop("NA values in 'y' not allowed")
  if (any(is.na(X))) stop("NA values in 'X' not allowed")
  if (any(y < 0)) stop("number of observed counts 'y' must be positive")
  if (length(y) != nrow(X)) stop("'y' and 'nrow(X)' must have same length")

  if (is.null(E)) E <- rep(1, length(y))
  if (!is.vector(E)) stop("'E' must be vector")
  if (any(is.na(E))) stop("NA values in 'E' not allowed")
  if (any(E < 0)) stop("'E' must be positive")
  if (length(E) != length(y))  stop("'y' and 'E' must have same length")

  if (is.null(W)) W <- X
  if (!is.matrix(W)) stop("'W' must be a matrix")
  if (any(is.na(W))) stop("NA values in 'W' not allowed")
  
  if (!(method %in% c("val", "infprior"))){
    stop("invalid 'method' argument: only 'val' and 'infprior' are supported")
  }
  
  y <- as.integer(y)
  E <- as.integer(E)
  if (!all(X[,1] == 1)) X <- cbind(rep(1, dim(X)[1]), X)    
  if (!all(W[,1] == 1)) W <- cbind(rep(1, dim(W)[1]), W) 
  
  I <- length(y)
  dB <- ncol(X) - 1
  dA <- ncol(W) - 1

  colnames(X) <- paste("X", seq(0,dB), sep = "")
  colnames(W) <- paste("W", seq(0,dA), sep = "")
  

  ## validation data
  if (method == "val"){
    if (is.null(validation)){
      stop(paste(strwrap(paste("specify 'validation' data or use method 'infprior'."), 
                         exdent = 1), collapse = "\n"))
    }
    if (!is.data.frame(validation)){
      if (is.matrix(validation)){
        validation <- as.data.frame(validation)
      } 
      if (is.list(validation)){
        validation <- data.frame(v = validation$v, m = validation$m)
      } else stop("'validation' must be a data.frame")
    }

    if (any(is.na(validation))) stop("NA values in 'validation' not allowed")
    if (ncol(validation) != 2){
      stop("'validation must be da data frame with 2 columns'")
    }
    if (any(!(colnames(validation) %in% c("m", "v")))){
      if (any(validation[,1] > validation[, 2])){
        colnames(validation) <- c("m", "v")
      } else {
        colnames(validation) <- c("v", "m")
      }
    }

    if (!all(validation$v <= validation$m)){
      stop("invalid specification of 'validation' (column m must be >/= v)")
    }
    
    if (nrow(validation) != nrow(W)){
      stop("invalid specification of 'validation' sample or design matrix 'W'")
    }

    rownames(validation) <- as.character(seq_len(nrow(validation)))
    
  } 
  
  #### --- Check and update parameter lists
  
  #### (Pogit) MODEL structure
  defaultModel <- list(
    deltaBetafix  = matrix(0,1,dB), # effects subject to selection in Poisson sub-model
    deltaAlphafix = matrix(0,1,dA), # effects subject to selection in logit sub-model  
    gammaBetafix  = 0,              # random intercept selection in Poisson sub-model
    gammaAlphafix = 0,              # random intercept selection in logit sub-model
    riBeta        = 0,              # random intercept in Poisson sub-model
    riAlpha       = 0,              # random intercept in logit sub-model
    clBetaID      = NULL,           # cluster index c = 1,...,C (Poisson sub-model)
    clAlphaID     = NULL,           # cluster index c = 1,...,C (logit sub-model)
    subcat        = as.factor(seq(1,I)) 
    # factor variable to specify for which subcategory validation data is available
  )
  model <- modifyList(defaultModel, as.list(model), keep.null = TRUE)
  
  ## check parameters
  model$deltaBetafix  <- as.matrix(model$deltaBetafix)
  model$deltaAlphafix <- as.matrix(model$deltaAlphafix)
  model$gammaBetafix  <- as.matrix(model$gammaBetafix) 
  model$gammaAlphafix <- as.matrix(model$gammaAlphafix) 
  model$riBeta        <- as.integer(model$riBeta)
  model$riAlpha       <- as.integer(model$riAlpha)
  
  if (dB > 0){
    if (length(model$deltaBetafix) != dB || !all(model$deltaBetafix %in% c(0, 1))){
      stop("invalid specification of effects subject to selection in Poisson sub-model")
    }
  } else {
    if (length(model$deltaBetafix) != dB){
      stop("invalid 'deltaBetafix' argument")
    }
  }
  if (length(model$gammaBetafix)>1 || !(model$gammaBetafix %in% c(0, 1)) || !(model$riBeta %in% c(0, 1)))
    stop("invalid specification of random intercept selection in Poisson sub-model")
  
  if (dA > 0){
    if (length(model$deltaAlphafix) != dA || !all(model$deltaAlphafix %in% c(0, 1))){
      stop("invalid specification of effects subject to selection in logit sub-model")
    }
  } else {
    if (length(model$deltaAlphafix) != dA){
      stop("invalid 'deltaAlphafix' argument")
    }
  }
  if (length(model$gammaAlphafix) > 1 || !(model$gammaAlphafix %in% c(0, 1)) || !(model$riAlpha %in% c(0, 1)))
    stop("invalid specification of random intercept selection in logit sub-model")  
  

  ## validation data available only for sub-category 
  if (is.null(model$subcat)){
    stop(paste(strwrap(paste("need factor variable 'subcat' to specify for which
                               sub-category validation data is available"), 
                         exdent = 1), collapse = "\n"))
  }  
  if (!is.factor(model$subcat)){
    if (is.numeric(model$subcat)){
      model$subcat <- as.factor(model$subcat)
    } else stop("'invalid' specification of 'subcat'")
  }
  if (length(model$subcat) != I) 
    stop("invalid specification of 'subcat'")
  
  if (method == "infprior"){
    if (!is.null(model$subcat)){
      validation <- as.data.frame(matrix(0, nlevels(model$subcat), 2))
    } else {
      validation <- as.data.frame(matrix(0, I, 2))
    }
    colnames(validation) <- c("v", "m")
  }
  
  ycat <- aggregate(y, list(model$subcat), sum)[, 2]
  if (length(ycat) != nrow(validation)){
    stop("invalid specification of 'subcat'")
  } else {
    ytilde <- ycat + validation$v
  }

  
  ## (Poisson) sub-model structure 
  model.pois <- list(
    deltafix  = model$deltaBetafix,
    gammafix  = model$gammaBetafix,
    ri        = model$riBeta,
    clusterID = model$clBetaID,
    d         = dB,
    family    = "pogit", 
    Zp        = seq_len(I)
  )
  
  ## (Logit) sub-model structure 
  model.logit <- list(
    deltafix  = model$deltaAlphafix,
    gammafix  = model$gammaAlphafix,
    ri        = model$riAlpha,
    clusterID = model$clAlphaID,
    d         = dA,
    family    = "pogit",
    Zl        = seq_len(length(ytilde))
    )
  

  if (model.pois$ri == 1){
    if (!is.vector(model.pois$clusterID)){
      stop(paste(strwrap(paste("you must specify a cluster ID for random
                                intercept selection in the Poisson sub-model"), 
                          exdent = 1), collapse = "\n"))
    }
    if (min(model.pois$clusterID) != 1 || !(is.integer(model.pois$clusterID))){
      stop("specify the cluster ID as c = 1,...,C")
    }
    # random intercept model specification
    if (length(model.pois$clusterID) != I){
      stop("'y' and 'clusterID' must have same length")
    }
    model.pois$Zp <- model.pois$clusterID
    nCp <- max(model.pois$Zp)
    btilde <- rnorm(nCp, 0, 1)
    Hp <- matrix(0, I, nCp)
    for (i in 1:nCp){
      Hp[which(model.pois$Zp==i), i] <- 1
    }
  } else {
    btilde <- Hp <- NULL        
  }

  # specify function call for mcmc
  fun_select <- "select_poisson"
  txt_fun <- "Pogit"
  if (model.pois$ri == 1){
    if (nCp == length(y)){
      fun_select <- "select_poissonOD"
      txt_fun <- "overdispersed Pogit"
    }
  }

  if (model.logit$ri == 1){
    if (method == "infprior"){
      warning(simpleWarning(paste(strwrap(paste("random intercept term in the
                                                logit sub-model is not supported
                                                for method='infprior'"), 
                         exdent = 1), collapse = "\n")))
      model.logit$ri <- 0
    }
    if (!is.vector(model.logit$clusterID)){
      stop(paste(strwrap(paste("you must specify a cluster ID for random
                             intercept selection in the logit sub-model"), 
                       exdent = 1), collapse = "\n"))
    }
    if (min(model.logit$clusterID) != 1 || !(is.integer(model.logit$clusterID))){
      stop("specify the cluster ID as c = 1,...,C")
    }
    # random intercept model specification
    if (length(model.logit$clusterID) != length(ytilde)){
      stop("'y' and 'clusterID' must have same length")
    }
    model.logit$Zl <- model.logit$clusterID
    nCl <- max(model.logit$Zl)
    atilde <- rnorm(nCl, 0, 1)
    Hl <- matrix(0, length(ytilde), nCl)
    for (i in 1:nCl){
      Hl[which(model.logit$Zl==i), i] <- 1
    }
  } else {
    atilde <- Hl <- NULL        
  }  

  dfB   <- dB + 1 
  deffB <- dB + model.pois$ri
  dallB <- deffB + 1

  dfA   <- dA + 1 
  deffA <- dA + model.logit$ri
  dallA <- deffA + 1  


#### PRIOR settings and hyper-parameters controlling the priors
  defaultPrior <- list(
    slabP = "Student",        # slab distribution (Poisson sub-model)
    slabL = "Student",        # slab distribution (logit sub-model)
    psi.nuP = 5,              # hyper-parameter of slab variance (Poisson sub-model)
    psi.nuL = 5,              # hyper-parameter of slab variance (logit sub-model)
    VP = 5,                   # slab variance (Poisson sub-model)
    VL = 5,                   # slab variance (logit sub-model)
    M0b = 100,                # prior variance of intercept parameter
    M0a = 100,                # prior variance of intercept parameter
    m0b = 0,                  # prior mean of beta_0 (only relevant if slab="Normal")  
    m0a = 0,                  # prior mean of alpha_0 (only relevant if slab="Normal")  
    bj0 = rep(0, deffB),      # prior mean of regression effects beta_j 
    aj0 = rep(0, deffA),      # prior mean of regression effects alpha_j     
    wBeta   = c(wa0=1,wb0=1), # w_beta ~ Beta(wa0,wb0)
    wAlpha  = c(wa0=1,wb0=1), # w_alpha ~ Beta(wa0,wb0)
    piBeta  = c(pa0=1,pb0=1), # pi_beta ~ Beta(pa0,pb0)
    piAlpha = c(pa0=1,pb0=1)  # pi_alpha ~ Beta(pa0,pb0)
  )

  prior <- modifyList(defaultPrior, as.list(prior))


  ## (Poisson) sub-model prior structure 
  prior.pois <- list(
    slab = prior$slabP,
    psi.nu = prior$psi.nuP,
    V = prior$VP,
    M0 = prior$M0b,
    m0 = prior$m0b,
    aj0 = prior$bj0,
    w = prior$wBeta,
    pi = prior$piBeta
  )

  if (!all(prior.pois$aj0 == 0)){
    prior.pois$slab <- "Normal"
    model.pois$deltafix <- matrix(1, 1, model.pois$d)
    if (model.pois$ri) model.pois$gammafix <- 1
  }
  
  if (prior.pois$slab == "Normal"){
    prior.pois$psi.Q  <- prior.pois$V
    prior.pois$psi.nu <- NULL     
  } else {
    prior.pois$psi.Q <- prior.pois$V*(prior.pois$psi.nu - 1)
    prior.pois$aj0   <- rep(0, deffB) 
  }
  prior.pois$a0 <- matrix(c(prior.pois$m0, prior.pois$aj0), ncol = 1)
  prior.pois$invM0 <- solve(prior.pois$M0)      

  ## check and update paramaters
  with(prior.pois, 
       stopifnot(all(w > 0),
                 all(pi > 0),
                 V > 0,
                 M0 > 0, 
                 psi.Q > 0,
                 slab %in% c("Student", "Normal")     
                 )
       )

  if (prior.pois$slab == "Student") with(prior.pois, stopifnot(psi.nu > 0))
  if (length(prior.pois$aj0) != deffB)
    stop(paste(strwrap(paste("invalid specification of prior means 'bj0' for 
                             regression effects"), exdent = 1), collapse = "\n"))
  if (model.pois$ri==1) prior.pois$invB0 <- diag(nCp)


  ## (logit) sub-model prior structure 
  prior.logit <- list(
    slab = prior$slabL,
    psi.nu = prior$psi.nuL,
    V = prior$VL,
    M0 = prior$M0a,
    m0 = prior$m0a,
    aj0 = prior$aj0,
    w = prior$wAlpha,
    pi = prior$piAlpha
    )

  # no validation data available for Pogit model
  if (method == "infprior"){
    prior.logit$M0 <- prior.logit$V
    if (prior.logit$slab == "Student"){
      warning(simpleWarning(paste("prior information is encoded in a 'N(a0,A0)' distribution")))  
      prior.logit$slab <- "Normal"
    }
    #if (any(c(model.logit$deltafix, model.logit$gammafix)==0)){
    #  warning(simpleWarning(paste("BVS is not intended for method='infprior'")))
      model.logit$deltafix <- matrix(1, 1, model.logit$d)
      model.logit$gammafix <- 1                     
    #}
    
    if (!(prior.logit$V < 0.01) || all(prior.logit$aj0 == 0)){
      warning(simpleWarning(paste(strwrap(paste("This method requires an informative prior 
                                                to achieve parameter identification if
                                                no validation data is available."),
                                          exdent = 1), collapse = "\n")))
    }
  } 

  # prior information + validation data
  if (method == "validation" &&  !all(prior.logit$aj0 == 0)){
    prior.logit$slab <- "Normal"
    model.logit$deltafix <- matrix(1, 1, model.logit$d)
    if (model.logit$ri) model.logit$gammafix <- 1
  }

  if (prior.logit$slab == "Normal"){
    prior.logit$psi.Q  <- prior.logit$V
    prior.logit$psi.nu <- NULL      # removes psi.nu from list
  } else {
    prior.logit$psi.Q <- prior.logit$V*(prior.logit$psi.nu - 1)
    prior.logit$aj0   <- rep(0, deffA) 
  }
  prior.logit$a0 <- matrix(c(prior.logit$m0, prior.logit$aj0), ncol=1)
  prior.logit$invM0 <- solve(prior.logit$M0)      

  ## check and update paramaters
  with(prior.logit, 
       stopifnot(all(w > 0),
                 all(pi > 0),
                 V > 0, 
                 M0 > 0, 
                 psi.Q > 0,
                 slab %in% c("Student", "Normal")     
                 )
       )

  if (prior.logit$slab == "Student") with(prior.logit, stopifnot(psi.nu > 0))
  if (length(prior.logit$aj0) != deffA){
    stop(paste(strwrap(paste("invalid specification of prior means 'aj0' for 
                             regression effects"), exdent = 1), collapse = "\n"))
  }
  if (model.logit$ri == 1) prior.logit$invB0 <- diag(nCl)


  #### MCMC sampling options
  defaultMCMC <- list(
    M = 8000,         # number of MCMC iterations after burn-in
    burnin = 2000,    # number of MCMC iterations discared as burn-in
    startsel = 1000,  # number of MCMC iterations drawn from the unrestricted model
    thin = 1,         # saves only every thin'th draw
    verbose = 500,    # progress report in each verbose'th iteration step (if 0, no output is generated)
    msave = FALSE     # returns additional output with variable selection details
  )
  
  mcmc <- modifyList(defaultMCMC, as.list(mcmc))
  
  ## check parameters
  mcmc$M        <- as.integer(mcmc$M)
  mcmc$burnin   <- as.integer(mcmc$burnin)
  mcmc$startsel <- as.integer(mcmc$startsel)
  mcmc$thin     <- as.integer(mcmc$thin)
  mcmc$verbose  <- as.integer(mcmc$verbose)
  mcmc$nmc      <- with(mcmc, M + burnin)
  
  with(mcmc,
       stopifnot(all(c(M, burnin, startsel, thin, verbose) >= c(1, 0, 1, 1, 0)),
                 typeof(msave) == "logical")
  )

 
  if (typeof(BVS) != "logical") stop("invalid 'BVS' argument")  
  if (BVS && mcmc$startsel > mcmc$nmc){
    warning(simpleWarning(paste("invalid 'BVS' or 'startsel' argument:\n",
                                "'startsel' is set to ", mcmc$burnin/2, sep = "")))
    mcmc$startsel <- mcmc$burnin/2 
  }

  selP <- with(model.pois, sum(sum(deltafix) + sum(gammafix)) == sum(d + ri))
  selL <- with(model.logit, sum(sum(deltafix) + sum(gammafix)) == sum(d + ri))
  if (BVS && method=="val" && any(c(selP, selL))){ 
    warning(simpleWarning(paste("invalid 'BVS' argument or selections specification:\n",
                                "'BVS' is set to ", FALSE, sep = "")))    
    BVS <- FALSE
  }
  if (BVS && method == "infprior" && selP){
    warning(simpleWarning(paste("invalid 'BVS' argument or selections specification:\n",
                                "'BVS' is set to ", FALSE, sep = "")))    
    BVS <- FALSE    
  }
  if (!BVS){
    mcmc$startsel  <- mcmc$M + mcmc$burnin + 1
    model.pois$deltafix <- matrix(1, 1, model.pois$d)    
    model.logit$deltafix <- matrix(1, 1, model.logit$d)    
    if (model.pois$ri)  model.pois$gammafix <- 1 
    if (model.logit$ri) model.logit$gammafix <- 1 
    txt.verbose   <- ""
    mcmc$msave    <- FALSE
  } else {
    txt.verbose <- " with variable selection" 
  }
  

  #### STARTING values
  if (is.null(start$beta)){
    betaStart <- glm(y ~ X[, -1, drop=FALSE] + offset(log(E)), family="poisson")
    beta <- as.numeric(betaStart$coefficients)
  } else beta <- NULL

  if (is.null(start$alpha)){
    if (method == "val"){
        alphaStart <- glm(cbind(v, m - v) ~ W[, -1, drop=FALSE], 
                          family = binomial(), data = validation)
        alpha <- as.numeric(alphaStart$coefficients)
    } else { # sample starting values from the informative prior
      alpha <- rnorm(dallA, prior.logit$a0, sqrt(prior.logit$V)) 
    }
  } else alpha <- NULL
  

  #### STARTING values (optional)
  defaultStart <- list(
    beta     = beta,
    alpha    = alpha,
    firth    = FALSE # use of Firth's penalized maximum likelihood (separation)
  )

  start <- modifyList(defaultStart, as.list(start))

  # separation in the logit model
  # uses Firth's penalized-likelihood logistic regression from package 'logistf'
  if (start$firth && method == "val"){ 
    if (length(unique(validation$m)) == 1){
      bin <- c(apply(validation, 1, function(x){
        c(rep(1, x[1]), rep(0, x[2] - x[1]))
      }))  
    } else {
      bin <- stack(apply(validation, 1, function(x){
        c(rep(1, x[1]), rep(0, x[2] - x[1]))
      }))
    }
    Wbin <- apply(W, 2, rep, times = validation$m)
    alphaStart <- logistf::logistf(values ~ Wbin[, -1, drop = FALSE], 
                                   firth = TRUE, data = bin)
    start$alpha <- as.numeric(alphaStart$coefficients)
  } else start$firth <- FALSE

  ## check starting values
  if(length(start$alpha) != (model.logit$d + 1) || !(all(is.numeric(start$alpha)))){
   stop("invalid specification of starting values 'alpha'") 
  }
  if(length(start$beta) != (model.pois$d + 1) || !(all(is.numeric(start$beta)))){
    stop("invalid specification of starting values 'beta'") 
  }

  ## fixed starting settings
  start$deltaBeta  <- start$pdeltaBeta <- matrix(1, 1, model.pois$d)
  start$deltaAlpha <- start$pdeltaAlpha <- matrix(1, 1, model.logit$d)
  start$psiBeta    <- matrix(prior.pois$V, 1, deffB)  
  start$psiAlpha   <- matrix(prior.logit$V, 1, deffA)  
  start$omegaBeta <- start$omegaAlpha <- 1
  if(model.pois$ri == 1){
    start$gammaBeta <- start$pgammaBeta <- matrix(1, 1, model.pois$ri)   
    start$thetaBeta <- rnorm(model.pois$ri, 0, 0.1)
    start$piBeta    <- matrix(1, 1, model.pois$ri)
  } else {
    start$gammaBeta <- start$pgammaBeta <- start$thetaBeta <- start$piBeta <- NULL 
  }
  if(model.logit$ri == 1){
    start$gammaAlpha <- start$pgammaAlpha <- matrix(1, 1, model.logit$ri)   
    start$thetaAlpha <- rnorm(model.logit$ri, 0, 0.1)
    start$piAlpha    <- matrix(1, 1, model.logit$ri)
  } else {
    start$gammaAlpha <- start$pgammaAlpha <- start$thetaAlpha <- start$piAlpha <- NULL 
  }


  lambda <- exp(X%*%start$beta)
  p  <- exp(W%*%start$alpha)/(1 + exp(W%*%start$alpha))
  pS <- p[model$subcat]

  par.pois <- list(
    beta      = start$beta,
    delta     = start$deltaBeta,
    pdelta    = start$pdeltaBeta,
    omega     = start$omegaBeta, 
    psi       = start$psiBeta,
    btilde    = btilde,
    theta     = start$thetaBeta,
    gamma     = start$gammaBeta,
    pgamma    = start$pgammaBeta,
    pi        = start$piBeta)

  par.logit <- list(
    alpha     = start$alpha,
    delta     = start$deltaAlpha,
    pdelta    = start$pdeltaAlpha,
    omega     = start$omegaAlpha, 
    psi       = start$psiAlpha,
    atilde    = atilde,
    theta     = start$thetaAlpha,
    gamma     = start$gammaAlpha,
    pgamma    = start$pgammaAlpha,
    pi        = start$piAlpha)


  #### SAVE --- matrices to save the MCMC draws
  ## Poisson part
  beta <- matrix(0, mcmc$nmc, model.pois$d + 1)
  colnames(beta) <- paste("beta", seq(0, model.pois$d), sep = ".")

  if (model.pois$d > 0){
    pdeltaBeta <- matrix(0, mcmc$nmc, model.pois$d)
    colnames(pdeltaBeta) <- paste("pdeltaB", seq_len(model.pois$d), sep = ".")
  } else {
    pdeltaBeta <- NULL
  }
  if (model.pois$d > 0 && prior.pois$slab == "Student"){
    psiBeta   <- matrix(0, mcmc$nmc, deffB)
    colnames(psiBeta)   <- paste("psiB", seq_len(deffB), sep = ".")  
  } else {
    psiBeta <- NULL
  }

  if (mcmc$msave && model.pois$d > 0){
    omegaBeta <- matrix(0, mcmc$nmc, 1)
    colnames(omegaBeta) <- "omegaBeta"
    deltaBeta <- matrix(0, mcmc$nmc, model.pois$d) 
    colnames(deltaBeta) <- paste("deltaB", seq_len(model.pois$d), sep = ".")
  }
  if (model.pois$ri == 1){
    pgammaBeta <- thetaBeta <- matrix(0, mcmc$nmc, model.pois$ri)
    colnames(pgammaBeta) <- "pgammaB"
    colnames(thetaBeta)  <- "thetaB"
    bi <- matrix(0, mcmc$nmc, nCp)
    colnames(bi) <- paste("b", seq_len(nCp), sep = ".")
    if (mcmc$msave){
      piBeta <- gammaBeta <- matrix(0, mcmc$nmc, model.pois$ri)
      colnames(piBeta) <- "piB"
      colnames(gammaBeta) <- "gammaB"
    } else piBeta <- gammaBeta <- NULL
  } else {
    pgammaBeta <- thetaBeta <- bi <- piBeta <- gammaBeta <- NULL
  }

  ## Logit part
  alpha <- matrix(0, mcmc$nmc, model.logit$d + 1)
  colnames(alpha) <- paste("alpha", seq(0,model.logit$d), sep = ".")

  if (model.logit$d > 0){
    pdeltaAlpha <- matrix(0, mcmc$nmc, model.logit$d)
    colnames(pdeltaAlpha) <- paste("pdeltaA", seq_len(model.logit$d), sep = ".")
  } else {
    pdeltaAlpha <- NULL
  }
  if (model.logit$d > 0 && prior.logit$slab == "Student"){
    psiAlpha   <- matrix(0, mcmc$nmc, deffA)
    colnames(psiAlpha)   <- paste("psiA", seq_len(deffA), sep = ".")  
  } else {
    psiAlpha <- NULL
  }

  if (mcmc$msave && method == "val" && model.pois$d > 0){
    omegaAlpha <- matrix(0, mcmc$nmc, 1)
    colnames(omegaAlpha) <- "omegaAlpha"
    deltaAlpha <- matrix(0, mcmc$nmc, model.logit$d) 
    colnames(deltaAlpha) <- paste("deltaA", seq_len(model.logit$d), sep = ".")
  }
  if (model.logit$ri == 1){
    pgammaAlpha <- thetaAlpha <- matrix(0, mcmc$nmc, model.logit$ri)
    colnames(pgammaAlpha) <- "pgammaA"
    colnames(thetaAlpha)  <- "thetaA"
    ai <- matrix(0, mcmc$nmc, nCl)
    colnames(ai) <- paste("a", seq_len(nCl), sep = ".")
    if (mcmc$msave && method=="val"){
      piAlpha <- gammaAlpha <- matrix(0, mcmc$nmc, model.logit$ri)
      colnames(piAlpha) <- "piA"
      colnames(gammaAlpha) <- "gammaA"
    } else piAlpha <- gammaAlpha <- NULL
  } else {
    pgammaAlpha <- thetaAlpha <- ai <- piAlpha <- gammaAlpha <- NULL
  }

  ## further Pogit model components
  urc <- matrix(0, mcmc$nmc, I)
  colnames(urc) <- paste("d", seq_len(I), sep = ".")


  ## read mixture components of Gaussian mixture approximation
  mcomp <- mixcomp_poisson()

  c1  <- rep(1,I)%*%t(log(mcomp$w[1, ]) - 0.5*log(mcomp$v[1, ]))
  cm1 <- list(comp = list(m = mcomp$m[1, ], v = mcomp$v[1, ], w = mcomp$w[1, ]), 
              c1 = c1)


  #### ------------------------------ MCMC ------------------------------
  for (imc in 1:mcmc$nmc){
  
    if (mcmc$verbose > 0){
      if (imc == 1) cat(paste("\nMCMC for the ", txt_fun, " model", txt.verbose, ":\n\n", sep = ""))
      if (is.element(imc, c(1:5, 10, 20, 50, 100, 200, 500))){
        cat("it =", imc, "/--- duration of MCMC so far:", 
            round(timediff <- proc.time()[3] - starttime, 2), "sec.,  expected time to end:", 
            round((timediff/(imc - 1) * mcmc$nmc - timediff)/60, 2), " min. \n")
        flush.console()
      } else if (imc %% mcmc$verbose == 0 && imc < mcmc$nmc){
        cat("it =", imc, "/--- duration of MCMC so far:", 
            round(timediff <- proc.time()[3] - starttime, 2), "sec.,  expected time to end:", 
            round((timediff/(imc - 1) * mcmc$nmc - timediff)/60, 2), " min. \n")
        flush.console()
      }
      else if (imc == mcmc$nmc) {
        cat("it =", imc, "/--- duration of MCMC (total):", 
            round(timediff <- proc.time()[3] - starttime, 2), "sec. \n \n")
        flush.console() 
      }
    }# end(verbose)

    if (imc == (mcmc$burnin + 1)) starttimeM <- proc.time()[3]
    
    ## step (1) --- sample the number of unreported cases
    d <- rpois(I, E*lambda*(1 - pS))
    n <- y + d
    ncat <- aggregate(n, list(model$subcat), sum)[, 2]
    ntilde <- ncat + validation$m
    
    
    ## step (2) --- update logit part
    par.logit <- select_logit(ytilde, ntilde, W, Hl, compmix.bin = NULL, model.logit, 
                           prior.logit, mcmc, par.logit, imc)     

    muL <- W%*%par.logit$alpha
    if (model.logit$ri == 1){
      linp <- muL + par.logit$atilde[model.logit$Zl]*par.logit$theta  
    } else linp <- muL
    
    p  <- exp(linp)/(1+exp(linp))
    pS <- p[model$subcat]
    
    
    ## step (3) --- update Poisson part
    par.pois <- do.call(fun_select, args=list(y=n, X=X, offset=E, H=Hp, 
                                              mcomp=mcomp, 
                                              compmix.pois = NULL, 
                                              cm1=cm1, model=model.pois,
                                              prior=prior.pois, mcmc=mcmc,
                                              param=par.pois, imc=imc)
    )
      
    muP <- X%*%par.pois$beta
    if (model.pois$ri == 1){
      linp <- muP + par.pois$btilde[model.pois$Zp]*par.pois$theta
    } else linp <- muP
    
    lambda <- exp(linp)
    
    
    ## saving MCMC draws
    urc[imc,] <- d

    alpha[imc,] <- par.logit$alpha
    if (model.logit$d > 0){
      pdeltaAlpha[imc,] <- par.logit$pdelta
      if (prior.logit$slab == "Student") psiAlpha[imc,] <- par.logit$psi
    }
    
    if (mcmc$msave && model.logit$d > 0){
      omegaAlpha[imc,] <- par.logit$omega
      deltaAlpha[imc,] <- par.logit$delta
    }
    
    if (model.logit$ri == 1){
      pgammaAlpha[imc] <- par.logit$pgamma
      thetaAlpha[imc]  <- par.logit$theta
      ai[imc,]         <- t(par.logit$atilde*par.logit$theta)
      
      if (mcmc$msave && method == "val"){
        piAlpha[imc] <- par.logit$pi
        gammaAlpha[imc] <- par.logit$gamma
      }
    }

    beta[imc,] <- par.pois$beta
    if (model.pois$d > 0){
      pdeltaBeta[imc, ] <- par.pois$pdelta
      if (prior.pois$slab=="Student") psiBeta[imc, ] <- par.pois$psi
    }
    
    if (mcmc$msave && model.pois$d > 0){
      omegaBeta[imc,] <- par.pois$omega
      deltaBeta[imc,] <- par.pois$delta
    }
  
    if (model.pois$ri == 1){
      pgammaBeta[imc] <- par.pois$pgamma
      thetaBeta[imc]  <- par.pois$theta
      bi[imc,]        <- t(par.pois$btilde*par.pois$theta)
      
      if (mcmc$msave){
        piBeta[imc] <- par.pois$pi
        gammaBeta[imc] <- par.pois$gamma
      }
    }
  
  }# end(MCMC)

  finish <- proc.time()[3]
  durT <- finish - starttime
  durM <- finish - starttimeM
  dur  <- list(total = durT, durM = durM) 

  samplesL <- if (mcmc$msave){
    list(alpha = alpha, pdeltaAlpha = pdeltaAlpha, psiAlpha = psiAlpha, 
         omegaAlpha = omegaAlpha, deltaAlpha = deltaAlpha, 
         thetaAlpha = thetaAlpha, pgammaAlpha = pgammaAlpha, ai = ai,
         piAlpha = piAlpha, gammaAlpha = gammaAlpha)
    } else {
      list(alpha = alpha, pdeltaAlpha = pdeltaAlpha, psiAlpha = psiAlpha, 
           thetaAlpha = thetaAlpha, pgammaAlpha = pgammaAlpha, ai = ai)
    }

  samplesP <- if (mcmc$msave){
    list(beta = beta, pdeltaBeta = pdeltaBeta, psiBeta = psiBeta, 
         omegaBeta = omegaBeta, deltaBeta = deltaBeta, 
         thetaBeta = thetaBeta, pgammaBeta = pgammaBeta, bi = bi, 
         piBeta = piBeta, gammaBeta = gammaBeta, urc = urc)
  } else {
    list(beta = beta, pdeltaBeta = pdeltaBeta, psiBeta = psiBeta,
         thetaBeta = thetaBeta, pgammaBeta = pgammaBeta, bi = bi)
  }

  ret <- list(samplesL = samplesL, samplesP = samplesP, 
              data = list(y = y, E = E, X = X, W = W, val = validation, 
                          subcat = model$subcat), model.logit = model.logit, 
              model.pois = model.pois, mcmc = mcmc, prior.logit = prior.logit, 
              prior.pois = prior.pois, dur = dur, BVS = BVS, method = method, 
              start = start, family = "pogit", call = cl, fun = fun_select) 
  class(ret) <- "pogit"
  return(ret)
}

