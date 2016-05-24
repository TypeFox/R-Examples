#' Bayesian variable selection for the Poisson model
#' 
#' This function performs Bayesian variable selection for Poisson regression models 
#' via spike and slab priors. A cluster- (or observation-) specific 
#' random intercept can be included in the model to account for within-cluster 
#' dependence (or overdispersion) with variance selection of the random intercept. 
#' For posterior inference, a MCMC sampling scheme is used
#' which relies on data augmentation and involves only Gibbs sampling steps.
#' 
#' The method provides a Bayesian framework for variable selection in regression
#' modelling of count data using mixture priors with a spike and a slab 
#' component to identify regressors with a non-zero effect. More specifically, a
#' Dirac spike is used, i.e. a point mass at zero, and (by default), the slab 
#' component is specified as a scale mixture of normal distributions, resulting 
#' in a Student-t distribution with 2\code{psi.nu} degrees of freedom. 
#' In the more general random intercept model, variance selection of the random 
#' intercept is based on the non-centered parameterization of the model, where 
#' the signed standard deviation \eqn{\theta_\beta} of the random intercept term 
#' appears as a further regression effect in the model equation. For further 
#' details, see Wagner and Duller (2012).
#' 
#' The implementation of Bayesian variable selection further relies on the 
#' representation of the Poisson model as a Gaussian regression model in 
#' auxiliary variables. Data augmentation is based on the auxiliary mixture 
#' sampling algorithm of Fruehwirth-Schnatter et al. (2009), where the
#' inter-arrival times of an assumed Poisson process are introduced as latent
#' variables. The error distribution, a negative log-Gamma distribution, 
#' in the auxiliary model is approximated by a finite mixture of normal 
#' distributions where the mixture parameters of the matlab package 
#' \code{bayesf}, Version 2.0 of Fruehwirth-Schnatter (2007) are used.  
#' See Fruehwirth-Schnatter et al. (2009) for details.
#' 
#' For details concerning the sampling algorithm, see Dvorzak and Wagner (2016) 
#' and Wagner and Duller (2012). 
#' 
#' Details for the model specification (see arguments): 
#' \describe{ 
#'  \item{\code{model}}{\describe{\item{}{A list:}  
#'    \item{\code{deltafix}}{an indicator vector of length 
#'    \code{ncol(X)-1} specifying which regression effects are subject to selection
#'    (i.e., 0 = subject to selection, 1 = fix in the model); defaults to a vector 
#'    of zeros.} 
#'    \item{\code{gammafix}}{an indicator for variance selection of 
#'    the random intercept term (i.e., 0 = with variance selection (default), 1 = no 
#'    variance selection); only used if a random intercept is includued in the model 
#'    (see \code{ri}).} 
#'    \item{\code{ri}}{logical. If \code{TRUE}, a cluster- (or observation-) specific 
#'    random intercept is included in the model; defaults to \code{FALSE}.} 
#'    \item{\code{clusterID}}{a numeric vector of length equal to the number
#'    of observations containing the cluster ID c = 1,...,C for each observation 
#'    (required if \code{ri=TRUE}). Note that \code{seq_along(y)} specifies an
#'    overdispersed Poisson model with observation-specific (normal) random intercept 
#'    (see note).} 
#'  }}
#'  \item{\code{prior}}{\describe{\item{}{A list:}      
#'    \item{\code{slab}}{distribution of the slab component, i.e. "\code{Student}" 
#'    (default) or "\code{Normal}".} 
#'    \item{\code{psi.nu}}{hyper-parameter of the Student-t slab (used for a 
#'    "\code{Student}" slab); defaults to 5.} 
#'    \item{\code{m0}}{prior mean for the intercept parameter; defaults to 0.} 
#'    \item{\code{M0}}{prior variance for the intercept parameter; defaults to 100.} 
#'    \item{\code{aj0}}{a vector of prior means for the regression effects 
#'    (which is encoded in a normal distribution, see note); defaults to vector of 
#'    zeros.} 
#'    \item{\code{V}}{variance of the slab; defaults to 5.} 
#'    \item{\code{w}}{hyper-parameters of the Beta-prior for the mixture weight
#'    \eqn{\omega}; defaults to \code{c(wa0=1, wb0=1)}, i.e. a uniform distribution.} 
#'    \item{\code{pi}}{hyper-parameters of the Beta-prior for the mixture weight 
#'    \eqn{\pi}; defaults to \code{c(pa0=1, pb0=1)}, i.e. a uniform distribution.} 
#'  }}
#'  \item{\code{mcmc}}{\describe{\item{}{A list:}
#'    \item{\code{M}}{number of MCMC iterations after the burn-in phase; defaults
#'    to 8000.} 
#'    \item{\code{burnin}}{number of MCMC iterations discarded as burn-in; defaults 
#'    to 2000.} 
#'    \item{\code{thin}}{thinning parameter; defaults to 1.} 
#'    \item{\code{startsel}}{number of MCMC iterations drawn from the unrestricted 
#'    model (e.g., \code{burnin/2}); defaults to 1000.} 
#'    \item{\code{verbose}}{MCMC progress report in each \code{verbose}-th iteration 
#'    step; defaults to 500. If \code{verbose=0}, no output is generated.} 
#'    \item{\code{msave}}{returns additional output with variable selection details 
#'    (i.e. posterior samples for \eqn{\omega}, \eqn{\delta}, \eqn{\pi}, 
#'    \eqn{\gamma}); defaults to \code{FALSE}.}
#'  }}}
#' 
#' 
#' @param y an integer vector of count data
#' @param offset an (optional) offset term; should be \code{NULL} or an integer 
#'   vector of length equal to the number of counts.
#' @param X a design matrix (including an intercept term)
#' @param model an (optional) list specifying the structure of the model (see 
#'   details)
#' @param prior an (optional) list of prior settings and hyper-parameters 
#'   controlling the priors (see details)
#' @param mcmc an (optional) list of MCMC sampling options (see details)
#' @param start an (optional), numeric vector containing starting values for the
#'   regression effects (including an intercept term); defaults to \code{NULL} 
#'   (i.e. a vector of zeros is used).
#' @param BVS if \code{TRUE} (default), Bayesian variable selection is performed
#'   to identify regressors with a non-zero effect; otherwise, an unrestricted 
#'   model is estimated (without variable selection).
#'   
#' @return The function returns an object of class "\code{pogit}" with methods
#'  \code{\link{print.pogit}}, \code{\link{summary.pogit}} and 
#'  \code{\link{plot.pogit}}. 
#' 
#'  The returned object is a list containing the following elements: 
#' 
#'  \item{\code{samplesP}}{a named list containing the samples from the posterior
#'    distribution of the parameters in the Poisson model 
#'    (see also \code{msave}): 
#'    \describe{
#'    \item{\code{beta, thetaBeta}}{regression coefficients \eqn{\beta} and
#'    \eqn{\theta_\beta}}
#'    \item{\code{pdeltaBeta}}{P(\eqn{\delta_\beta}=1)}
#'    \item{\code{psiBeta}}{scale parameter \eqn{\psi_\beta} of the slab component}
#'    \item{\code{pgammaBeta}}{P(\eqn{\gamma_\beta}=1)}
#'    \item{\code{bi}}{cluster- (or observation-) specific random intercept}
#'    }}
#'  \item{\code{data}}{a list containing the data \code{y}, \code{offset} and \code{X}}
#'  \item{\code{model.pois}}{a list containing details on the model specification, 
#'  see details for \code{model}}
#'  \item{\code{mcmc}}{see details for \code{mcmc}}
#'  \item{\code{prior.pois}}{see details for \code{prior}}
#'  \item{\code{dur}}{a list containing the total runtime (\code{total}) 
#'    and the runtime after burn-in (\code{durM}), in seconds}
#'  \item{\code{BVS}}{see arguments}
#'  \item{\code{start}}{a list containing starting values, see arguments}
#'  \item{\code{family}}{"poisson"}
#'  \item{\code{call}}{function call}
#'  
#'   
#' @note If prior information on the regression parameters is available, this 
#'   information is encoded in a normal distribution instead of the spike 
#'   and slab prior (consequently, \code{BVS} is set to \code{FALSE}).
#' @note This function can also be used to accommodate overdispersion in 
#'   count data by specifying an observation-specific random intercept
#'   (see details for \code{model}). The resulting model is an alternative 
#'   to the negative binomial model, see \code{\link{negbinBvs}}. 
#'   Variance selection of the random intercept may be useful to examine 
#'   whether overdispersion is present in the data. 
#'   
#' @seealso \code{\link{pogitBvs}}, \code{\link{negbinBvs}}
#' @keywords models
#'  
#' @references Dvorzak, M. and Wagner, H. (2016). Sparse Bayesian modelling
#'  of underreported count data. \emph{Statistical Modelling}, \strong{16}(1),
#'  24 - 46, \url{http://dx.doi.org/10.1177/1471082x15588398}.
#' @references Fruehwirth-Schnatter, S. (2007). Matlab package \code{bayesf} 2.0
#'   on \emph{Finite Mixture and Markov Switching Models}, Springer. 
#'   \url{http://statmath.wu.ac.at/~fruehwirth/monographie/}.
#' @references Fruehwirth-Schnatter, S., Fruehwirth, R., Held, L. and Rue, H. 
#'   (2009). Improved auxiliary mixture sampling for hierarchical models of 
#'   non-Gaussian data. \emph{Statistics and Computing}, \strong{19}, 479 - 492.
#' @references Wagner, H. and Duller, C. (2012). Bayesian model selection for 
#'   logistic regression models with random intercept. \emph{Computational 
#'   Statistics and Data Analysis}, \strong{56}, 1256 - 1274.
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
#' ## Examples below should take about 1-2 minutes.
#' 
#' ## ------ (use simul_pois1) ------
#' # load simulated data set 'simul_pois1'
#' data(simul_pois1)
#' y <- simul_pois1$y
#' X <- as.matrix(simul_pois1[, -1])
#' 
#' # Bayesian variable selection for simulated data set
#' m1 <- poissonBvs(y = y, X = X)
#' 
#' # print, summarize and plot results
#' print(m1)
#' summary(m1)
#' plot(m1, maxPlots = 4)
#' plot(m1, burnin = FALSE, thin = FALSE, maxPlots = 4)
#' plot(m1, type = "acf")
#' 
#' # MCMC sampling without BVS with specific MCMC and prior settings
#' m2 <- poissonBvs(y = y, X = X, prior = list(slab = "Normal"), 
#'                  mcmc = list(M = 6000, thin = 10), BVS = FALSE)
#' print(m2)
#' summary(m2, IAT = TRUE)
#' plot(m2)
#' # show traceplots disregarding thinning
#' plot(m2, thin = FALSE)
#' 
#' # specification of an overdispersed Poisson model with observation-specific 
#' # (normal) random intercept
#' cID <- seq_along(y)
#' m3  <- poissonBvs(y = y, X = X, model = list(ri = TRUE, clusterID = cID))
#' 
#' # print, summarize and plot results
#' print(m3)
#' summary(m3) 
#' # note that variance selection of the random intercept indicates that 
#' # overdispersion is not present in the data
#' plot(m3, burnin = FALSE, thin = FALSE)
#' 
#' ## ------ (use simul_pois2) ------
#' # load simulated data set 'simul_pois2'
#' data(simul_pois2)
#' y <- simul_pois2$y
#' X <- as.matrix(simul_pois2[, -c(1,2)])
#' cID <- simul_pois2$cID
#' 
#' # BVS for a Poisson model with cluster-specific random intercept
#' m4 <- poissonBvs(y = y, X = X, model = list(ri = TRUE, clusterID = cID),
#'                  mcmc = list(M = 4000, burnin = 2000))
#' print(m4)
#' summary(m4)
#' plot(m4)
#'                                
#' # similar to m4, but without variance selection of the random intercept term
#' model <- list(gammafix = 1, ri = 1, clusterID = cID)
#' m5 <- poissonBvs(y = y, X = X, model = model, mcmc = list(M = 4000, thin = 5))
#' print(m5)       
#' summary(m5)          
#' plot(m5)
#' 
#' # MCMC sampling without BVS for clustered observations
#' m6 <- poissonBvs(y = y, X = X, model = list(ri = 1, clusterID = cID), 
#'                  BVS = FALSE)
#' print(m6)         
#' summary(m6)        
#' plot(m6, maxPlots = 4)
#' }

poissonBvs <- function(y, offset = NULL, X, model = list(), mcmc  = list(),
                       prior = list(), start = NULL, BVS = TRUE){
  
  starttime <- proc.time()[3]
  cl <- match.call()
  
  #### --- Check data and get dimensions
  if (is.null(y) || is.null(X)) stop("need 'y' and 'X' argument")
  if (!is.vector(y)) stop("number of counts 'y' must be a vector")
  if (!is.matrix(X)) stop("'X' must be a matrix")
  if (any(is.na(y))) stop("NA values in 'y' not allowed")
  if (any(is.na(X))) stop("NA values in 'X' not allowed")  
  if (any(y < 0)) stop("'y' must be positive")
  if (length(y)!=nrow(X)) stop("'y' and 'nrow(X)' must have same length")
  
  if (is.null(offset)) offset <- rep(1, length(y))
  if (!is.vector(offset)) stop("'offset' must be vector")
  if (any(is.na(offset))) stop("NA values in 'offset' not allowed")
  if (any(offset < 0)) stop("'offset' must be positive")
  if (length(offset) != length(y))  stop("'y' and 'offset' must have same length")
  
  y <- as.integer(y)
  offset <- as.integer(offset)
  n <- length(y)
  
  if (!all(X[,1]==1)) X <- cbind(rep(1, dim(X)[1]), X)    
  d <- ncol(X)-1
  colnames(X) <- paste("X", seq(0, d), sep="")
  
  
  #### --- Check and update parameter lists
  
  #### MODEL structure
  defaultModel <- list(
    deltafix  = matrix(0,1,d),    # effects subject to selection
    gammafix  = 0,                # random intercept selection
    ri        = 0,                # RI model for clustered observations
    clusterID = NULL              # cluster index c = 1,...,C
  )
  model <- modifyList(defaultModel, as.list(model))
  
  ## check parameters
  model$deltafix <- as.matrix(model$deltafix)
  model$gammafix <- as.matrix(model$gammafix) 
  model$ri       <- as.integer(model$ri)
  model$d        <- d 
  model$family   <- "poisson"
  
  if (model$d > 0){
    if (length(model$deltafix) != model$d || !all(model$deltafix %in% c(0, 1))){
      stop("invalid specification of effects subject to selection")
    }
  } else {
    if (length(model$deltafix) != model$d){
      stop("invalid 'deltafix' argument")
    }
  }
  
  if (length(model$gammafix) > 1 || !(model$gammafix %in% c(0, 1)) || !(model$ri %in% c(0, 1)))
    stop("invalid specification of random intercept selection")
  
  if (model$ri == 1){
    if (!is.vector(model$clusterID)){
      stop("you must specify a cluster ID for random intercept selection")
    }
    if (min(model$clusterID) != 1 || !(is.numeric(model$clusterID))){
      stop("specify the cluster ID as c = 1,...,C")
    }
    # random intercept model specification
    if (length(model$clusterID) != n){
      stop("'y' and 'clusterID' must have same length")
    }
    model$Zp <- model$clusterID
    nC <- max(model$Zp)
    btilde <- rnorm(nC, 0, 1)
    H <- matrix(0, n, nC)
    for (i in 1:nC){
      H[which(model$Zp==i), i] <- 1
    }
  } else {
    btilde <- H <- NULL       
  }
  
  # specify function call for mcmc
  fun_select <- "select_poisson"
  txt_fun <- "Poisson"
  if (model$ri == 1){
    if (nC == length(y)){
      fun_select <- "select_poissonOD"
      txt_fun <- "overdispersed Poisson"
    }
  }
  
  df   <- d + 1 
  deff <- d + model$ri
  dall <- deff + 1
  
  
  #### PRIOR settings and hyper-parameters controlling the priors
  defaultPrior <- list(
    slab = "Student",         # slab distribution
    psi.nu = 5,               # hyper-parameter of slab variance
    V = 5,                    # slab variance
    M0 = 100,                 # prior variance of intercept parameter
    m0 = 0,                   # only relevant if slab="Normal"  
    aj0 = rep(0,deff),        # only relevant if slab="Normal"    
    w = c(wa0=1,wb0=1),       # w ~ Beta(wa0,wb0)
    pi = c(pa0=1,pb0=1)       # pi ~ Beta(pa0,pb0)
  )
  
  prior <- modifyList(defaultPrior, as.list(prior))
  if (!all(prior$aj0 == 0)){
    prior$slab <- "Normal"
    if (BVS){
      warning(simpleWarning(paste(strwrap(paste("BVS is not performed if prior
                                                information on the regression
                                                parameters is available:\n",
                                                "'BVS' was set to FALSE", sep = ""), 
                                          exdent = 1), collapse = "\n")))      
      
      BVS <- FALSE
    }
  }
  
  if (prior$slab == "Normal"){
    prior$psi.Q  <- prior$V
    prior$psi.nu <- NULL      
  } else {
    prior$psi.Q <- prior$V*(prior$psi.nu - 1)
  }
  prior$a0 <- matrix(c(prior$m0, prior$aj0), ncol = 1)
  prior$invM0 <- solve(prior$M0)      
  
  ## check and update paramaters
  with(prior, 
       stopifnot(all(w > 0),
                 all(pi > 0),
                 V > 0, 
                 M0 > 0, 
                 psi.Q > 0,
                 slab %in% c("Student", "Normal")     
       )
  )
  
  if (prior$slab == "Student") with(prior, stopifnot(psi.nu > 0))
  if (length(prior$aj0) != deff)
    stop("invalid specification of prior means for regression effects")
  if (model$ri == 1) prior$invB0 <- diag(nC)
  
  
  #### MCMC sampling options
  defaultMCMC <- list(
    M = 8000,         # number of MCMC iterations after burn-in
    burnin = 2000,    # number of MCMC iterations discared as burn-in
    startsel = 1000,  # number of MCMC iterations drawn from the unrestricted model
    thin = 1,         # saves only every thin'th draw
    verbose = 500,    # progress report in each verbose'th iteration step 
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
                                "'startsel' was set to ", mcmc$burnin/2, sep = "")))
    mcmc$startsel <- mcmc$burnin/2 
  }
  if (BVS && sum(sum(model$deltafix) + sum(model$gammafix)) == sum(model$d + model$ri)){
    warning(simpleWarning(paste("invalid 'BVS' argument:\n",
                                "'BVS' was set to ", FALSE, sep = "")))
    BVS <- FALSE
  }
  if (!BVS){
    mcmc$startsel  <- mcmc$M + mcmc$burnin + 1
    model$deltafix <- matrix(1, 1, model$d)
    if (model$ri) model$gammafix <- 1 
    txt.verbose   <- ""
    mcmc$msave    <- FALSE
  } else {
    txt.verbose <- " with variable selection" 
  }
  
  #### STARTING values (optional)
  if (is.null(start)){
    start <- list(beta = rep(0, model$d + 1))
  } else {
    start <- list(beta = start)
  }
  
  ## check starting values
  if(length(start$beta) != (model$d + 1) || !(is.numeric(c(start$beta)))){
    stop("invalid specification of starting values") 
  }
  
  ## fixed starting settings
  start$delta <- start$pdelta <- matrix(1, 1, model$d)
  start$psi   <- matrix(prior$V, 1, deff)  
  start$omega <- 1
  if(model$ri == 1){
    start$gamma <- start$pgamma <- matrix(1, 1, model$ri)   
    start$theta <- rnorm(model$ri, 0, 0.1)
    start$pi    <- matrix(1, 1, model$ri)
  } else {
    start$gamma <- start$pgamma <- start$theta <- start$pi <- NULL 
  }
  
  par.pois <- list(
    beta     = start$beta,
    delta     = start$delta,
    pdelta    = start$pdelta,
    omega     = start$omega, 
    psi       = start$psi,
    btilde    = btilde,
    theta     = start$theta,
    gamma     = start$gamma,
    pgamma    = start$pgamma,
    pi        = start$pi)
  
  
  ## mixture components of Gaussian mixture approximation
  mcomp <- mixcomp_poisson()
  
  c1  <- rep(1,n)%*%t(log(mcomp$w[1, ]) - 0.5*log(mcomp$v[1, ]))
  cm1 <- list(comp = list(m = mcomp$m[1, ], v = mcomp$v[1, ], w = mcomp$w[1, ]), 
              c1 = c1)
  
  
  #### MIXTURE COMPONENTS, given y
  compmix.pois <- get_mixcomp(y, mcomp) 
  
  
  #### SAVE --- matrices to save the MCMC draws
  beta <- matrix(0, mcmc$nmc, model$d + 1)
  colnames(beta) <- paste("beta", seq(0, model$d), sep = ".")
  
  if (model$d > 0){
    pdeltaBeta <- matrix(0, mcmc$nmc, model$d)
    colnames(pdeltaBeta) <- paste("pdeltaB", seq_len(model$d), sep = ".")
  } else {
    pdeltaBeta <- NULL
  }
  if (model$d > 0 && prior$slab=="Student"){
    psiBeta   <- matrix(0, mcmc$nmc, deff)
    colnames(psiBeta)   <- paste("psiB", seq_len(deff), sep = ".")  
  } else {
    psiBeta <- NULL
  }
  
  if (mcmc$msave && model$d > 0){
    omegaBeta <- matrix(0, mcmc$nmc, 1)
    colnames(omegaBeta) <- "omegaBeta"
    deltaBeta <- matrix(0, mcmc$nmc, model$d) 
    colnames(deltaBeta) <- paste("deltaB", seq_len(model$d), sep = ".")
  }
  if (model$ri == 1){
    pgammaBeta <- thetaBeta <- matrix(0, mcmc$nmc, model$ri)
    colnames(pgammaBeta) <- "pgammaB"
    colnames(thetaBeta)  <- "thetaB"
    bi <- matrix(0, mcmc$nmc, nC)
    colnames(bi) <- paste("b", seq_len(nC), sep = ".")
    if (mcmc$msave){
      piBeta <- gammaBeta <- matrix(0, mcmc$nmc, model$ri)
      colnames(piBeta) <- "piB"
      colnames(gammaBeta) <- "gammaB"
    } else piBeta <- gammaBeta <- NULL
  } else {
    pgammaBeta <- thetaBeta <- bi <- piBeta <- gammaBeta <- NULL
  }
  

  
  #### ------------------------------ MCMC ------------------------------
  for (imc in 1:mcmc$nmc){
    
    if (mcmc$verbose > 0){
      if (imc == 1) cat(paste0("\nMCMC for the ", txt_fun, " model", txt.verbose, ":\n\n"))
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
    
    
    ## time since burn-in phase
    if (imc == (mcmc$burnin + 1)) starttimeM <- proc.time()[3] 
    
    ## update 
    par.pois <- do.call(fun_select, args=list(y=y, X=X, offset=offset, H=H, 
                                              mcomp=mcomp, 
                                              compmix.pois = compmix.pois, 
                                              cm1=cm1, model=model,
                                              prior=prior, mcmc=mcmc,
                                              param=par.pois, imc=imc)
    )
    
    ## saving MCMC draws
    beta[imc,] <- par.pois$beta
    if (model$d > 0){
      pdeltaBeta[imc,] <- par.pois$pdelta
      if (prior$slab=="Student") psiBeta[imc,]  <- par.pois$psi
    }
    
    if (mcmc$msave && model$d > 0){
      omegaBeta[imc,] <- par.pois$omega
      deltaBeta[imc,] <- par.pois$delta
    }
    
    if(model$ri == 1){
      pgammaBeta[imc] <- par.pois$pgamma
      thetaBeta[imc]  <- par.pois$theta
      bi[imc,]         <- t(par.pois$btilde*par.pois$theta)
      
      if(mcmc$msave){
        piBeta[imc] <- par.pois$pi
        gammaBeta[imc] <- par.pois$gamma
      }
    }
  } # end(MCMC)
  
  finish <- proc.time()[3]
  durT <- finish - starttime
  durM <- finish - starttimeM
  dur  <- list(total = durT, durM = durM) 
  
  samplesP <- if(mcmc$msave){
    list(beta = beta, pdeltaBeta = pdeltaBeta, psiBeta = psiBeta, 
         omegaBeta = omegaBeta, deltaBeta = deltaBeta, 
         thetaBeta = thetaBeta, pgammaBeta = pgammaBeta, bi = bi,
         piBeta = piBeta, gammaBeta = gammaBeta)
  } else {
    list(beta = beta, pdeltaBeta = pdeltaBeta, psiBeta = psiBeta, 
         thetaBeta = thetaBeta, pgammaBeta = pgammaBeta, bi = bi)
  }
  
  ret <- list(samplesP = samplesP, data = list(y = y, offset = offset, X = X), 
              model.pois = model, mcmc = mcmc, prior.pois = prior, dur = dur, 
              BVS = BVS, start = start, family = model$family, call = cl,
              fun = fun_select) 
  class(ret) <- "pogit"
  return(ret)
}

