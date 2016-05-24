#' Bayesian variable selection for the negative binomial model
#' 
#' This function performs Bayesian regression modelling of overdispersed count 
#' data including variable selection via spike and slab priors. 
#' Posterior inference is based on MCMC sampling techniques.
#' 
#' The method provides Bayesian variable selection in regression modelling 
#' of overdispersed count data. The negative binomial 
#' distribution is derived marginally from a Poisson-Gamma (mixture) model, 
#' which can be interpreted as an overdispersed Poisson model with 
#' observation-specific random intercept log \eqn{\gamma}, where 
#' \eqn{\gamma|\rho \sim \Gamma(\rho,\rho)}. A hyper-prior for \eqn{\rho} is
#' specified as \eqn{\rho \sim \Gamma(c_0,C_0)}, see details for \code{prior} 
#' below. 
#' By default, variable selection is incorporated in the model based on mixture 
#' priors with a spike and a slab component for the regression effects 
#' \eqn{\beta}. More specifically, a Dirac spike is used, i.e. a point mass at 
#' zero, and (by default), the slab component is specified as a scale mixture 
#' of normal distributions, resulting in a Student-t distribution with 
#' 2\code{psi.nu} degrees of freedom. 
#' 
#' The MCMC sampling scheme relies on the representation of the conditional 
#' Poisson model as a Gaussian regression model in auxiliary variables, 
#' as described in \code{\link{poissonBvs}}. Data augmentation is based
#' on the auxiliary mixture sampling algorithm of Fruehwirth-Schnatter 
#' et al. (2009). For details concerning the algorithm, see 
#' Dvorzak and Wagner (2016b), available on request from the authors.   
#' 
#' Details for the model specification (see arguments): 
#' \describe{ 
#'  \item{\code{model}}{\describe{\item{}{A list:}  
#'    \item{\code{deltafix}}{an indicator vector of length 
#'    \code{ncol(X)-1} specifying which regression effects are subject to selection
#'    (i.e., 0 = subject to selection, 1 = fix in the model); defaults to a vector 
#'    of zeros.} 
#'  }}
#'  \item{\code{prior}}{\describe{\item{}{A list:}      
#'    \item{\code{slab}}{distribution of the slab component, i.e. "\code{Student}" 
#'    (default) or "\code{Normal}".} 
#'    \item{\code{psi.nu}}{hyper-parameter of the Student-t slab (used for a 
#'    "\code{Student}" slab); defaults to 5.} 
#'    \item{\code{m0}}{prior mean for the intercept parameter; defaults to 0.} 
#'    \item{\code{M0}}{prior variance for the intercept parameter; defaults to 100.} 
#'    \item{\code{aj0}}{a vector of prior means for the regression effects 
#'    (which is encoded in a normal distribution, see notes); defaults to vector of 
#'    zeros.} 
#'    \item{\code{V}}{variance of the slab; defaults to 5.} 
#'    \item{\code{w}}{hyper-parameters of the Beta-prior for the mixture weight
#'    \eqn{\omega}; defaults to \code{c(wa0=1, wb0=1)}, i.e. a uniform distribution.} 
#'    \item{\code{c0, C0}}{scale and rate of the gamma prior for the 
#'    hyper-parameter \eqn{\rho}; defaults to 2 and 1.} 
#'    \item{\code{eps}}{tuning parameter in the MH-step to sample \eqn{\rho}; 
#'    defaults to 0.05.} 
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
#'    (i.e. posterior samples for \eqn{\omega}, \eqn{\delta}); defaults to \code{FALSE}.}
#'  }}}
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
#'  \item{\code{samplesNB}}{a named list containing the samples from the posterior
#'    distribution of the parameters in the negative binomial model 
#'    (see also \code{msave}): 
#'    \describe{
#'    \item{\code{beta, rho}}{regression coefficients \eqn{\beta} and \eqn{\rho}}
#'    \item{\code{pdeltaBeta}}{P(\eqn{\delta_\beta}=1)}
#'    \item{\code{psiBeta}}{scale parameter \eqn{\psi_\beta} of the slab component}
#'    }}
#'  \item{\code{data}}{a list containing the data \code{y}, \code{offset} and \code{X}}
#'  \item{\code{model.nb}}{a list containing details on the model specification, 
#'  see details for \code{model}}
#'  \item{\code{mcmc}}{see details for \code{mcmc}}
#'  \item{\code{prior.nb}}{see details for \code{prior}}
#'  \item{\code{dur}}{a list containing the total runtime (\code{total}) 
#'    and the runtime after burn-in (\code{durM}), in seconds}
#'  \item{\code{acc.rho}}{acceptance rate of parameter \eqn{\rho}}
#'  \item{\code{BVS}}{see arguments}
#'  \item{\code{start}}{a list containing starting values, see arguments}
#'  \item{\code{family}}{"negbin"}
#'  \item{\code{call}}{function call}
#'  
#' @note Alternatively, a Poisson model with observation-specific normal 
#' random intercept (i.e., a Poisson-log-normal mixture model) can be used 
#' to deal with overdispersion of count data, which is provided in the 
#' function \code{\link{poissonBvs}}.
#' @note If prior information on the regression parameters is available, this 
#' information is encoded in a normal distribution instead of the spike 
#' and slab prior (consequently, \code{BVS} is set to \code{FALSE}).
#'
#' @seealso \code{\link{poissonBvs}}
#' @keywords models
#' 
#' @references Dvorzak, M. and Wagner, H. (2016b). Bayesian inference for 
#' overdispersed count data subject to underreporting - An application to
#' norovirus illness in Germany. (Unpublished) working paper. 
#' @references Fruehwirth-Schnatter, S., Fruehwirth, R., Held, L. and Rue, H. 
#'   (2009). Improved auxiliary mixture sampling for hierarchical models of 
#'   non-Gaussian data. \emph{Statistics and Computing}, \strong{19}, 479 - 492.
#'    
#' @author Michaela Dvorzak <m.dvorzak@@gmx.at>
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
#' data(simul_pois1)
#' y <- simul_pois1$y
#' X <- as.matrix(simul_pois1[, -1])
#' 
#' # Bayesian variable selection for simulated data set
#' m1 <- negbinBvs(y = y, X = X)
#' 
#' # print results (check acceptance rate for 'rho')
#' print(m1)
#' 
#' # re-run with adapted tuning parameter 'eps'
#' m2 <- negbinBvs(y = y, X = X, prior = list(eps = 0.4)) 
#' 
#' # print and summarize results
#' print(m2)
#' summary(m2)
#' 
#' # alternatively, compare results to overdispersed Poisson model with 
#' # normal random intercept (subject to selection), provided in 'poissonBvs' 
#' 
#' # specify observation-specific random intercept
#' cID <- seq_along(y)
#' m3  <- poissonBvs(y = y, X = X, model = list(ri = TRUE, clusterID = cID))
#' 
#' # print, summarize and plot results
#' print(m3)
#' summary(m3) 
#' # note that thetaB is not selected (!)
#' 
#' plot(m3, burnin = FALSE, thin = FALSE)
#' 
#' 
#' ## ------ (use data set "azdrg112" from package "COUNT") ------
#' 
#' if (!requireNamespace("COUNT", quietly = TRUE)){
#'  stop("package 'COUNT' is needed for this example to work. 
#'        Please install it.")
#' }
#' 
#' library(COUNT)
#' # load data set 'azdrg112' 
#' # (Arizona Medicare data for DRG (Diagnostic Related Group) 112)
#' data(azdrg112) 
#' 
#' y <- as.numeric(azdrg112$los) # hospital length of stay: 1-53 days
#' X <- as.matrix(azdrg112[,-1]) # covariates (gender, type1, age75)
#' m4 <- negbinBvs(y = y, X = X, mcmc = list(M = 4000))  
#' 
#' # print results (check acceptance rate for 'rho')
#' print(m4)
#' summary(m4)
#' plot(m4, burnin = FALSE)
#' 
#' # adapte tuning parameter eps (and set BVS to FALSE)
#' prior <- list(eps = 0.1)
#' m5 <- negbinBvs(y = y, X = X, mcmc = list(M = 4000), prior = prior, 
#'                 BVS = FALSE)  
#'
#' # print, summarize and plot results
#' print(m5)
#' summary(m5)
#' plot(m5, burnin = FALSE, thin = FALSE)
#' plot(m5, type = "acf", lag.max = 50)
#' }

negbinBvs <- function(y, offset = NULL, X, model = list(), prior  = list(),
                       mcmc = list(), start = NULL, BVS = TRUE){

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
  d <- ncol(X) - 1
  colnames(X) <- paste("X", seq(0, d), sep="")

  #### --- Check and update parameter lists

  #### MODEL structure
  defaultModel <- list(
    deltafix  = matrix(0,1,d)    # effects subject to selection
  )
  model <- modifyList(defaultModel, as.list(model))

  ## check parameters
  model$deltafix <- as.matrix(model$deltafix)
  model$d      <- d 
  model$family <- "negbin"

  if (model$d > 0){
    if (length(model$deltafix) != model$d || !all(model$deltafix %in% c(0, 1))){
      stop("invalid specification of effects subject to selection")
    }
  } else {
    if (length(model$deltafix) != model$d){
      stop("invalid 'deltafix' argument")
    }
  }

  df   <- d + 1 
  deff <- d 
  dall <- deff + 1


  #### PRIOR settings and hyper-parameters controlling the priors
  defaultPrior <- list(
    slab = "Student",         # slab distribution
    psi.nu = 5,               # hyper-parameter of slab variance
    V = 5,                    # slab variance
    M0 = 100,                 # prior variance of intercept parameter
    m0 = 0,                   # only relevant if slab="Normal"  
    aj0 = rep(0, d),          # only relevant if slab="Normal"    
    w = c(wa0=1,wb0=1),       # w ~ Beta(wa0,wb0)
    c0 = 2,                   # gamma shape of hyper-prior
    C0 = 1,                   # gamma rate of hyper-prior
    eps = 0.05                # tuning parameter for MH step (rho)
  )

  prior <- modifyList(defaultPrior, as.list(prior))
  if (!all(prior$aj0==0)){
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

  if (prior$slab=="Normal"){
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
                 slab %in% c("Student", "Normal"),
                 c0 > 0,
                 C0 > 0
                 )
       )

  if (prior$slab=="Student") with(prior, stopifnot(psi.nu > 0))
  if (length(prior$aj0) != deff)
    stop("invalid specification of prior means for regression effects")

  if (!(prior$eps > 0)){
    stop("invalid specification of tuning parameter 'eps'")
  }


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
                 typeof(msave)=="logical")
  )

  if (typeof(BVS) != "logical") stop("invalid 'BVS' argument")    
  if (BVS && mcmc$startsel > mcmc$nmc){
    warning(simpleWarning(paste("invalid 'BVS' or 'startsel' argument:\n",
                              "'startsel' was set to ", mcmc$burnin/2, sep = "")))
    mcmc$startsel <- mcmc$burnin/2 
  }
  if (BVS && sum(model$deltafix) == model$d){
    warning(simpleWarning(paste("invalid 'BVS' argument:\n",
                                "'BVS' was set to ", FALSE, sep = "")))
    BVS <- FALSE
  }
  if (!BVS){
    mcmc$startsel  <- mcmc$M + mcmc$burnin + 1
    model$deltafix <- matrix(1, 1, model$d)
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
  if (length(start$beta) != (model$d + 1) || !(is.numeric(c(start$beta)))){
    stop("invalid specification of starting values") 
  }

  ## fixed starting settings
  start$delta <- start$pdelta <- matrix(1, 1, model$d)
  start$psi   <- matrix(prior$V, 1, deff)  
  start$omega <- 1
  start$rho <- rgamma(1, prior$c0, prior$C0)
  start$g   <- rgamma(n, start$rho, start$rho)  


  par.nb <- list(
    beta      = start$beta,
    delta     = start$delta,
    pdelta    = start$pdelta,
    omega     = start$omega, 
    psi       = start$psi,
    g         = start$g,
    rho       = start$rho,
    rho.acc   = NULL)


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

  rho <- rho.acc <- matrix(0, mcmc$nmc, 1)
  colnames(rho) <- "rho"
  colnames(rho.acc) <- "rho.acc"   

  g <- matrix(0, mcmc$nmc, n)


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


  #### ------------------------------ MCMC ------------------------------
  for (imc in 1:mcmc$nmc){
  
    if (mcmc$verbose > 0){
      if (imc==1) cat(paste("\nMCMC for the negative binomial model", txt.verbose, ":\n\n", sep = ""))
      if (is.element(imc, c(1:5, 10, 20, 50, 100, 200, 500))){
        cat("it =", imc, "/--- duration of MCMC so far:", 
            round(timediff <- proc.time()[3] - starttime, 2), "sec.,  expected time to end:", 
            round((timediff/(imc - 1) * mcmc$nmc - timediff)/60, 2), " min. \n")
        flush.console()
      } else if (imc %% mcmc$verbose==0 && imc < mcmc$nmc){
        cat("it =", imc, "/--- duration of MCMC so far:", 
            round(timediff <- proc.time()[3] - starttime, 2), "sec.,  expected time to end:", 
            round((timediff/(imc - 1) * mcmc$nmc - timediff)/60, 2), " min. \n")
        flush.console()
      }
      else if (imc==mcmc$nmc) {
        cat("it =", imc, "/--- duration of MCMC (total):", 
            round(timediff <- proc.time()[3] - starttime, 2), "sec. \n \n")
        flush.console() 
      }
    }# end(verbose)
  
    ## time since burn-in phase
    if (imc == (mcmc$burnin + 1)) starttimeM <- proc.time()[3] 
  
    ## update 
    par.nb <- select_negbin(y, X, offset, mcomp, compmix.pois, cm1, model, prior, 
                            mcmc, par.nb, imc)
  
  
    ## saving MCMC draws
    beta[imc,] <- par.nb$beta
    if (model$d > 0){
      pdeltaBeta[imc,] <- par.nb$pdelta
      if (prior$slab=="Student") psiBeta[imc,]  <- par.nb$psi
    }
  
    rho[imc]     <- par.nb$rho
    rho.acc[imc] <- par.nb$rho.acc
    g[imc,]   <- par.nb$g
  
    if (mcmc$msave && model$d > 0){
      omegaBeta[imc,] <- par.nb$omega
      deltaBeta[imc,] <- par.nb$delta
    }
  } # end(MCMC)

  finish <- proc.time()[3]
  durT <- finish - starttime
  durM <- finish - starttimeM
  dur  <- list(total = durT, durM = durM) 
  
  # acceptance rate for rho
  acc.rate <- round(mean(rho.acc[-c(1:mcmc$burnin)])*100, 2)

  samplesNB <- if(mcmc$msave){
    list(beta = beta, pdeltaBeta = pdeltaBeta, psiBeta = psiBeta, 
         omegaBeta = omegaBeta, deltaBeta = deltaBeta, 
         rho = rho, g = g)
  } else {
    list(beta = beta, pdeltaBeta = pdeltaBeta, psiBeta = psiBeta, 
         rho = rho, g = g)
  }

  ret <- list(samplesNB = samplesNB, data = list(y = y, offset = offset, X = X), 
            model.nb = model, mcmc = mcmc, prior.nb = prior, dur = dur,
            acc.rho = acc.rate, BVS = BVS, start = start, 
            family = model$family, call = cl) 
  class(ret) <- "pogit"
  return(ret)
}


