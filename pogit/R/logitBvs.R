#' Bayesian variable selection for the binomial logit model
#' 
#' This function performs Bayesian variable selection for binomial logit regression
#' models via spike and slab priors. A cluster-specific random intercept can be included 
#' in the model to account for within-cluster dependence with variance selection of 
#' the random intercept to determine whether there is between-cluster variation 
#' in the model. For posterior inference, a MCMC sampling algorithm is used which 
#' is based on data augmentation. 
#' 
#' The method provides Bayesian variable selection for binomial logit models 
#' using mixture priors with a spike and a slab component to identify regressors 
#' with a non-zero effect. More specifically, a Dirac spike is used, i.e. a 
#' point mass at zero and (by default), the slab component is specified as a scale 
#' mixture of normal distributions, resulting in a Student-t distribution with 
#' 2\code{psi.nu} degrees of freedom. 
#' In the more general random intercept model, variance selection of the random
#' intercept is based on the non-centered parameterization of the model, where
#' the signed standard deviation \eqn{\theta_\alpha} of the random intercept term 
#' appears as a further regression effect in the model equation. 
#' For details, see Wagner and Duller (2012). 
#'
#' The implementation of Bayesian variable selection further relies on the 
#' representation of the binomial logit model as a Gaussian regression model 
#' in auxiliary variables. Data augmentation is based on Fussl et 
#' al. (2013), who show that the binomial logit model can be represented as a 
#' linear regression model in the latent variable, which has an interpretation as 
#' the difference of aggregated utilities. The error distribution in the auxiliary 
#' model is approximated by a finite scale mixture of normal distributions, where
#' the mixture parameters are taken from the R package \code{binomlogit}. 
#' See Fussl (2014) for details. 
#' 
#' For details concerning the sampling algorithm see Dvorzak and Wagner (2016) 
#' and Wagner and Duller (2012). 
#' 
#' 
#' Details for the model specification (see arguments): 
#' \describe{ 
#'  \item{\code{model}}{\describe{\item{}{A list:}
#'    \item{\code{deltafix}}{an indicator vector of length \code{ncol(X)-1} 
#'    specifying which regression effects are subject to selection (i.e., 0 =
#'    subject to selection, 1 = fix in the model); defaults to a vector of zeros.} 
#'    \item{\code{gammafix}}{an indicator for variance selection of the random 
#'    intercept term (i.e., 0 = with variance selection (default), 1 = no
#'    variance selection); only used if a random intercept is includued in the 
#'    model (see \code{ri}).} 
#'    \item{\code{ri}}{logical. If \code{TRUE}, a cluster-specific
#'    random intercept is included in the model; defaults to \code{FALSE}.} 
#'    \item{\code{clusterID}}{a numeric vector of length equal to the number
#'    of observations containing the cluster ID c = 1,...,C for each observation 
#'    (required if \code{ri=TRUE}).} 
#'  }}
#'  \item{\code{prior}}{\describe{\item{}{A list:}
#'   \item{\code{slab}}{distribution of the slab component, i.e. "\code{Student}" 
#'    (default) or "\code{Normal}".}
#'    \item{\code{psi.nu}}{hyper-parameter of the Student-t slab (used for a 
#'    "\code{Student}" slab); defaults to 5.}
#'    \item{\code{m0}}{prior mean for the intercept parameter; defaults to 0.}
#'    \item{\code{M0}}{prior variance for the intercept parameter; defaults to 100.}
#'    \item{\code{aj0}}{a vector of prior means for the regression effects (which
#'    is encoded in a normal distribution, see notes); defaults to vector of zeros.}
#'    \item{\code{V}}{variance of the slab; defaults to 5.}
#'    \item{\code{w}}{hyper-parameters of the Beta-prior for the mixture weight
#'    \eqn{\omega}; defaults to \code{c(wa0=1, wb0=1)}, i.e. a uniform
#'    distribution.}
#'    \item{\code{pi}}{hyper-parameters of the Beta-prior for the mixture weight 
#'    \eqn{\pi}; defaults to \code{c(pa0=1, pb0=1)}, i.e. a uniform 
#'    distribution.}
#'  }}
#'  \item{\code{mcmc}}{\describe{\item{}{A list:}
#'    \item{\code{M}}{number of MCMC iterations after the burn-in phase; defaults 
#'    to 8000.} 
#'    \item{\code{burnin}}{number of MCMC iterations discarded as burn-in; 
#'    defaults to 2000.} 
#'    \item{\code{thin}}{thinning parameter; defaults to 1.} 
#'    \item{\code{startsel}}{number of MCMC iterations drawn from the unrestricted
#'    model (e.g., \code{burnin/2}); defaults to 1000.} 
#'    \item{\code{verbose}}{MCMC progress report in each \code{verbose}-th 
#'    iteration step; defaults to 500. If \code{verbose=0}, no output is 
#'    generated.} 
#'    \item{\code{msave}}{returns additional output with variable 
#'    selection details (i.e. posterior samples for \eqn{\omega}, 
#'    \eqn{\delta}, \eqn{\pi}, \eqn{\gamma}); defaults to \code{FALSE}.} 
#' }}}
#' 
#' 
#' @param y an integer vector of binomial counts
#' @param N an integer vector containing the number of trials
#' @param X a design matrix (including an intercept term)
#' @param model an (optional) list specifying the structure of the model 
#'  (see details)
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
#'  \item{\code{samplesL}}{a named list containing the samples from the posterior
#'    distribution of the parameters in the binomial logit model 
#'    (see also \code{msave}): 
#'    \describe{
#'    \item{\code{alpha, thetaAlpha}}{regression coefficients \eqn{\alpha} and
#'    \eqn{\theta_\alpha}}
#'    \item{\code{pdeltaAlpha}}{P(\eqn{\delta_\alpha}=1)}
#'    \item{\code{psiAlpha}}{scale parameter \eqn{\psi_\alpha} of the slab component}
#'    \item{\code{pgammaAlpha}}{P(\eqn{\gamma_\alpha}=1)}
#'    \item{\code{ai}}{cluster-specific random intercept}
#'    }}
#'  \item{\code{data}}{a list containing the data \code{y}, \code{N} and \code{X}}
#'  \item{\code{model.logit}}{a list containing details on the model specification, 
#'  see details for \code{model}}
#'  \item{\code{mcmc}}{see details for \code{mcmc}}
#'  \item{\code{prior.logit}}{see details for \code{prior}}
#'  \item{\code{dur}}{a list containing the total runtime (\code{total}) 
#'    and the runtime after burn-in (\code{durM}), in seconds}
#'  \item{\code{BVS}}{see arguments}
#'  \item{\code{start}}{a list containing starting values, see arguments}
#'  \item{\code{family}}{"logit"}
#'  \item{\code{call}}{function call}
#' 
#' @note If prior information on the regression parameters is available, this
#'  information is encoded in a normal distribution instead of the 
#'  spike and slab prior (\code{BVS} is set to \code{FALSE}). 
#'  
#'  For binary observations, a vector of ones for the number of trials \code{N}
#'  is required. 
#' 
#' @seealso \code{\link{pogitBvs}}
#' @keywords models
#'  
#' @references Dvorzak, M. and Wagner, H. (2016). Sparse Bayesian modelling
#'  of underreported count data. \emph{Statistical Modelling}, \strong{16}(1),
#'  24 - 46, \url{http://dx.doi.org/10.1177/1471082x15588398}.
#' @references Fussl, A., Fruehwirth-Schnatter, S. and Fruehwirth, R. (2013). 
#'   Efficient MCMC for Binomial Logit Models. \emph{ACM Transactions on 
#'   Modeling and Computer Simulation}, 23, \strong{1}, Article 3, 1-21.
#' @references Fussl, A. (2014). \code{binomlogit}: Efficient MCMC for Binomial 
#'   Logit Models. R package version 1.2, 
#'   \url{http://CRAN.R-project.org/package=binomlogit}.
#' @references Wagner, H. and Duller, C. (2012). Bayesian model selection for 
#'   logistic regression models with random intercept. \emph{Computational 
#'   Statistics and Data Analysis}, \strong{56}, 1256-1274.
#'   
#' @author Michaela Dvorzak <m.dvorzak@@gmx.at>, Helga Wagner
#' @import stats
#' @import utils
#' @export
#' 
#' @examples
#' \dontrun{
#' ## Examples below should take about 1-2 minutes.
#' 
#' # load simulated data set 'simul_binomial'
#' data(simul_binomial)
#' y <- simul_binomial$y
#' N <- simul_binomial$N
#' X <- as.matrix(simul_binomial[, -c(1, 2)])
#' 
#' # Bayesian variable selection for simulated data set
#' m1 <- logitBvs(y = y, N = N, X = X)
#' 
#' # print, summarize and plot results
#' print(m1)
#' summary(m1)
#' plot(m1)
#' 
#' # MCMC sampling without BVS with specific MCMC and prior settings
#' m2 <- logitBvs(y = y, N = N, X = X, prior = list(slab = "Normal"), 
#'                mcmc = list(M = 4000, burnin = 1000, thin = 5),
#'                BVS = FALSE)
#' print(m2)    
#' summary(m2)
#' plot(m2, maxPlots = 4) 
#' 
#' # BVS with specification of regression effects subject to selection
#' m3 <- logitBvs(y = y, N = N, X = X, mcmc = list(M = 4000, burnin = 1000), 
#'                model = list(deltafix = c(1, 1, 1, 0, 0, 0, 1, 0, 0)))   
#' print(m3)   
#' summary(m3)
#' plot(m3, burnin = FALSE, maxPlots = 4)
#' plot(m3, type = "acf", maxPlots = 4)       
#' }

logitBvs <- function(y, N, X, model = list(), prior = list(),  mcmc = list(), 
                     start = NULL, BVS = TRUE){
  
  starttime <- proc.time()[3]
  cl <- match.call()
    
  #### --- Check data and get dimensions
  if (is.null(y) || is.null(X)) stop("need 'y' and 'X' argument")
  if (!is.matrix(X)) stop("'X' must be a matrix")
  if (!is.vector(y)) stop("'y' must be a vector")
  if (!is.vector(N)) stop("'N' must be a vector")
  if (any(is.na(X))) stop("NA values in 'X' not allowed")
  if (any(is.na(y))) stop("NA values in 'y' not allowed")
  if (any(is.na(N))) stop("NA values in 'N' not allowed")
  if (any(y < 0)) stop("'y' must be positive")
  if (any(N < 0)) stop("'N' must be positive")
  if (length(y) != length(N)) stop("'y' and 'N' must have same length")
  if (length(y) != nrow(X)) stop("'y' and 'nrow(X)' must have same length")
    
  y <- as.integer(y)
  N <- as.integer(N) 
  if (!all(X[, 1] == 1)) X <- cbind(rep(1, dim(X)[1]), X)    
  if (any((N - y < 0))) stop("number of trials 'N' < number of binomial counts 'y'")  
    
  if (any(N==0)){   # only non-zero trials
    iN0 <- which(N==0)
    X <- X[-iN0, , drop = FALSE]
    y <- y[-iN0]
    N <- N[-iN0]
  }
  n <- length(y)
    
  d <- ncol(X) - 1
  W <- X    
  colnames(W) <- paste("W", seq(0,d), sep = "")
     
    
  #### --- Check and update parameter lists ---
    
  #### MODEL structure
  defaultModel <- list(
    deltafix  = matrix(0,1,d),    # effects subject to selection
    gammafix  = 0,                # random intercept selection
    ri        = FALSE,            # RI model for clustered observations
    clusterID = NULL              # cluster index c = 1,...,C
    )
  model <- modifyList(defaultModel, as.list(model))
  
  ## check parameters
  model$deltafix <- as.matrix(model$deltafix)
  model$gammafix <- as.matrix(model$gammafix) 
  model$ri       <- as.integer(model$ri)
  model$d        <- d 
  model$family   <- "logit"
  model$Zl       <- seq_len(n) 
  
  if (model$d > 0){
    if (length(model$deltafix) != model$d || !all(model$deltafix %in% c(0, 1))){
      stop("invalid specification of effects subject to selection")
    }
  } else {
    if (length(model$deltafix) != model$d){
      stop("invalid 'deltafix' argument")
    }
  }

  if (length(model$gammafix)>1 || !(model$gammafix %in% c(0, 1)) || !(model$ri %in% c(0, 1)))
      stop("invalid specification of random intercept selection")
      
  if (model$ri == 1){
      if (!is.vector(model$clusterID)){
        stop("you must specify a cluster ID for random intercept selection")
      }
      if (min(model$clusterID) != 1 || !(is.numeric(model$clusterID))){
        stop("specify the cluster ID as c = 1,...,C")
      }
      # random intercept model specification
      if (any(N == 0)) model$clusterID <- model$clusterID[-iN0]
      if (length(model$clusterID) != n){
        stop("'y' and 'clusterID' must have same length")
      }
      model$Zl <- model$clusterID
      nC <- max(model$Zl)
      atilde <- rnorm(nC, 0, 1)
      H <- matrix(0, n, nC)
      for (i in 1:nC){
        H[which(model$Zl == i), i] <- 1
      }
  } else {
      atilde <- H <- NULL      
  }
      
  df   <- d + 1 
  deff <- d + model$ri
  dall <- deff + 1
    
    
  #### PRIOR settings and hyper-parameters controlling the priors
  defaultPrior <- list(
    slab = "Student",         # slab distribution
    psi.nu = 5,               # hyper-parameter of slab variance
    m0 = 0,                   # only relevant if slab="Normal"  
    M0 = 100,                 # prior variance of intercept parameter
    aj0 = rep(0, deff),       # only relevant if slab="Normal"    
    V = 5,                    # slab variance
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
  if (model$ri==1) prior$invB0 <- diag(nC)
     
    
    
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
                                "'startsel' was set to ", mcmc$burnin/2, sep = "")))
    mcmc$startsel <- mcmc$burnin/2 
  }
  if (BVS && sum(sum(model$deltafix) + sum(model$gammafix)) == sum(model$d + model$ri)){
    warning(simpleWarning(paste("invalid 'BVS' argument:\n",
                                "'BVS' was set to ", FALSE, sep="")))
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
    start <- list(alpha = rep(0, model$d + 1))
  } else {
    start <- list(alpha = start)
  }
    
  ## check starting values
  if(length(start$alpha) != (model$d + 1) || !(is.numeric(c(start$alpha)))){
    stop("invalid specification of starting values") 
  }
    
  ## fixed starting settings 
  start$delta <- start$pdelta <- matrix(1, 1, model$d)
  start$psi   <- matrix(prior$V, 1, deff)  
  start$omega <- 1
  if(model$ri==1){
    start$gamma <- start$pgamma <- matrix(1, 1, model$ri)   
    start$theta <- rnorm(model$ri, 0, 0.1)
    start$pi    <- matrix(1, 1, model$ri)
  } else {
    start$gamma <- start$pgamma <- start$theta <- start$pi <- NULL 
  }

  
  par.logit <- list(
    alpha     = start$alpha,
    delta     = start$delta,
    pdelta    = start$pdelta,
    omega     = start$omega, 
    psi       = start$psi,
    atilde    = atilde,
    theta     = start$theta,
    gamma     = start$gamma,
    pgamma    = start$pgamma,
    pi        = start$pi)
    

    
  #### DATA AUGMENTATION (part I) for binomial logit models taken from "binomlogit"
  compmix.bin <- dataug_binom_dRUM1(y, N) 
    

  #### SAVE --- matrices to save the MCMC draws
  alpha <- matrix(0, mcmc$nmc, model$d + 1)
  colnames(alpha) <- paste("alpha", seq(0, model$d), sep = ".")
  if (model$d > 0){
    pdeltaAlpha <- matrix(0, mcmc$nmc, model$d)
    colnames(pdeltaAlpha) <- paste("pdeltaA", seq_len(model$d), sep = ".")
  } else {
    pdeltaAlpha <- NULL
  }
  if (model$d > 0 && prior$slab=="Student"){
    psiAlpha   <- matrix(0, mcmc$nmc, deff)
    colnames(psiAlpha)   <- paste("psiA", seq_len(deff), sep = ".")  
  } else {
    psiAlpha <- NULL
  }

  if (mcmc$msave && model$d > 0){
    omegaAlpha <- matrix(0, mcmc$nmc, 1)
    colnames(omegaAlpha) <- "omega"
    deltaAlpha <- matrix(0, mcmc$nmc, model$d) 
    colnames(deltaAlpha) <- paste("deltaA", seq_len(model$d), sep = ".")
  }
  if (model$ri == 1){
    pgammaAlpha <- thetaAlpha <- matrix(0, mcmc$nmc, model$ri)
    colnames(pgammaAlpha) <- "pgammaA"
    colnames(thetaAlpha)  <- "theta"
    ai <- matrix(0, mcmc$nmc, nC)
    colnames(ai) <- paste("a", seq_len(nC), sep = ".")
    if (mcmc$msave){
      piAlpha <- gammaAlpha <- matrix(0, mcmc$nmc, model$ri)
      colnames(piAlpha) <- "piA"
      colnames(gammaAlpha) <- "gammaA"
    } else piAlpha <- gammaAlpha <- NULL
  } else {
    pgammaAlpha <- thetaAlpha <- ai <- piAlpha <- gammaAlpha <- NULL
  }


    
  #### ------------------------------ MCMC ------------------------------
  for (imc in 1:mcmc$nmc){

    if (mcmc$verbose > 0){
      if (imc == 1) cat(paste("\nMCMC for the binomial logit model", txt.verbose, ":\n\n", sep = ""))
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
    
    ## update 
    par.logit <- select_logit(y, N, W, H, compmix.bin, model, prior, mcmc, par.logit, imc) 
      
      
    ## saving MCMC draws
    alpha[imc,] <- par.logit$alpha
    if (model$d > 0){
      pdeltaAlpha[imc,] <- par.logit$pdelta
      if (prior$slab=="Student") psiAlpha[imc,] <- par.logit$psi
    }
        
    if (mcmc$msave && model$d > 0){
      omegaAlpha[imc,] <- par.logit$omega
      deltaAlpha[imc,] <- par.logit$delta
    }
      
    if (model$ri == 1){
      pgammaAlpha[imc] <- par.logit$pgamma
      thetaAlpha[imc]  <- par.logit$theta
      ai[imc,]         <- t(par.logit$atilde*par.logit$theta)
        
      if (mcmc$msave){
        piAlpha[imc] <- par.logit$pi
        gammaAlpha[imc] <- par.logit$gamma
      }
    }
  } # end(MCMC)
      
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
  
  
  ret <- list(samplesL = samplesL, data = list(y = y, N = N, X = W), 
              model.logit = model, mcmc=mcmc, prior.logit = prior, 
              dur = dur, BVS = BVS, start = start, family = model$family, 
              call = cl) 
  class(ret) <- "pogit"
  return(ret)
}


    

