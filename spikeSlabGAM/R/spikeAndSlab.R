#' @include utils.R
{}

#' Set up and sample a spike-and-slab prior model.
#'
#' This function sets up a spike-and-slab model for variable selection and model
#' choice in generalized additive models and samples its posterior. It uses a
#' blockwise Metropolis-within-Gibbs sampler and the redundant multiplicative
#' parameter expansion described in the reference. This routine is not meant to
#' be called directly by the user -- \code{\link{spikeSlabGAM}} provides a
#' formula-based interface for specifying models and takes care of (most of) the
#' housekeeping. Sampling of the chains is done in parallel using package
#' \code{parallel}. A "SOCK" cluster is set up under Windows to do so (and
#' closed after computations are done, I try to clean up after myself), see
#' \code{\link[parallel]{makeCluster}} etc. Use \code{options(mc.cores =<foo>)}
#' to set the (maximal) number of processes forked by the parallelization. If
#' \code{options()$mc.cores} is unspecified, it is set to 2.
#'
#' Details for model specification: \describe{
#' \item{\code{hyperparameters}}{\describe{\item{}{A list:}
#' \item{\code{w}}{hyperparameters for the \eqn{Beta}-prior for \eqn{w};
#' defaults to \code{c(alphaW = 1, betaW = 1)}, i.e. a uniform distribution.}
#' \item{\code{tau2}}{hyperparameters for the \eqn{\Gamma^{-1}}-prior of the
#' hypervariances \eqn{\tau^2}; defaults to \code{c(a1 = 5, a2 = 25)}}
#' \item{\code{gamma}}{sets \eqn{v_0}, the ratio between the spike and slab
#' variances, defaults to \code{c(v0 = 0.00025)}}
#' \item{\code{sigma2}}{hyperparameters for \eqn{\Gamma^{-1}}-prior for error
#' variance; defaults to \code{c(b1 = 1e-4, b2 = 1e-4)}. Only relevant for Gaussian
#' response.} \item{\code{varKsi}}{variance for prior of \eqn{\xi}, defaults to
#' 1} \item{\code{ksiDF}}{defaults to 0 for a gaussian prior for \eqn{\xi}, else
#' induces a t-prior for \eqn{\xi}} with \code{ksiDF} degrees of freedom.}}
#' \item{\code{model}}{\describe{\item{}{A list:}
#' \item{\code{groupIndicators}}{a factor that maps the columns of X to the
#' different model terms} \item{\code{H}}{a matrix containing the hierarchy of
#' the penalized model terms} \item{\code{n}}{number of observations}
#' \item{\code{q}}{length of \eqn{\beta}} \item{\code{scale}}{scale/weights of
#' the response, defaults to \code{rep(1, n)}, use this to specify number of
#' trials for binomial data} \item{\code{offset}}{defaults to \code{rep(0,
#' n)}}}} \item{\code{mcmc}}{\describe{\item{}{A list:}
#' \item{\code{nChains}}{how many parallel chains to run: defaults to 3}
#' \item{\code{chainLength}}{how many samples should be generated per chain,
#' defaults to 500} \item{\code{burnin}}{how many initial iterations should be
#' discarded, defaults to 100} \item{\code{thin}}{save only every \code{thin}-th
#' iteration, defaults to 5} \item{\code{verbose}}{verbose output and  report
#' progress? defaults to TRUE} \item{\code{returnSamples}}{defaults to TRUE}
#' \item{\code{sampleY}}{generate samples of y and its conditional expectation
#' from posterior predictive? defaults to FALSE}
#' \item{\code{useRandomStart}}{use random draw or ridge estimate for beta as
#' starting value? defaults to TRUE, i.e. random starting values.}
#' \item{\code{blocksize}}{approx. blocksizes of the updates for \eqn{\alpha,
#' \xi}. Defaults to 50 for gaussian responses and 5/15 for non-gaussian
#' responses.} \item{\code{scalemode}}{how to do term-wise rescaling of
#' subvectors of \eqn{\xi} in each iteration:  0 means no rescaling, 1 means
#' rescaling s.t. each mean\eqn{(|\xi_g|) = 1},  2 means rescaling s.t. each
#' max\eqn{(|\xi_g|) = 1}} \item{\code{modeSwitching}}{probability to do P-IWLS
#' with the mode of the proposal set to the current value, which is useful if
#' the chain gets stuck. Defaults to \eqn{0.05}. Increase if accceptance rates
#' are too low.} \item{\code{reduceRet}}{don't return data and samples for
#' \eqn{\alpha, \xi, \tau^2}? defaults to FALSE}}}
#' \item{\code{start}}{\describe{\item{}{A list:} \item{\code{beta}}{starting
#' values for \eqn{\beta}. Defaults to a modified approximate ridge-penalized ML
#' estimate. See vignette for details on default specification.}
#' \item{\code{gamma}}{starting values for \eqn{\gamma}. Defaults to a vector of
#' 1's if \code{mcmc$useRandomStart} is \code{FALSE}, otherwise drawn from the
#' prior.} \item{\code{tau2}}{starting values for \eqn{\tau^2}. Defaults to the
#' mode of the prior if \code{mcmc$useRandomStart} is \code{FALSE}, otherwise
#' drawn from the prior.} \item{\code{sigma2}}{starting values for
#' \eqn{\sigma^2}. Only relevant for Gaussian response. Defaults to the variance
#' of the response divided by the number of covariates if
#' \code{mcmc$useRandomStart} is \code{FALSE}, otherwise drawn from the prior.}
#' \item{\code{w}}{starting value for \eqn{w}. Defaults to the mean of the prior
#' if \code{mcmc$useRandomStart} is \code{FALSE}, otherwise drawn from the
#' prior.} \item{\code{seed}}{Sets RNG seed for reproducible results. Parallel
#' chains are seeded with this seed incremented by the number of the chain.}}}}
#'
#' @param y response
#' @param X design matrix
#' @param family (character) the family of the response, defaults to
#'   normal/Gaussian response
#' @param hyperparameters a list of hyperparameters controlling the priors (see
#'   details)
#' @param model a list with information about the grouping structure of the
#'   model (see details)
#' @param mcmc (optional) list setting arguments for the sampler (see details)
#' @param start (optional) list containing the starting values for \eqn{\beta,
#'   \gamma, \tau^2, \sigma^2, w} and, optionally, the random seed
#' @return a list with components:\describe{ \item{\code{formula}}{see
#'   arguments} \item{\code{data}}{see arguments} \item{\code{family}}{see
#'   arguments} \item{\code{y}}{see arguments} \item{\code{X}}{see arguments}
#'   \item{\code{hyperparameters}}{see arguments} \item{\code{model}}{see
#'   arguments} \item{\code{mcmc}}{see arguments} \item{\code{start}}{see
#'   arguments} \item{\code{posteriorPred}}{a list with entries \code{mu} and
#'   \code{y} containing samples of the expected values and realizations of the
#'   response from the posterior predictive} \item{\code{postMeans}}{a list
#'   containing the posterior means of the parameters: \describe{
#'   \item{\code{beta}}{the regression coefficients} \item{\code{alpha}}{}
#'   \item{\code{ksi}}{} \item{\code{tau}}{hypervariances of the penalized model
#'   terms} \item{\code{gamma}}{inclusion indicator variables of the model
#'   terms} \item{\code{pV1}}{\eqn{P(\gamma = 1)}} \item{\code{w}}{hyperparameter
#'   for \code{gamma}} \item{\code{sigma2}}{error variance (for Gaussian data)}
#'   \item{\code{logLik}}{log likelihood} \item{\code{logPost}}{log of
#'   (unnormalized) posterior}}} \item{\code{samples}}{a list containing the
#'   posterior samples of the parameters, see above for explanation of the
#'   entries} \item{\code{DIC}}{a vector with \eqn{DIC, pD, \bar{D},\hat{D}}.
#'   Usually doesn't make much sense for this kind of model because of the
#'   posterior's multimodality.} \item{\code{fitted}}{a matrix with the
#'   posterior mean of the linear predictor in the first column and the
#'   posterior mean of the expected response in the second.}
#'   \item{\code{runTime}}{of the sampler, in seconds}}
#' @author Fabian Scheipl, Daniel Sabanes Bove
#' @references Scheipl, F. (2010) Normal-Mixture-of-Inverse-Gamma Priors for
#'   Bayesian Regularization and Model Selection in Structured Additive
#'   Regression Models. \emph{LMU Munich, Department of Statistics}: Technical
#'   Reports, No.84 (\url{http://epub.ub.uni-muenchen.de/11785/})
#' @export
#' @import parallel
#' @importFrom MASS ginv
#' @importFrom MCMCpack dinvgamma
#' @importFrom MCMCpack rinvgamma
#' @importFrom coda mcmc
#' @importFrom coda mcmc.list
#' @importFrom mvtnorm rmvnorm
#' @useDynLib spikeSlabGAM sampler
spikeAndSlab <- function(
  y,    # response (n x 1)
  X,    # design matrix with covariates (n x q)
  family = c("gaussian","binomial", "poisson"),
  hyperparameters = list(),    # prior hyperparameters
  model = list(),       # model structure
  mcmc = list(),        # MCMC sampler options
  start = list()        # start values for the sampler
)
{

  startTime <- Sys.time()

  if(is.null(X) || is.null(y)) stop("Need X and y.")

  familystr <- match.arg(family)
  family <- as.integer(switch(familystr,
    "gaussian" = 0,
    "binomial" = 1,
    "poisson" = 2))
  if((family!= 0)&&(any(y<0)))
    stop("non-gaussian responses must be non-negative.\n")
  if((family == 1)&&(any(y>1)))
    stop("binomial responses must be between 0 and 1.\n")
  if((family == 2)&&(any(y%%1 != 0)))
    stop("poisson reponses must be integers.")

  ### Check data and get dimensions
  y <- as.matrix(y)
  X <- as.matrix(X)

  n <- nrow(y)
  if(nrow(X) != n)
    stop(simpleError("X doesn't have same length as y."))

  q <- ncol(X)

  #### Complete and check parameter lists
  #### MCMC
  defaultMcmc <-        # Gibbs sampling options:
    list(
      nChains = 3,			#number of chains to run
      chainLength = 500,     # length of the output markov chain
      burnin = 100,    # how many starting samples are
      # discarded?
      thin = 5,   # only every thin'th sample is retained for the output chain
      verbose = TRUE,    # should information on progress etc been given?
      returnSamples = TRUE,      # return the samples as a list of mcmc objects?
      sampleY = FALSE,		   # sample Y from posterior predictive?
      useRandomStart = TRUE,    # use ridge estimate for beta as starting value or use random draw?
      blocksize = if(family == 0) {
        50
      } else c(5, 15),  # blocksizes  for QR-updates
      scalemode = 1, #rescale ksi in each iteration s.t. 0: no rescaling; 1: mean(|ksi_g|) = 1;  2: max(|ksi_g|) = 1
      allKsi = TRUE,  #use redundant parameterization for grps with d = 1?
      modeSwitching = 0.05,
      reduceRet = FALSE # return smaller object (without data, X, alpha/ksi/gamma samples)
    )
  mcmc <- modifyList(defaultMcmc, mcmc)
  if(length(mcmc$blocksize)== 2) {
    mcmc$blocksizeAlpha <- mcmc$blocksize[1]
    mcmc$blocksizeKsi <- mcmc$blocksize[2]
  } else {
    mcmc$blocksizeAlpha <- mcmc$blocksizeKsi <- mcmc$blocksize
  }


  ## check parameters
  mcmc$chainLength <- as.integer(mcmc$chainLength)
  mcmc$burnin <- as.integer(mcmc$burnin)
  mcmc$thin <- as.integer(mcmc$thin)


  with(mcmc,
    stopifnot(all(c(chainLength, burnin, thin) >= c(1, 0, 1)),
      scalemode %in% c(-1, 0, 1, 2)))


  mcmc$totalLength <- with(mcmc, burnin + thin * chainLength)


  ### MODEL

  H <- diag(q-1)
  colnames(H) <- rownames(H) <-  seq_len(q - 1)
  defaultModel <-
    list(
      groupIndicators = c("u", seq_len(q - 1)), #intercept, single coefficients
      H = H,
      n = n,
      q = q,
      scale = rep(1, n),
      offset = rep(0, n)
    )
  model <- modifyList(defaultModel, as.list(model))



  ## adapt length of groupIndicators
  if(length(model$groupIndicators) != model$q)
  {
    oldGroupIndicators <- model$groupIndicators
    model$groupIndicators <- rep(model$groupIndicators, length.out = q)

    warning(simpleWarning(paste("Length of model$groupIndicators did not equal",
      "the number of columns in design matrix X.\n",
      "Expanded to\n",
      paste(model$groupIndicators, collapse =" "),
      "\nfrom \n",
      paste(oldGroupIndicators, collapse =" "))))
  }
  model$groupIndicatorsOrig <- model$groupIndicators

  # expand model$groupIndicators s.t. each unpen param is in a separate block
  unpen <- which(model$groupIndicators =="u")
  pUnpen <- length(unpen)
  unpenLevels <- if(pUnpen) paste("u.", 1:length(unpen), sep ="") else c()
  model$groupIndicators[unpen] <-  unpenLevels

  #re-order model$groupIndicators and X  s.t. unpen grps and grps with
  #d = 1 (ksi = 1) come first.
  d1Levels <- if(pUnpen) {
    model$groupIndicators[which(!(model$groupIndicators %in%
        c("u", names(table(model$groupIndicators[-unpen]))[
          table(model$groupIndicators[-unpen])>1])))]
  } else {
    model$groupIndicators[which(!(model$groupIndicators %in%
        c(names(table(model$groupIndicators))[table(model$groupIndicators)>1])))]
  }
  d1Levels <- setdiff(d1Levels, unpenLevels)
  grpLevels <- unique(model$groupIndicators[which(!(model$groupIndicators %in%
      c(unpenLevels, d1Levels)))])

  model$groupIndicators <- factor(model$groupIndicators,
    levels = c(unpenLevels, d1Levels, grpLevels), ordered = T)

  properOrder <- order(model$groupIndicators)
  properPenOrder <- order(model$groupIndicators[!(model$groupIndicators %in%
      unpenLevels)])
  reverseOrder <- order(properOrder)
  reversePenOrder <- order(properPenOrder)


  model$groupIndicators <- model$groupIndicators[properOrder]

  ## do indicator hierarchy
  model$H <- model$H[c(d1Levels, grpLevels), c(d1Levels, grpLevels)]
  model$hierarchy <- apply(model$H * upper.tri(model$H), 1, function(x) {
    ind <- max(0, which(x != 0)-1)
  })
  model$properOrder <- properOrder
  model$reverseOrder <- reverseOrder

  model$d <- table(model$groupIndicators)
  model$nGroups <- length(model$d) #p length(alpha)
  model$G <- if(nlevels(model$groupIndicators)>1) {
    model.matrix( ~ 0 + model$groupIndicators, contrasts.arg = contr.treatment)
  } else{
    matrix(1, nrow = length(model$groupIndicators), ncol = 1)
  }
  model$nUnpenGrps <- pUnpen
  model$nPenGrps <- model$nGroups - pUnpen


  if(!mcmc$allKsi) {
    model$updateKsiGrps <-	which(!(unique(model$groupIndicators) %in%
        names(model$d)[model$d == 1]))
    model$updateKsi <- which(!(model$groupIndicators %in%
        names(model$d)[model$d == 1]))
  } else {
    #use redundant params for better mixing
    model$updateKsiGrps <- 1:model$nGroups
    model$updateKsi <- 1:q
  }


  model$qKsiUpdate <- length(model$updateKsi)
  model$qKsiNoUpdate <- model$q - model$qKsiUpdate

  #re-order X
  X <- X[, model$properOrder]

  ###  HYPER
  defaultHyperparameters <-
    list(
      w = c(alphaW = 1, betaW = 1),   # w ~ Beta(alphaW, betaW)
      tau2 = c(a1 = 5, a2 = 25), 		# tau^{2} ~ IG(a1, b1)
      gamma = c(v0 = 0.00025), # ratio spike/slab variance
      # variance
      sigma2 = c(b1 = 1e-4, b2 = 1e-4), # sigma^2 ~ IG(b1, b2)
      varKsi = rep(1, model$qKsiUpdate), #only relevant if ksiDF>0
      ksiDF = 0  #ksiDF == 0 --> ksi ~ N(+/-1, varKsi), else ksi ~ N(+/-1, 1/Gamma(df/2, df/2))
    )
  hyperparameters <- modifyList(defaultHyperparameters, as.list(hyperparameters))
  if(length(hyperparameters$varKsi)== 1) {
    hyperparameters$varKsi <- rep(hyperparameters$varKsi, model$qKsiUpdate)
  }

  ## check the parameters
  with(hyperparameters,
    stopifnot(all(w > 0),
      all(tau2 > 0),
      gamma >= 0,
      sigma2 > 0),
    ksiDF >= 0,
    all(varKsi>0),
    length(varKsi)== model$qKsiUpdate)

  ## insert the 1 into gamma hyperparameter
  hyperparameters$gamma["v1"] <- 1


  ##name hyperparams as expected
  names(hyperparameters$w) <- c('alpha','beta')
  names(hyperparameters$tau2) <- c('a1','a2')
  names(hyperparameters$gamma) <- c('v0','v1')
  names(hyperparameters$sigma2) <- c('b1','b2')




  #### START
  ## set seed if there is none
  if(!length(start$seed)) start$seed <- as.integer(round(1e6 * runif(1)))
  set.seed(start$seed)


  gamma <- matrix(if(mcmc$useRandomStart) {
    replicate(mcmc$nChains, as.vector(sample(hyperparameters$gamma,
      size = model$nPenGrps,
      replace = TRUE,
      prob = c(1, 0) + c(-1, 1) *
        rbeta(1, hyperparameters$w[1],
          hyperparameters$w[2]))))

  }else{
    replicate(mcmc$nChains, rep(1, model$nPenGrps))
  }, ncol = mcmc$nChains)
  tau2 <- matrix(if(mcmc$useRandomStart) {
    replicate(mcmc$nChains, as.vector(1/rgamma(model$nPenGrps,
      hyperparameters$tau2[1], hyperparameters$tau2[2])))
  }else{
    replicate(mcmc$nChains, rep(hyperparameters$tau2[2]/hyperparameters$tau2[1],
      model$nPenGrps))
  }, ncol = mcmc$nChains)
  sigma2 <- if(mcmc$useRandomStart) {
    replicate(mcmc$nChains, as.vector(1/rgamma(1,
      model$n/2 + hyperparameters$sigma2[1],
      var(y)/2 + hyperparameters$sigma2[2])))
  }else{
    replicate(mcmc$nChains, var(y)/model$nPenGrps)
  }
  w <- if(mcmc$useRandomStart) {
    replicate(mcmc$nChains, as.vector(rbeta(1, hyperparameters$w[1],
      hyperparameters$w[2])))
  }else{
    replicate(mcmc$nChains, hyperparameters$w[1]/(hyperparameters$w[1] +
        hyperparameters$w[2]))
  }

  if(is.null(start$beta)|| (family > 0)) {
    betaM <- switch(as.character(family),
      "0" = as.vector(solve(crossprod(sqrt(model$scale)* X) + diag(q),
        crossprod(X, y))),
      "1" = iwls.start(X, y, 1, model$scale, model$offset),
      "2" = iwls.start(X, y, 2, model$scale, model$offset))

    beta <- matrix(if (mcmc$useRandomStart) {
      tmp <- replicate(mcmc$nChains, betaM + rnorm(q))
      gammaLong <- model$G %*% rbind(matrix(1, ncol = ncol(gamma),
        nrow = model$nUnpenGrps), gamma)
      sapply(1:mcmc$nChains, function(i) tmp[, i] * gammaLong[,
        i])
    }
      else {
        replicate(mcmc$nChains, betaM)
      }, ncol = mcmc$nChains)
  } else beta <- NULL




  defaultStart <-  list(   # start values for the Gibbs sampler:
    beta = beta,
    gamma = gamma,
    tau2 = tau2,
    sigma2 = sigma2,
    w = w,
    seed = as.integer(start$seed) # optional RNG seed (integer or .Random.seed)
  )
  ## allow start close to origin (start in exact null not possible because divide-by-zero)
  if(!is.null(start$beta)) {
    if(sum(abs(start$beta))== 0) start$beta <- matrix(runif(q * mcmc$nChains,
      -.01, .01), nrow = q)
  }

  start <- modifyList(defaultStart, start)
  ## check starting values
  with(start,
    stopifnot(nrow(beta) == model$q,
      nrow(gamma) == model$nPenGrps,
      nrow(tau2) == model$nPenGrps,
      all(gamma %in% hyperparameters$gamma),
      all(tau2 > 0),
      all(is.finite(tau2)),
      is.finite(sigma2),
      sigma2 > 0,
      all(is.integer(seed))))

  # initialize&rescale ksi, alpha
  if(length(model$updateKsiGrps)) {
    ksi <- matrix(start$beta, ncol = mcmc$nChains)
    alpha <- matrix(1, nrow = model$nGroups, ncol = mcmc$nChains)
    alpha[-model$updateKsiGrps, ] <- drop(t(t(ksi) %*%
        model$G[, -model$updateKsiGrps]))
    rescales <-  apply( abs(ksi[model$updateKsi, , drop = F]), 2,  function(x) {
      tapply(x, model$groupIndicators[model$updateKsi, drop = T], mean)
    })
    alpha[model$updateKsiGrps, ] <- alpha[model$updateKsiGrps, ]* rescales
    ksi[-model$updateKsi, ] <- 1
    ksi[model$updateKsi, ] <- if(length(model$updateKsiGrps)>1) {
      ksi[model$updateKsi, , drop = F]/apply(rescales, 2, rep,
        model$d[model$updateKsiGrps])
    }else{
      ksi[model$updateKsi, , drop = F]/rescales
    }
    start$ksi <- ksi
    start$alpha <- alpha
  } else {
    start$ksi <- matrix(1, nrow = nrow(start$beta), ncol = mcmc$nChains)
    start$alpha <- start$beta
  }
  with(start,
    stopifnot(nrow(alpha) == model$nGroups,
      nrow(ksi) == model$q))


  #
  #
  ##### prepare info for blockwise updates
  #alpha-updates (zero-based indices for C)
  blocksAlpha <- max(1, round(model$nGroups/mcmc$blocksizeAlpha))
  alphaIndABegin <- alphaIndAEnd <- numeric(blocksAlpha)
  startInd <- 0
  for(i in 1:blocksAlpha) {
    alphaIndABegin[i]<- startInd
    alphaIndAEnd[i] <- ifelse(
      i == blocksAlpha,
      model$nGroups-1,
      min(model$nGroups-1,(startInd + mcmc$blocksizeAlpha - 1)))
    startInd <- alphaIndAEnd[i] + 1
  }

  #ksi-updates: define s.t. update blks don't split pen grps if possible
  blocksizes <- model$d[model$updateKsiGrps]
  tooLargeBlock <- blocksizes > 1.5 * mcmc$blocksizeKsi
  if(any(tooLargeBlock)) {
    blocksizes <- as.list(blocksizes)
    for(b in seq(along = blocksizes)) {
      if(tooLargeBlock[b]) {
        blcks <- round(blocksizes[[b]]/mcmc$blocksizeKsi)
        blcksSz <- blocksizes[[b]]%/%blcks
        if(sum(rep(blcksSz, blcks)) == blocksizes[[b]]) {
          blocksizes[[b]] <- rep(blcksSz, blcks)
        } else {
          blocksizes[[b]] <- c(rep(blcksSz, blcks-1), blcksSz +
              blocksizes[[b]]%%blcks)
        }
      }
    }
    blocksizes <- unlist(blocksizes)
  }
  ksiIndABegin <- ksiIndAEnd <- numeric()
  ksiIndABegin[1] <- ksiIndAEnd[1] <- 1
  block <- 0
  i <- 1
  while(block < length(blocksizes)) {
    if(block == max(which(cumsum(blocksizes) - ksiIndABegin[i] <=
        mcmc$blocksizeKsi))) {
      block <- block + 1
    }
    else block <- max(which(cumsum(blocksizes) - ksiIndABegin[i] <=
        mcmc$blocksizeKsi))
    ksiIndAEnd[i] <- cumsum(blocksizes)[block]
    if(block <  length(blocksizes)) {
      ksiIndABegin[i + 1] <- ksiIndAEnd[i] + 1
      i <- i + 1
    }
  }
  #zero-based indices for C:
  ksiIndAEnd <- ksiIndAEnd - 1
  ksiIndABegin <- ksiIndABegin - 1
  blocksKsi <- length(ksiIndAEnd)

  if(mcmc$verbose) {
    cat("\nModel has ", model$q, " coefficients in ",
      model$nGroups, "model terms.\n")
    cat("Blockwise sampling: alpha:", blocksAlpha, " block(s); xi:",
      blocksKsi, "block(s).\n")
  }


  ##

  pcts <- round(quantile(mcmc$burnin:mcmc$totalLength, seq(.1,.9, by =.1)))

  betaMat <- ksiMat <- matrix(0, nrow = mcmc$chainLength, ncol = model$q)
  tau2Mat <- gammaMat <- probV1Mat <- matrix(0, nrow = mcmc$chainLength,
    ncol = model$nPenGrps)
  alphaMat <- matrix(0, nrow = mcmc$chainLength, ncol = model$nGroups)
  wMat <- likMat <- logPostMat <- sigma2Mat <- matrix(0,
    nrow = mcmc$chainLength, ncol = 1)


  parallel <- if(.Platform$OS.type != "windows") {
    "parallel"
  } else {
    "snow"
  }


  if((is.null(options()$mc.cores) || is.na(options()$mc.cores))) {
    if(interactive()) {
      options(mc.cores = as.integer(readline(prompt =
          paste0("Setting up parallel computation:\n",
        "How many processes do you want to run? "))))
      stopifnot(!is.null(options()$mc.cores), !is.na(options()$mc.cores),
        options()$mc.cores>0)
    } else {
      message("Using 2 parallel processes.\n",
        "Use 'options(mc.cores = <YourNumberHere>)' to override next time.")
      options(mc.cores = 2)
    }
  }



  if(mcmc$verbose) {
    cat('\nstarting chain(s):\n')
    if(parallel =="parallel") {
      cat(paste(paste(rep("b", mcmc$nChains), collapse =""),
        "0",
        paste(rep("-", mcmc$nChains * 10-2), collapse =""),
        "100%\n", sep ="", collapse =""))
    }
    if(parallel == "snow") {
      cat("using <parallel> in SOCKET mode -- no progress info available.\n")
    }
  }

  do1Chain <- function(i) {
    # set seed for each chain s.t. results are reproducible
    set.seed(as.integer(start$seed + i))
    res <- .C("sampler",
      a1 = as.double(hyperparameters$tau2['a1']),
      a2 = as.double(hyperparameters$tau2['a2']),
      b1 = as.double(hyperparameters$sigma2['b1']),
      b2 = as.double(hyperparameters$sigma2['b2']),
      alphaW = as.double(hyperparameters$w['alpha']),
      betaW = as.double(hyperparameters$w['beta']),
      v0 = as.double(hyperparameters$gamma['v0']),
      varKsi = as.double(hyperparameters$varKsi),

      q =	as.integer(model$q),	qKsiUpdate = as.integer(model$qKsiUpdate),
      p =	as.integer(model$nGroups), pPen = as.integer(model$nPenGrps),
      n = as.integer(model$n),
      d = as.integer(model$d),


      beta = as.double(start$beta[, i]), alpha = as.double(start$alpha[, i]),
      ksi = as.double(start$ksi[, i]), tau2 = as.double(start$tau2[, i]),
      gamma = as.double(start$gamma[, i]), sigma2 = as.double(start$sigma2[i]),
      w = as.double(start$w[i]),

      y = as.double(y),
      X = as.double(X),
      G = as.double(model$G),
      scale = as.double(model$scale),
      offset = as.double(model$offset),

      blocksAlpha	= as.integer(blocksAlpha),
      indA1Alpha = as.integer(alphaIndABegin),
      indA2Alpha = as.integer(alphaIndAEnd),

      blocksKsi = as.integer(blocksKsi),
      indA1Ksi = as.integer(ksiIndABegin),
      indA2Ksi = as.integer(ksiIndAEnd),

      pcts = as.integer(pcts),
      burnin = as.integer(mcmc$burnin),
      thin = as.integer(mcmc$thin),
      totalLength = as.integer(mcmc$totalLength),
      verbose = as.integer(mcmc$verbose),
      ksiDF = as.double(hyperparameters$ksiDF),
      scaleMode = as.integer(mcmc$scalemode),
      modeSwitching = as.double(mcmc$modeSwitching),
      family = as.integer(family),
      acceptKsi = as.double(rep(0, blocksKsi)),
      acceptAlpha = as.double(rep(0, blocksAlpha)),

      betaMat	=  as.double(betaMat),
      alphaMat =  as.double(alphaMat),
      ksiMat =  as.double(ksiMat),
      gammaMat =  as.double(gammaMat),
      probV1Mat =  as.double(probV1Mat),
      tau2Mat =  as.double(tau2Mat),
      sigma2Mat =  as.double(sigma2Mat),
      wMat =  as.double(wMat),
      likMat = as.double(likMat),
      logPostMat = as.double(logPostMat), PACKAGE ="spikeSlabGAM")


    beta <- matrix(res$betaMat, mcmc$chainLength, model$q)[, model$reverseOrder]
    colnames(beta) <- paste(model$groupIndicators,
      unlist(sapply(model$d, function(x) {
        return(1:x)
      })), sep =".")[model$reverseOrder]

    alpha <- matrix(res$alphaMat, mcmc$chainLength,  model$nGroups)
    colnames(alpha) <-  paste("alpha", levels(model$groupIndicators), sep =".")

    ksi <- matrix(res$ksiMat, mcmc$chainLength,  model$q)
    colnames(ksi) <-  paste("ksi", model$groupIndicators,
      unlist(sapply(model$d, function(x) {
        return(1:x)
      })), sep =".")


    tau <- matrix(res$tau2Mat, mcmc$chainLength, model$nPenGrps)
    gamma <- matrix(res$gammaMat, mcmc$chainLength, model$nPenGrps)
    pV1 <- matrix(res$probV1Mat, mcmc$chainLength, model$nPenGrps)
    colnames(tau) <- colnames(gamma) <-
      colnames(pV1) <-
      unique(model$groupIndicatorsOrig[properOrder[!(properOrder%in%unpen)]])
    reorderedPenGroupNames <- if(length(unpen)) {
      unique(model$groupIndicatorsOrig[-grep("u(.[0-9]+)?$",
        model$groupIndicatorsOrig)])
    } else unique(model$groupIndicatorsOrig)
    tau <- tau[, reorderedPenGroupNames, drop = F]
    gamma <- gamma[, reorderedPenGroupNames, drop = F]
    pV1 <- pV1[, reorderedPenGroupNames, drop = F ]
    colnames(tau) <- paste("tau.", colnames(tau), sep ="")
    colnames(gamma) <- paste("gamma.", colnames(gamma), sep ="")
    colnames(pV1) <- paste("p1.", colnames(pV1), sep ="")


    w <- matrix(res$wMat, mcmc$chainLength, 1)
    colnames(w) <- "w"
    sigma2 <- matrix(res$sigma2Mat, mcmc$chainLength, 1)
    colnames(sigma2) <- "sigma2"
    logLik <- matrix(res$likMat, mcmc$chainLength, 1)
    colnames(logLik) <- "logLik"
    logPost <- matrix(res$logPostMat, mcmc$chainLength, 1)
    colnames(logPost) <- "uLogPost"

    samples <- if(mcmc$reduceRet) {
      list(beta  = beta, gamma = gamma,
        pV1 = pV1, w = w, sigma2 = sigma2, logLik = logLik, logPost = logPost)
    } else {
      list(beta  = beta, alpha = alpha,
        ksi = ksi, tau = tau, gamma = gamma,
        pV1 = pV1, w = w, sigma2 = sigma2, logLik = logLik, logPost = logPost)
    }
    samples <- lapply(samples, mcmc, start = mcmc$burnin + 1,
      end = mcmc$totalLength, thin = mcmc$thin)

    accept <- if(family != 0) {
      list(alpha = res$acceptAlpha/mcmc$totalLength,
        ksi = res$acceptKsi/mcmc$totalLength)
    } else NULL
    posteriorPred <- if(mcmc$sampleY) {
      #if(mcmc$verbose) cat("\nsampling y from posterior predictive...\n")
      mu <- X %*% t(samples$beta) + model$offset
      yPred <- switch(familystr,
        gaussian = mu + t(as.vector(sqrt(samples$sigma2)) *
            matrix(rnorm(n * mcmc$chainLength),
              nrow = mcmc$chainLength)),
        binomial = t(matrix(rbinom(n * mcmc$chainLength,
          model$scale, plogis(mu)),
          nrow = mcmc$chainLength)),
        poisson = t(matrix(rpois(n * mcmc$chainLength, exp(mu)),
          nrow = mcmc$chainLength))
      )
      list(mu = mu, y = yPred)
    } else NULL

    restart <- list(beta = beta[mcmc$chainLength, , drop = F],
      tau = tau[mcmc$chainLength, , drop = F],
      gamma = gamma[mcmc$chainLength, , drop = F],
      w = w[mcmc$chainLength, , drop = F],
      sigma2 = sigma2[mcmc$chainLength, , drop = F])


    return(list(samples = samples, posteriorPred = posteriorPred,
      accept = accept, restart = restart))
  }

  if(parallel =="parallel") {
    res <- mclapply(1:mcmc$nChains, do1Chain)
  }
  if(parallel == "snow") {
    clusterExportLocal <-  function (cl, list) {
      for (name in list) {
        clusterCall(cl, assign, name, get(name, pos = -1))
      }
    }
    cl <- makeCluster(spec = options()$mc.cores, type ="PSOCK")
    clusterExportLocal(cl,
      c("hyperparameters","model","start","mcmc","y","X",
        "blocksAlpha","alphaIndABegin","alphaIndAEnd",
        "blocksKsi","ksiIndABegin","ksiIndAEnd","pcts","family",
        "betaMat", "alphaMat", "ksiMat", "gammaMat", "probV1Mat",
        "tau2Mat", "sigma2Mat", "wMat", "likMat", "logPostMat","unpen",
        "familystr"))
    res <- parLapply(cl, as.list(1:mcmc$nChains), do1Chain)
    stopCluster(cl)
  }
  if(parallel =="none") {
    res <- lapply(1:mcmc$nChains, function(x) {
      ret <- do1Chain(x)
      if (mcmc$verbose)
        cat("\n")
      return(ret)
    })
  }

  ##build the return object
  ret <- mget(names(formals()),
    envir = as.environment(-1))
  if(family != 0) {
    acceptance <- list(alpha = sapply(lapply(res, "[[", "accept"),
      "[[", "alpha", simplify = TRUE),
      ksi = sapply(lapply(res, "[[", "accept"),
        "[[", "ksi", simplify = TRUE))
    acceptanceWarning <- .15
    if(any(unlist(acceptance) < acceptanceWarning)) {
      cat("\nLow acceptance rates detected:\n")
      print(acceptance, digits = 2)

    }
    ret$mcmc$accept <- c(
      alpha = mean(acceptance$alpha),
      ksi = mean(acceptance$ksi))
    if(mcmc$verbose) {
      cat("\nMean acceptance rates:\n")
      print(ret$mcmc$accept, digits = 2)
    }
  }

  if(mcmc$sampleY) {
    ret$posteriorPred <- lapply(res, "[[", "posteriorPred")
  }

  ret$samples <- {
    tmp <- 	lapply(res, "[[", "samples")
    smpls <- vector(length(tmp[[1]]), mode ="list")
    names(smpls) <- names(tmp[[1]])
    for(n in names(smpls)) {
      smpls[[n]] <- do.call(mcmc.list, lapply(tmp, "[[", n))
    }
    rm(tmp)
    smpls
  }; rm(smpls)


  ret$postMeans <- if(mcmc$nChains>1) {
    lapply(ret$samples, function(x) {
      tmp <- sapply(x, colMeans)
      if(NCOL(tmp) > 1) {
        return(rowMeans(tmp))
      } else return(mean(tmp))
    })
  } else {
    lapply(ret$samples, function(x) {
      if(NCOL(x[[1]]) > 1) {
        return(sapply(x, colMeans))
      } else  {
        return(mean(x[[1]]))
      }
    })
  }


  ret$fitted <- cbind(eta = ret$X %*% ret$postMeans$beta)
  ret$fitted <- cbind(eta = ret$fitted, mu = switch(familystr,
    gaussian = ret$fitted,
    binomial = plogis(ret$fitted),
    poisson =  exp(ret$fitted)))


  ret$DIC <- {
    Dbar <- -2 * mean(unlist(ret$samples$logLik), na.rm = T)
    Dhat <- -2 * switch(familystr,
      gaussian = sum(dnorm(y - X %*% ret$postMeans$beta, mean = 0,
        sd = sqrt(ret$postMeans$sigma2), log = T)),
      binomial = sum(dbinom(y * model$scale, model$scale,
        plogis(X %*% ret$postMeans$beta), log = T)),
      poisson =  sum(dpois(y, exp(X %*% ret$postMeans$beta), log = T)))
    pD <- Dbar - Dhat
    DIC <- pD + Dbar
    c(DIC = DIC, pD = pD, Dbar = Dbar, Dhat = Dhat)
  }

  ret$X <- ret$X[, model$reverseOrder]

  if(mcmc$reduceRet) {
    ret$data <- NULL
  }

  runTime <- difftime(Sys.time(), startTime, units = "secs")
  ret$runTime <- runTime

  return(ret)
}
