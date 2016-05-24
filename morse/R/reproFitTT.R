#' Fits a Bayesian exposure-response model for target-time reproduction analysis
#'
#' This function estimates a model of the cumulated reproduction outputs of a
#' population in a given time period in presence of mortality.
#'
#' Because some individuals may die during the observation period, the
#' reproduction rate alone is not sufficient to account for the observed number
#' of offspring. In addition, we need the time individuals have stayed alive
#' during the experiment. The \code{reproFitTT} function estimates the number
#' of individual-days in an experiment between its start and the target time.
#' This covariable is then used to estimate a relation between the toxicant
#' concentration and the reproduction rate \emph{per individual-day}.
#'
#' The \code{reproFitTT} function fits two models, one where inter-individual
#' variability is neglected ("Poisson" model) and one where it is taken into
#' account ("gamma-Poisson" model). When setting \code{stoc.part} to
#' \code{"bestfit"}, a model comparison procedure is used to choose between
#' them. More details are presented in the vignette accompanying the package.
#'
#' @param data an object of class \code{reproData}
#' @param stoc.part stochastic part of the model. Possible values are \code{"bestfit"},
#' \code{"poisson"} and \code{"gammapoisson"}
#' @param target.time defines the observation period. By default the last time point
#' @param ecx desired values of \eqn{x} (in percent) for which to compute
#' \eqn{EC_{x}}{ECx}
#' @param n.chains number of MCMC chains. The minimum required number of chains is
#' @param quiet if \code{TRUE}, does not print messages and progress bars from JAGS
#'
#'
#' @return The function returns an object of class \code{reproFitTT} which is a list
#' of the following objects:
#' \item{DIC}{DIC value of the selected model}
#' \item{estim.ECx}{a table of the estimated 5, 10, 20 and 50 \% effective
#' concentrations (by default) and their 95 \% credible intervals}
#' \item{estim.par}{a table of the estimated parameters as medians and 95 \%
#' credible intervals}
#' \item{mcmc}{an object of class \code{mcmc.list} with the posterior distributions}
#' \item{model}{a JAGS model object}
#' \item{parameters}{a list of the parameters names used in the model}
#' \item{n.chains}{an integer value corresponding to the number of chains used
#' for the MCMC computation.}
#' \item{n.iter}{a list of two indices indicating the beginning and
#' the end of monitored iterations}
#' \item{n.thin}{a numerical value corresponding to the thinning interval}
#'
# FIXME
# \describe{
#
# Credible limits: For 100 values of concentrations regularly spread within
# the range of tested concentrations the joint posterior distribution of
# parameters is used to simulate 5000 values of \eqn{f_{ij}}, the number of
# offspring per individual-day for various replicates. For each concentration,
# 2.5, 50 and 97.5 percentiles of simulated values are calculated, from which
# there is a point estimate and a 95 \% credible interval (Delignette-Muller
# et al., 2014).
#
# DIC: The Deviance Information Criterium (DIC) as defined by Spiegelhalter et
# al. (2002) is provided by the \code{dic.samples} function. The DIC is a
# goodness-of-fit criterion penalized by the complexity of the model
# (Delignette-Muller et al., 2014).
#
# Raftery and Lewis's diagnostic: The \code{raftery.diag} is a run length
# control diagnostic based on a criterion that calculates the appropriate
# number of iterations required to accurately estimate the parameter
# quantiles. The Raftery and Lewis's diagnostic value used in the
# \code{reproFitTT} function is the \code{resmatrix} object. See the
# \code{\link[coda]{raftery.diag}} help for more details.
#
# Model selection: When \code{stoc.part = "bestfit"}, the \code{reproFitTT}
# function chooses itself between the Poisson and the Gamma-Poisson model
# depending on the number of MCMC samples and on the DIC values.  The minimum
# number of MCMC samples for the pilot run is provided by the Raftery and
# Lewis's diagnostic (Raftery and Lewis 1992). If this number is less than 100
# 000 for the Poisson and the Gamma-Poisson model and if the DIC difference
# between Poisson and Gamma-poisson models is small (typically less than 10),
# then the Poisson model is selected. If this number is more than 100 000 for
# only one model, the other one is selected.
# }
# @seealso \code{\link[rjags]{rjags}}, \code{\link[rjags]{coda.samples}},
# \code{\link[rjags]{dic.samples}}, \code{\link[coda]{raftery.diag}}
# and \code{\link[ggplot2]{ggplot}}
#
# @references Delignette-Muller, M.L., Lopes, C., Veber, P. and Charles, S.
# (2014) Statistical handling of reproduction data for exposure-response
# modelling.
# \url{http://pubs.acs.org/doi/abs/10.1021/es502009r?journalCode=esthag}.
#
# Plummer, M. (2013) JAGS Version 4.0.0 user manual.
# \url{http://sourceforge.net/projects/mcmc-jags/files/Manuals/4.x/jags_user_manual.pdf/download}
#
# Raftery A.E. and Lewis, S.M. (1992) One long run with diagnostics:
# Implementation strategies for Markov chain Monte Carlo. \emph{Statistical
# Science}, 7, 493-497.
#
# Spiegelhalter, D., N. Best, B. Carlin, and A. van der Linde (2002) Bayesian
# measures of model complexity and fit (with discussion).  \emph{Journal of
# the Royal Statistical Society}, Series B 64, 583-639.
#
#'
#' @keywords estimation
#'
#' @examples
#'
#' # (1) Load the data
#' data(cadmium1)
#'
#' # (2) Create an object of class "reproData"
#' dat <- reproData(cadmium1)
#'
#' \dontrun{
#' # (3) Run the reproFitTT function with the log-logistic gamma-poisson model
#' out <- reproFitTT(dat, stoc.part = "gammapoisson",
#'                   ecx = c(5, 10, 15, 20, 30, 50, 80), quiet = TRUE)
#'
#' # (4) Summary look the estimated values (ECx and parameters)
#' out$estim.ECx
#' out$estim.par
#'
#' # (5) Plot the fitted curve with credible limits
#' plot(out, log.scale = TRUE, ci = TRUE,
#'      main = "log-logistic gamma-poisson model")
#'
#' # (6) Plot the fitted curve with ggplot style
#' require("ggplot2")
#' plot(out, xlab = expression("Concentration in" ~ mu~g.L^{-1}),
#'      fitcol = "blue", ci = TRUE, cicol = "blue", style = "ggplot",
#'      main = "Log-logistic response to concentration")
#' }
#'
#' @import rjags
#' 
#' @export
reproFitTT <- function(data,
                       stoc.part = "bestfit",
                       target.time = NULL,
                       ecx = c(5, 10, 20, 50),
                       n.chains = 3,
                       quiet = FALSE) {
  # test class object
  if (! is(data, "reproData"))
    stop("reproFitTT: object of class reproData expected")

  # stocastic verification
  stoc.partpossible <- c("poisson", "gammapoisson", "bestfit")

  if (!any(stoc.partpossible == stoc.part))
    stop("Invalid value for argument [stoc.part]")

  # check 0 Nreprocumul
  if (all(data$Nreprocumul == 0))
    stop("Nreprocumul contains only 0 values !")

  # parameters
  parameters <- list(poisson = c("d", "log10b", "log10e"),
                     gammapoisson = c("d", "log10b","log10e", "log10omega"))

  # select Data at target.time
  dataTT <- selectDataTT(data, target.time)

  # create priors parameters
  jags.data <- reproCreateJagsData(stoc.part, dataTT)

  # Poisson model only
  if (stoc.part == "poisson") {
    # Define model
    poisson.model <- reproLoadPoissonModel(model.program = llm.poisson.model.text,
                                           data = jags.data,
                                           n.chains, quiet)

    # Determine sampling parameters
    poisson.sampling.parameters <- modelSamplingParameters(poisson.model,
                                                           parameters$poisson,
                                                           n.chains, quiet)

    if (poisson.sampling.parameters$niter > 100000)
      stop("The model needs too many iterations to provide reliable parameter estimates !")

    # calcul DIC
    poisson.DIC <- calcDIC(poisson.model, poisson.sampling.parameters, quiet)

    # list of objet for the coda.sample function
    coda.arg <- list(model = poisson.model,
                     model.label = "P",
                     niter = poisson.sampling.parameters$niter,
                     thin = poisson.sampling.parameters$thin,
                     nburnin = poisson.sampling.parameters$burnin,
                     parameters = parameters$poisson,
                     DIC = poisson.DIC)
  }

  # Gamma-poisson model only
  if (stoc.part == "gammapoisson") {
    # Define model
    gammapoisson.model <- reproLoadGammapoissonModel(model.program = llm.gammapoisson.model.text,
                                                     data = jags.data,
                                                     n.chains, quiet)

    # Determine sampling parameters
    gammapoisson.sampling.parameters <- modelSamplingParameters(gammapoisson.model,
                                                                parameters$gammapoisson,
                                                                n.chains, quiet)

    if (gammapoisson.sampling.parameters$niter > 100000)
      stop("The model needs too many iterations to provide reliable parameter estimates !")

    # calcul DIC
    gammapoisson.DIC <- calcDIC(gammapoisson.model,
                                gammapoisson.sampling.parameters, quiet)

    # list of objet for the coda.sample function
    coda.arg <- list(model = gammapoisson.model,
                     model.label = "GP",
                     niter = gammapoisson.sampling.parameters$niter,
                     thin = gammapoisson.sampling.parameters$thin,
                     nburnin = gammapoisson.sampling.parameters$burnin,
                     parameters = parameters$gammapoisson,
                     DIC = gammapoisson.DIC)
  }

  # Model Selection by the DIC
  if (stoc.part == "bestfit") {
    # Define models
    poisson.model <- reproLoadPoissonModel(model.program = llm.poisson.model.text,
                                           data = jags.data,
                                           n.chains, quiet)

    gammapoisson.model <- reproLoadGammapoissonModel(model.program = llm.gammapoisson.model.text,
                                                     data = jags.data,
                                                     n.chains, quiet)
    # Determine sampling parameters
    poisson.sampling.parameters <- modelSamplingParameters(poisson.model,
                                                           parameters$poisson,
                                                           n.chains, quiet)

    gammapoisson.sampling.parameters <- modelSamplingParameters(gammapoisson.model,
                                                                parameters$gammapoisson,
                                                                n.chains, quiet)

    if (poisson.sampling.parameters$niter > 100000 && gammapoisson.sampling.parameters$niter > 100000)
      stop("The model needs too many iterations to provide reliable parameter estimates !")

    # calcul DIC
    poisson.DIC <- calcDIC(poisson.model, poisson.sampling.parameters, quiet)
    gammapoisson.DIC <- calcDIC(gammapoisson.model,
                                gammapoisson.sampling.parameters, quiet)

    if (gammapoisson.sampling.parameters$niter > 100000) {
      # list of object for the coda.sample function
      coda.arg <- list(model = poisson.model,
                       model.label = "P",
                       niter = poisson.sampling.parameters$niter,
                       thin = poisson.sampling.parameters$thin,
                       nburnin = poisson.sampling.parameters$burnin,
                       parameters = parameters$poisson,
                       DIC = poisson.DIC)
    }

    if (poisson.sampling.parameters$niter > 100000) {
      # list of object for the coda.sample function
      coda.arg <- list(model = gammapoisson.model,
                       model.label = "GP",
                       niter = gammapoisson.sampling.parameters$niter,
                       thin = gammapoisson.sampling.parameters$thin,
                       nburnin = gammapoisson.sampling.parameters$burnin,
                       parameters = parameters$gammapoisson,
                       DIC = gammapoisson.DIC)
    }
    if (poisson.sampling.parameters$niter <= 100000 && gammapoisson.sampling.parameters$niter <= 100000) {
      if (poisson.DIC <= (gammapoisson.DIC + 10)) {
        # list of objet for the coda.sample function
        coda.arg <- list(model = poisson.model,
                         model.label = "P",
                         niter = poisson.sampling.parameters$niter,
                         thin = poisson.sampling.parameters$thin,
                         nburnin = poisson.sampling.parameters$burnin,
                         parameters = parameters$poisson,
                         DIC = poisson.DIC)
      } else {
        # list of objet for the coda.sample function
        coda.arg <- list(model = gammapoisson.model,
                         model.label = "GP",
                         niter = gammapoisson.sampling.parameters$niter,
                         thin = gammapoisson.sampling.parameters$thin,
                         nburnin = gammapoisson.sampling.parameters$burnin,
                         parameters = parameters$gammapoisson,
                         DIC = gammapoisson.DIC)
      }
    }
  }

  # Sampling
  prog.b <- ifelse(quiet == TRUE, "none", "text")
  mcmc <- coda.samples(coda.arg$model,
                       coda.arg$parameters,
                       n.iter = coda.arg$niter,
                       thin = coda.arg$thin,
                       progress.bar = prog.b)

  # summarize estime.par et CIs
  # calculate from the estimated parameters
  estim.par <- reproPARAMS(mcmc, coda.arg$model.label)

  # ECx calculation  estimated ECx and their CIs 95%
  # vector of ECX
  estim.ECx <- estimXCX(mcmc, ecx, "EC")

  # check if the maximum measured concentration is in the EC50's range of
  # 95% percentile
  if (50 %in% ecx) {
    EC50 <- log10(estim.ECx["EC50", "median"])
    if (!(min(log10(data$conc)) < EC50 & EC50 < max(log10(data$conc))))
      warning("The EC50 estimation lies outsides the range of tested concentration and may be unreliable !")
  }

  # output
  OUT <- list(DIC = coda.arg$DIC,
              estim.ECx = estim.ECx,
              estim.par = estim.par,
              det.part = "loglogistic",
              mcmc = mcmc,
              model = coda.arg$model,
              model.label = coda.arg$model.label,
              parameters = coda.arg$parameters,
              n.chains = summary(mcmc)$nchain,
              n.iter = list(start = summary(mcmc)$start,
                            end = summary(mcmc)$end),
              n.thin = summary(mcmc)$thin,
              jags.data = jags.data,
              transformed.data = data,
              dataTT = dataTT)

  class(OUT) <- "reproFitTT"
  return(OUT)
}


#' @importFrom stats sd
reproCreateJagsData <- function(stoc.part, data) {
  # create the parameters to define the prior of the log-logistic model
  # for reproduction data analysis
  # INPUTS
  # stoc.part: model name
  # data: object of class reproData
  # OUTPUT
  # jags.data : list data require for the jags.model function


  # separate control data to the other
  # tab0: data at conc = 0
  tab0 <- data[data$conc == min(data$conc), ]
  # tab: data at conc != 0
  tab <- data[data$conc != min(data$conc), ]

  Nindtime <- tab$Nindtime
  NreprocumulIndtime0 <- tab0$Nreprocumul / tab0$Nindtime # cumulated number of
  # offspring / number of
  # individual-days
  conc <- tab$conc
  Ncumul <- tab$Nreprocumul
  n <- nrow(tab) # number of observation != from the control

  # Parameter calculation of concentration min and max
  concmin <- min(sort(unique(conc))[-1])
  concmax <- max(conc)

  # create priors parameters for the log logistic model

  # Params to define log10e
  meanlog10e <- (log10(concmin) + log10(concmax)) / 2
  sdlog10e <- (log10(concmax) - log10(concmin)) / 4
  taulog10e <- 1 / sdlog10e^2

  # Params to define d
  meand <- mean(NreprocumulIndtime0)
  SEd <- sd(NreprocumulIndtime0) / sqrt(length(unique(tab0$replicate)))
  taud <- 1 / (SEd)^2

  # Params to define b
  log10bmin <- -2
  log10bmax <- 2

  # list of data use by jags
  jags.data <- list(meanlog10e = meanlog10e,
                    taulog10e = taulog10e,
                    meand = meand,
                    taud = taud,
                    log10bmin = log10bmin,
                    log10bmax = log10bmax,
                    n = n,
                    xconc = conc,
                    Nindtime = Nindtime,
                    Ncumul = Ncumul)

  # Params to define overdispersion rate
  if (stoc.part == "bestfit" || stoc.part == "gammapoisson") {
    log10omegamin <- -4
    log10omegamax <- 4

    # list of data use by jags
    jags.data <- c(jags.data,
                   log10omegamin = log10omegamin,
                   log10omegamax = log10omegamax)
  }
  return(jags.data)
}

reproLoadPoissonModel <- function(model.program,
                                  data,
                                  n.chains,
                                  quiet = quiet) {
                                    # sub function to load jags poisson model
                                    reproLoadModel(model.program, F, data, n.chains, quiet = quiet)
                                  }

reproLoadGammapoissonModel <- function(model.program,
                                       data,
                                       n.chains,
                                       quiet = quiet) {
                                         # sub function to load jags gamma poisson model
                                         reproLoadModel(model.program, T, data, n.chains, quiet = quiet)
}

#' @import rjags
reproLoadModel <- function(model.program,
                           lr.bound.keep,
                           data,
                           n.chains,
                           Nadapt = 3000,
                           quiet = quiet) {
  # create the JAGS model object and called by reproLoadPoissonModel
  # and reproLoadGammapoissonModel
  # INPUTS:
  # - model.program: character string containing a jags model description
  # - lr.bound.keep: boolean value to use omega parameter or not
  # - data: list of data created by reproCreateJagsData
  # - nchains: Number of chains desired
  # - Nadapt: length of the adaptation phase
  # - quiet: silent option
  # OUTPUT:
  # - JAGS model

  # delisting of lr.bound because not used in the function
  if (!lr.bound.keep) {
    data[c("meanlog10omega", "taulog10omega", "log10omegamin",
           "log10omegamax")] <- NULL
    }

    # load model text in a temporary file
    model.file <- tempfile() # temporary file address
    fileC <- file(model.file) # open connection
    writeLines(model.program, fileC) # write text in temporary file
    close(fileC) # close connection to temporary file
    # creation of the jags model
    model <- jags.model(file = model.file, data = data, n.chains = n.chains,
                        n.adapt = Nadapt, quiet = quiet)
    unlink(model.file)
    return(model)
}

reproPARAMS <- function(mcmc, MODEL = "P") {
  # create the table of posterior estimated parameters
  # for the reproduction analyses
  # INPUT:
  # - mcmc:  list of estimated parameters for the model with each item representing
  # a chain
  # - MODEL: a position flag model with P: poisson model and GP: gammapoisson
  # model
  # OUTPUT:
  # - data frame with 3 columns (values, CIinf, CIsup) and 3-4rows (the estimated
  # parameters)

  # Retrieving parameters of the model
  res.M <- summary(mcmc)

  b <- 10^res.M$quantiles["log10b", "50%"]
  d <- res.M$quantiles["d", "50%"]
  e <- 10^res.M$quantiles["log10e", "50%"]
  binf <- 10^res.M$quantiles["log10b", "2.5%"]
  dinf <- res.M$quantiles["d", "2.5%"]
  einf <- 10^res.M$quantiles["log10e", "2.5%"]
  bsup <- 10^res.M$quantiles["log10b", "97.5%"]
  dsup <- res.M$quantiles["d", "97.5%"]
  esup <- 10^res.M$quantiles["log10e", "97.5%"]

  # Definition of the parameter storage and storage data

  # If Poisson Model
  if (MODEL == "P") {
    rownames <- c("b", "d", "e")
    params <- c(b, d, e)
    CIinf <- c(binf, dinf, einf)
    CIsup <- c(bsup, dsup, esup)
  }
  # If Gamma Poisson Model
  if (MODEL == "GP") {
    # Calculation of the parameter omega
    omega <- 10^res.M$quantiles["log10omega", "50%"]
    omegainf <- 10^res.M$quantiles["log10omega", "2.5%"]
    omegasup <- 10^res.M$quantiles["log10omega", "97.5%"]
    # Definition of the parameter storage and storage data
    rownames <- c("b", "d", "e", "omega")
    params <- c(b, d, e, omega)
    CIinf <- c(binf, dinf, einf, omegainf)
    CIsup <- c(bsup, dsup, esup, omegasup)
  }

  res <- data.frame(median = params, Q2.5 = CIinf, Q97.5 = CIsup,
                    row.names = rownames)

  return(res)
}

llm.poisson.model.text <- "\nmodel # Loglogistic Poisson model\n{\n#\nfor (j in 1:n) # loop on replicates\n{\n# Explicit writting of a Poisson law for each replicate\n# mean is given by the theoretical curve\nytheo[j] <- d / (1 + pow(xconc[j]/e, b))\nnbtheo[j] <- ytheo[j]*Nindtime[j]\nNcumul[j] ~ dpois(nbtheo[j])\n}\n# Prior distributions\nd ~ dnorm(meand, taud)T(0,)\nlog10b ~ dunif(log10bmin, log10bmax)\nlog10e ~ dnorm(meanlog10e, taulog10e)\n\nb <- pow(10,log10b)\ne <- pow(10,log10e)\n}\n"

llm.gammapoisson.model.text <- "\nmodel # Loglogisitc Gamma poisson model\n{\n#\nfor (j in 1:n) # loop on replicates\n{\n# Explicit writting of a gamma-Poisson law for each replicate\n# the mean is given by a gamma law centered on the theoretical curve\nrate[j] <- d / (1 + pow(xconc[j]/e, b)) / omega\np[j] <- 1 / (Nindtime[j] * omega + 1)\nNcumul[j] ~ dnegbin(p[j], rate[j])\n}\n# Prior distributions\nd ~ dnorm(meand, taud)T(0,)\nlog10b ~ dunif(log10bmin, log10bmax)\nlog10e ~ dnorm(meanlog10e, taulog10e)\nlog10omega ~ dunif(log10omegamin, log10omegamax)\n\nomega <- pow(10,log10omega)\nb <- pow(10,log10b)\ne <- pow(10,log10e)\n}\n"

