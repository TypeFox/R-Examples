#' Apply the (un-)logged minimum age model (MAM) after Galbraith et al. (1999)
#' to a given De distribution
#'
#' Function to fit the (un-)logged three or four parameter minimum dose model
#' (MAM-3/4) to De data.
#'
#' \bold{Parameters} \cr\cr This model has four parameters: \cr\cr
#' \tabular{rl}{ \code{gamma}: \tab minimum dose on the log scale \cr
#' \code{mu}: \tab mean of the non-truncated normal distribution \cr
#' \code{sigma}: \tab spread in ages above the minimum \cr \code{p0}: \tab
#' proportion of grains at gamma \cr } If \code{par=3} (default) the
#' 3-parametric minimum age model is applied, where \code{gamma=mu}. For
#' \code{par=4} the 4-parametric model is applied instead.\cr\cr
#' \bold{(Un-)logged model} \cr\cr In the original version of the
#' three-parameter minimum dose model, the basic data are the natural
#' logarithms of the De estimates and relative standard errors of the De
#' estimates. This model will be applied if \code{log=TRUE}. \cr\cr If
#' \code{log=FALSE}, the modified un-logged model will be applied instead. This
#' has essentially the same form as the original version.  \code{gamma} and
#' \code{sigma} are in Gy and \code{gamma} becomes the minimum true dose in the
#' population. \cr\cr While the original (logged) version of the mimimum dose
#' model may be appropriate for most samples (i.e. De distributions), the
#' modified (un-logged) version is specially designed for modern-age and young
#' samples containing negative, zero or near-zero De estimates (Arnold et al.
#' 2009, p. 323). \cr\cr \bold{Initial values & boundaries} \cr\cr The log
#' likelihood calculations use the \link{nlminb} function for box-constrained
#' optimisation using PORT routines.  Accordingly, initial values for the four
#' parameters can be specified via \code{init.values}. If no values are
#' provided for \code{init.values} reasonable starting values are estimated
#' from the input data.  If the final estimates of \emph{gamma}, \emph{mu},
#' \emph{sigma} and \emph{p0} are totally off target, consider providing custom
#' starting values via \code{init.values}. \cr In contrast to previous versions
#' of this function the boundaries for the individual model parameters can no
#' longer be specified. Appropriate boundary are now hard-coded and are valid
#' for all input data sets. \cr\cr \bold{Bootstrap} \cr\cr When
#' \code{bootstrap=TRUE} the function applies the bootstrapping method as
#' described in Wallinga & Cunningham (2012). By default, the minimum age model
#' produces 1000 first level and 3000 second level bootstrap replicates
#' (actually, the number of second level bootstrap replicates is three times
#' the number of first level replicates unless specified otherwise).  The
#' uncertainty on sigmab is 0.04 by default. These values can be changed by
#' using the arguments \code{bs.M} (first level replicates), \code{bs.N}
#' (second level replicates) and \code{sigmab.sd} (error on sigmab). With
#' \code{bs.h} the bandwidth of the kernel density estimate can be specified.
#' By default, \code{h} is calculated as \cr \deqn{h =
#' (2*\sigma_{DE})/\sqrt{n}} \cr \bold{Multicore support} \cr\cr This function
#' supports parallel computing and can be activated by \code{multicore=TRUE}.
#' By default, the number of available logical CPU cores is determined
#' automatically, but can be changed with \code{cores}. The multicore support
#' is only available when \code{bootstrap=TRUE} and spawns \code{n} R instances
#' for each core to get MAM estimates for each of the N and M boostrap
#' replicates. Note that this option is highly experimental and may or may not
#' work for your machine. Also the performance gain increases for larger number
#' of bootstrap replicates. Also note that with each additional core and hence
#' R instance and depending on the number of bootstrap replicates the memory
#' usage can significantly increase. Make sure that memory is always availabe,
#' otherwise there will be a massive perfomance hit.
#'
#' @param data \code{\linkS4class{RLum.Results}} or \link{data.frame}
#' (\bold{required}): for \code{data.frame}: two columns with De \code{(data[
#' ,1])} and De error \code{(values[ ,2])}
#' @param sigmab \code{\link{numeric}} (\bold{required}): spread in De values
#' given as a fraction (e.g. 0.2). This value represents the expected
#' overdispersion in the data should the sample be well-bleached (Cunningham &
#' Walling 2012, p. 100).
#' @param log \code{\link{logical}} (with default): fit the (un-)logged minimum
#' dose model to De data
#' @param par \code{\link{numeric}} (with default): apply the 3- or
#' 4-parametric minimum age model (\code{par=3} or \code{par=4}). The MAM-3 is
#' used by default.
#' @param bootstrap \code{\link{logical}} (with default): apply the recycled
#' bootstrap approach of Cunningham & Wallinga (2012).
#' @param init.values \code{\link{numeric}} (optional): a named list with
#' starting values for gamma, sigma, p0 and mu (e.g. \code{list(gamma=100
#' sigma=1.5, p0=0.1, mu=100)}). If no values are provided reasonable values
#' are tried to be estimated from the data.
#' @param plot \code{\link{logical}} (with default): plot output
#' (\code{TRUE}/\code{FALSE})
#' @param multicore \code{\link{logical}} (with default): enable parallel
#' computation of the bootstrap by creating a multicore SNOW cluster. Depending
#' on the number of available logical CPU cores this will drastically reduce
#' the computation time. Note that this option is highly experimental and not
#' work for all machines. (\code{TRUE}/\code{FALSE})
#' @param \dots (optional) further arguments for bootstrapping (\code{bs.M,
#' bs.N, bs.h, sigmab.sd}).  See details for their usage. Further arguments are
#' \code{verbose} to de-/activate console output (logical), \code{debug} for
#' extended console output (logical) and \code{cores} (integer) to manually
#' specify the number of cores to be used when \code{multicore=TRUE}.
#' @return Returns a plot (optional) and terminal output. In addition an
#' \code{\linkS4class{RLum.Results}} object is returned containing the
#' following elements:
#'
#' \item{summary}{\link{data.frame} summary of all relevant model results.}
#' \item{data}{\link{data.frame} original input data} \item{args}{\link{list}
#' used arguments} \item{call}{\link{call} the function call}
#' \item{mle}{\link{mle2} object containing the maximum log likelhood functions
#' for all parameters} \item{BIC}{\link{numeric} BIC score}
#' \item{confint}{\link{data.frame} confidence intervals for all parameters}
#' \item{profile}{\link{profile.mle2} the log likelihood profiles}
#' \item{bootstrap}{\link{list} bootstrap results}
#'
#' The output should be accessed using the function
#' \code{\link{get_RLum}}
#' @note The default starting values for \emph{gamma}, \emph{mu}, \emph{sigma}
#' and \emph{p0} may only be appropriate for some De data sets and may need to
#' be changed for other data. This is especially true when the un-logged
#' version is applied. \cr Also note that all R warning messages are suppressed
#' when running this function. If the results seem odd consider re-running the
#' model with \code{debug=TRUE} which provides extended console output and
#' forwards all internal warning messages.
#' @section Function version: 0.4.1
#' @author Christoph Burow, University of Cologne (Germany) \cr Based on a
#' rewritten S script of Rex Galbraith, 2010 \cr The bootstrap approach is
#' based on a rewritten MATLAB script of Alastair Cunningham. \cr Alastair
#' Cunningham is thanked for his help in implementing and cross-checking the
#' code.
#' @seealso \code{\link{calc_CentralDose}}, \code{\link{calc_CommonDose}},
#' \code{\link{calc_FiniteMixture}}, \code{\link{calc_FuchsLang2001}},
#' \code{\link{calc_MaxDose}}
#' @references Arnold, L.J., Roberts, R.G., Galbraith, R.F. & DeLong, S.B.,
#' 2009. A revised burial dose estimation procedure for optical dating of young
#' and modern-age sediments. Quaternary Geochronology 4, 306-325. \cr\cr
#' Galbraith, R.F. & Laslett, G.M., 1993. Statistical models for mixed fission
#' track ages. Nuclear Tracks Radiation Measurements 4, 459-470. \cr\cr
#' Galbraith, R.F., Roberts, R.G., Laslett, G.M., Yoshida, H. & Olley, J.M.,
#' 1999. Optical dating of single grains of quartz from Jinmium rock shelter,
#' northern Australia. Part I: experimental design and statistical models.
#' Archaeometry 41, 339-364. \cr\cr Galbraith, R.F., 2005. Statistics for
#' Fission Track Analysis, Chapman & Hall/CRC, Boca Raton. \cr\cr Galbraith,
#' R.F. & Roberts, R.G., 2012. Statistical aspects of equivalent dose and error
#' calculation and display in OSL dating: An overview and some recommendations.
#' Quaternary Geochronology 11, 1-27. \cr\cr \bold{Further reading} \cr\cr
#' Arnold, L.J. & Roberts, R.G., 2009. Stochastic modelling of multi-grain
#' equivalent dose (De) distributions: Implications for OSL dating of sediment
#' mixtures. Quaternary Geochronology 4, 204-230. \cr\cr Bailey, R.M. & Arnold,
#' L.J., 2006. Statistical modelling of single grain quartz De distributions
#' and an assessment of procedures for estimating burial dose. Quaternary
#' Science Reviews 25, 2475-2502. \cr\cr Cunningham, A.C. & Wallinga, J., 2012.
#' Realizing the potential of fluvial archives using robust OSL chronologies.
#' Quaternary Geochronology 12, 98-106. \cr\cr Rodnight, H., Duller, G.A.T.,
#' Wintle, A.G. & Tooth, S., 2006. Assessing the reproducibility and accuracy
#' of optical dating of fluvial deposits.  Quaternary Geochronology 1, 109-120.
#' \cr\cr Rodnight, H., 2008. How many equivalent dose values are needed to
#' obtain a reproducible distribution?. Ancient TL 26, 3-10. \cr\cr
#' @examples
#'
#'
#' ## Load example data
#' data(ExampleData.DeValues, envir = environment())
#'
#' # (1) Apply the minimum age model with minimum required parameters.
#' # By default, this will apply the un-logged 3-parametric MAM.
#' calc_MinDose(data = ExampleData.DeValues$CA1, sigmab = 0.1)
#'
#' # (2) Re-run the model, but save results to a variable and turn
#' # plotting of the log-likelihood profiles off.
#' mam <- calc_MinDose(data = ExampleData.DeValues$CA1,
#'                     sigmab = 0.1,
#'                     plot = FALSE)
#'
#' # Show structure of the RLum.Results object
#' mam
#'
#' # Show summary table that contains the most relevant results
#' res <- get_RLum(mam, "summary")
#' res
#'
#' # Plot the log likelihood profiles retroactively, because before
#' # we set plot = FALSE
#' plot_RLum.Results(mam)
#'
#' # Plot the dose distribution in an abanico plot and draw a line
#' # at the minimum dose estimate
#' plot_AbanicoPlot(data = ExampleData.DeValues$CA1,
#'                  main = "3-parameter Minimum Age Model",
#'                  line = mam,polygon.col = "none",
#'                  hist = TRUE,
#'                  rug = TRUE,
#'                  summary = c("n", "mean", "mean.weighted", "median", "in.ci"),
#'                  centrality = res$de,
#'                  line.col = "red",
#'                  grid.col = "none",
#'                  line.label = paste0(round(res$de, 1), "\U00B1",
#'                                      round(res$de_err, 1), " Gy"),
#'                  bw = 0.1,
#'                  ylim = c(-25, 18),
#'                  summary.pos = "topleft",
#'                  mtext = bquote("Parameters: " ~
#'                                   sigma[b] == .(get_RLum(mam, "args")$sigmab) ~ ", " ~
#'                                   gamma == .(round(log(res$de), 1)) ~ ", " ~
#'                                   sigma == .(round(res$sig, 1)) ~ ", " ~
#'                                   rho == .(round(res$p0, 2))))
#'
#' # (3) Run the minimum age model with bootstrap
#' # NOTE: Bootstrapping is computationally intensive, which is why the
#' # following example is commented out. To run the examples just
#' # uncomment the code.
#' # (3.1) run the minimum age model with default values for bootstrapping
#' #calc_MinDose(data = ExampleData.DeValues$CA1,
#' #             sigmab = 0.15,
#' #             bootstrap = TRUE)
#'
#' # (3.2) Bootstrap control parameters
#' #mam <- calc_MinDose(data = ExampleData.DeValues$CA1,
#' #                    sigmab = 0.15,
#' #                    bootstrap = TRUE,
#' #                    bs.M = 300,
#' #                    bs.N = 500,
#' #                    bs.h = 4,
#' #                    sigmab.sd = 0.06,
#' #                    plot = FALSE)
#'
#' # Plot the results
#' #plot_RLum(mam)
#'
#' # save bootstrap results in a separate variable
#' #bs <- get_RLum(mam, "bootstrap")
#'
#' # show structure of the bootstrap results
#' #str(bs, max.level = 2, give.attr = FALSE)
#'
#' # print summary of minimum dose and likelihood pairs
#' #summary(bs$pairs$gamma)
#'
#' # Show polynomial fits of the bootstrap pairs
#' #bs$poly.fits$poly.three
#'
#' # Plot various statistics of the fit using the generic plot() function
#' #par(mfcol=c(2,2))
#' #plot(bs$poly.fits$poly.three, ask = FALSE)
#'
#' # Show the fitted values of the polynomials
#' #summary(bs$poly.fits$poly.three$fitted.values)
#'
#' @export
calc_MinDose <- function(
  data,
  sigmab,
  log = TRUE,
  par = 3,
  bootstrap = FALSE,
  init.values,
  plot = TRUE,
  multicore = FALSE,
  ...
){

  ##============================================================================##
  ## ... ARGUMENTS
  ##============================================================================##

  extraArgs <- list(...)

  ## check if this function is called by calc_MaxDose()
  if ("invert" %in% names(extraArgs)) {
    invert <- extraArgs$invert
    if (!log) {
      log <- TRUE # overwrite user choice as max dose model currently only supports the logged version
      cat(paste("\n[WARNING] The maximum dose model only supports the logged version.",
                "'log' was automatically changed to TRUE.\n\n"))
    }
  } else {
    invert <- FALSE
  }

  ## console output
  if ("verbose" %in% names(extraArgs)) {
    verbose <- extraArgs$verbose
  } else {
    verbose <- TRUE
  }

  ## bootstrap replications
  # first level bootstrap
  if ("bs.M" %in% names(extraArgs)) {
    M <- as.integer(extraArgs$bs.M)
  } else {
    M <- 1000
  }

  # second level bootstrap
  if ("bs.N" %in% names(extraArgs)) {
    N <- as.integer(extraArgs$bs.N)
  } else {
    N <- 3*M
  }

  # KDE bandwith
  if ("bs.h" %in% names(extraArgs)) {
    h <- extraArgs$bs.h
  } else {
    h <- (sd(data[ ,1])/sqrt(length(data[ ,1])))*2
  }

  # standard deviation of sigmab
  if ("sigmab.sd" %in% names(extraArgs)) {
    sigmab.sd <- extraArgs$sigmab.sd
  } else {
    sigmab.sd <- 0.04
  }

  if ("debug" %in% names(extraArgs)) {
    debug <- extraArgs$debug
  } else {
    debug <- FALSE
  }

  if ("cores" %in% names(extraArgs)) {
    cores <- extraArgs$cores
  } else {
    cores <- parallel::detectCores()
    if (multicore)
      message(paste("Logical CPU cores detected:", cores))
  }

  ## WARNINGS ----
  if (!debug)
    options(warn = -1)

  ##============================================================================##
  ## START VALUES
  ##============================================================================##

  if (missing(init.values)) {
    start <- list(gamma = ifelse(log, log(quantile(data[ ,1], probs = 0.25)),
                                 quantile(data[ ,1], probs = 0.25)),
                  sigma = 1.2,
                  p0 = 0.01,
                  mu = ifelse(log, log(quantile(data[ ,1], probs = 0.25)),
                              mean(data[ ,1])))
  } else {
    start <- list(gamma = init.values$gamma,
                  sigma = init.values$sigma,
                  p0 = init.values$p0,
                  mu = init.values$mu)
  }

  ##============================================================================##
  ## ESTIMATE BOUNDARY PARAMETERS
  ##============================================================================##

  gamma.xlb <- min(data[ ,1]/10)
  gamma.xub <- max(data[ ,1]*1.1)
  sigma.xlb <- 0
  sigma.xub <- 5
  mu.xlb <- min(data[ ,1])/10
  mu.xub <- max(data[ ,1]*1.1)

  # combine lower and upper boundary values to vectors
  if (log) {
    xlb <- c(log(gamma.xlb), sigma.xlb, 0)
    xub <- c(log(gamma.xub), sigma.xub, 1)
  } else {
    xlb <- c(gamma.xlb, sigma.xlb, 0)
    xub <- c(gamma.xub, exp(sigma.xub), 1)
  }
  if (par==4) {
    xlb <- c(xlb, ifelse(log, log(mu.xlb), mu.xlb))
    xub <- c(xub, ifelse(log, log(mu.xub), mu.xub))
  }

  ##============================================================================##
  ## AUXILLARY FUNCTIONS
  ##============================================================================##

  # THIS FUNCTION CALCULATES THE NEGATIVE LOG LIKELIHOOD OF THE DATA
  Neglik_f <- function(gamma, sigma, p0, mu, data) {
    # this calculates the negative of the log likelihood of the
    # data (data) for a given set of parameters (gamma, sigma, p0)
    # data is a 2x2 matrix of data: De, rel_error (including sigma_b)

    # recover the data
    zi <- data[ ,1]
    si <- data[ ,2]
    n <- length(zi)

    # in the MAM-3 gamma and mu are assumed to be equal
    if (par == 3)
      mu <- gamma

    # calculate sigma^2 + seld^2, mu0 and sigma0
    s2 <- sigma^2 + si^2
    sigma0 <- 1/sqrt(1/sigma^2 + 1/si^2)
    mu0 <- (mu/sigma^2 + zi/si^2)/(1/sigma^2 + 1/si^2)

    # calculate the log-likelihood
    logsqrt2pi <- 0.5*log(2*pi)
    res0 <- (gamma - mu0)/sigma0
    res1 <- (gamma - mu)/sigma
    lf1i <- log(p0) - log(si) - 0.5*((zi-gamma)/si)^2   - logsqrt2pi
    lf2i <- log(1-p0) - 0.5*log(s2) - 0.5*(zi-mu)^2/s2  - logsqrt2pi
    lf2i <- lf2i + log(1-pnorm(res0)) - log(1-pnorm(res1))
    llik <- log( exp(lf1i) + exp(lf2i) )
    negll <- -sum(llik)

    return(negll)
  }

  # THIS MAXIMIZES THE Neglik_f LIKELIHOOD FUNCTION AND RETURNS AN MLE OBJECT
  Get_mle <- function(data) {
    # TODO: PROPER ERROR HANDLING
    tryCatch({
      mle <- bbmle::mle2(data = list(data = data),
                         optimizer = "nlminb",
                         lower = c(gamma = -Inf, sigma = 0, p0 = 0, mu = -Inf),
                         upper = c(gamma = Inf, sigma = Inf, p0 = 1, mu = Inf),
                         minuslogl = Neglik_f,
                         control = list(iter.max = 1000L),
                         start = start)
    }, error = function(e) {
      stop(paste("Sorry, seems like I encountered an error...:", e), call. = FALSE)
    })
    return(mle)
  }

  ##============================================================================##
  ## MAIN PROGRAM
  ##============================================================================##

  # combine errors
  if (log) {
    if (invert) {
      lcd <- log(data[ ,1])*-1
      x.offset <- abs(min(lcd))
      lcd <- lcd+x.offset
    } else {
      lcd <- log(data[ ,1])
    }
    lse <- sqrt((data[ ,2]/data[ ,1])^2 + sigmab^2)
  } else {
    lcd <- data[ ,1]
    lse <- sqrt(data[ ,2]^2 + sigmab^2)
  }

  # create new data frame with DE and combined relative error
  dat <- cbind(lcd, lse)

  # get the maximum likelihood estimate
  ests <- Get_mle(dat)

  # check if any standard errors are NA or NaN
  coef_err <- t(as.data.frame(summary(ests)@coef[ ,2]))

  if (debug)
    print(summary(ests))

  if (any(is.nan(coef_err)))
    coef_err[which(is.nan(coef_err))] <- t(as.data.frame(ests@coef))/100
  if (any(is.na(coef_err)))
    coef_err[which(is.na(coef_err))] <- t(as.data.frame(ests@coef))/100

  if (par == 3)
    which <- c("gamma", "sigma", "p0")
  if (par == 4)
    which <- c("gamma", "sigma", "p0", "mu")

  # calculate profile log likelihoods
  prof <- bbmle::profile(ests,
                         which = which,
                         std.err = as.vector(coef_err),
                         #try_harder = TRUE,
                         quietly = TRUE,
                         tol.newmin = Inf,
                         skiperrs = TRUE,
                         prof.lower=c(gamma = -Inf, sigma = 0, p0 = 0, mu = -Inf),
                         prof.upper=c(gamma = Inf, sigma = Inf, p0 = 1, mu = Inf)
  )
  # Fallback when profile() returns a 'better' fit
  maxsteps <- 100
  cnt <- 1
  while (!inherits(prof, "profile.mle2")) {
    message(paste0("## Trying to find a better fit (", cnt, "/10) ##"))
    if (maxsteps == 0L)
      stop(paste("Sorry, but I can't find a converging fit for the profile log-likelihood."),
           call. = FALSE)

    prof <- profile(ests,
                    which = which,
                    std.err = as.vector(coef_err),
                    try_harder = TRUE,
                    quietly = TRUE,
                    maxsteps = maxsteps,
                    tol.newmin = Inf,
                    skiperrs = TRUE,
                    prof.lower = xlb,
                    prof.upper = xub
    )
    maxsteps <- maxsteps - 10
    cnt <- cnt + 1
  }

  ## TODO: reduce the redundant code
  ## DELETE rows where z = -Inf/Inf
  prof@profile$gamma <-  prof@profile$gamma[which(prof@profile$gamma["z"] != Inf), ]
  prof@profile$gamma <-  prof@profile$gamma[which(prof@profile$gamma["z"] != -Inf), ]
  prof@profile$sigma <-  prof@profile$sigma[which(prof@profile$sigma["z"] != Inf), ]
  prof@profile$sigma <-  prof@profile$sigma[which(prof@profile$sigma["z"] != -Inf), ]
  prof@profile$p0 <-  prof@profile$p0[which(prof@profile$p0["z"] != Inf), ]
  prof@profile$p0 <-  prof@profile$p0[which(prof@profile$p0["z"] != -Inf), ]

  if (par == 4) {
    prof@profile$mu <-  prof@profile$mu[which(prof@profile$mu["z"] != Inf), ]
    prof@profile$mu <-  prof@profile$mu[which(prof@profile$mu["z"] != -Inf), ]
  }

  # calculate Bayesian Information Criterion (BIC)
  BIC <- BIC(ests)

  # retrieve results from mle2-object
  pal <- if (log) {
    if (invert) {
      exp((bbmle::coef(ests)[["gamma"]]-x.offset)*-1)
    } else {
      exp(bbmle::coef(ests)[["gamma"]])
    }
  } else {
    bbmle::coef(ests)[["gamma"]]
  }
  sig <- bbmle::coef(ests)[["sigma"]]
  p0end <- bbmle::coef(ests)[["p0"]]

  if (par == 4) {
    muend <- ifelse(log, exp(bbmle::coef(ests)[["mu"]]), bbmle::coef(ests)[["mu"]])
  } else {
    muend <- NA
  }

  ##============================================================================##
  ## ERROR CALCULATION

  #### METHOD 1: follow the instructions of Galbraith & Roberts (2012) ####
  # "If the likelihood profile is symmetrical about the parameter, an approximate standard error
  #  can be calculated by dividing the length of this interval by 3.92"
  conf <- as.data.frame(bbmle::confint(prof, tol.newmin = Inf, quietly = TRUE))

  if (invert) {
    conf[1, ] <- (conf[1, ]-x.offset)*-1
    t <- conf[1,1]
    conf[1,1] <- conf[1,2]
    conf[1,2] <- t
  }
  gamma_err <- if (log) {
    (exp(conf["gamma",2])-exp(conf["gamma",1]))/3.92
  } else {
    (conf["gamma",2]-conf["gamma",1])/3.92
  }

  ##============================================================================##
  ## AGGREGATE RESULTS
  summary <- data.frame(de=pal,
                        de_err=gamma_err,
                        "ci_lower"=ifelse(log, exp(conf["gamma",1]), conf["gamma",1]),
                        "ci_upper"=ifelse(log, exp(conf["gamma",2]), conf["gamma",2]),
                        par=par,
                        sig=sig,
                        p0=p0end,
                        mu=muend,
                        Lmax=-ests@min,
                        BIC=BIC)
  call <- sys.call()
  args <- list(log=log, sigmab=sigmab, bootstrap=bootstrap,
               init.values=start,
               bs.M=M, bs.N=N, bs.h=h, sigmab.sd=sigmab.sd)

  ##============================================================================##
  ## BOOTSTRAP
  ##============================================================================##
  if (bootstrap) {

    ## BOOTSTRAP FUNCTIONS ----
    # Function that draws N+M sets of integer values from 1:n and returns
    # both the indices and frequencies
    draw_Freq <- function() {
      f <- R <- matrix(0L, N+M, n)
      for (i in seq_len(N+M)) {
        R[i, ] <- sample(x = n, size = n, replace = TRUE)
        f[i, ] <- tabulate(R, n)
      }
      return(list(R = R, freq = f))
    }

    # Function that adds the additional error sigmab to each individual DE error
    combine_Errors <- function(d, e) {
      if (log) {
        d[ ,2] <- sqrt((d[ ,2]/d[ ,1])^2 + e^2)
        d[ ,1] <- log(d[ ,1])
      } else {
        d[ ,2] <- sqrt(d[ ,2]^2 + e^2)
      }
      return(d)
    }

    # Function that produces N+M replicates from the original data set using
    # randomly sampled indices with replacement and adding a randomly drawn
    # sigmab error
    create_Replicates <- function(f, s) {
      d <- apply(f$R, 1, function(x) data[x, ])
      r <- mapply(function(x, y) combine_Errors(x, y), d, s, SIMPLIFY = FALSE)
      return(r)
    }

    # Function to extract the estimate of gamma from mle2 objects and converting
    # it back to the 'normal' scale
    save_Gamma <- function(d) {
      if (log) {
        if (invert) {
          m <- exp((bbmle::coef(d)[["gamma"]]-x.offset)*-1)
        } else {
          m <- exp(bbmle::coef(d)[["gamma"]])
        }
      } else {
        m <- bbmle::coef(d)[["gamma"]]
      }
      return(m)
    }

    # Function that takes each of the N replicates and produces a kernel density
    # estimate of length n. The normalised values are then returned as a matrix
    # with dimensions [N, n]
    get_KDE <- function(d) {
      f <- approx(density(x=d[ ,1], kernel="gaussian", bw = h), xout = d[ ,1])
      pStarTheta <- as.vector(f$y / sum(f$y))
      x <- matrix(t(pStarTheta/(1/n)), N, n, byrow = TRUE)
      return(x)
    }

    # Function that calculates the product term of the recycled bootstrap
    get_ProductTerm <- function(Pmat, b2Pmatrix) {
      prodterm <- apply(Pmat^b2Pmatrix$freq[1:N, ], 1, prod)
      return(prodterm)
    }

    # Function that calculates the pseudo likelihoods for M replicates and
    # returns the dose-likelihood pairs
    make_Pairs <- function(theta, b2mamvec, prodterm) {
      pairs <- matrix(0, M, 2)
      for (i in seq_len(M)) {
        thetavec <- matrix(theta[i], N, 1)
        kdthis <- (thetavec-b2mamvec)/h
        kd1 <- dnorm(kdthis)

        kd2 <- kd1*prodterm[[i]]
        kd <- sum(kd2)
        likelihood <- (1/(N*h))*kd
        pairs[i, ] <- c(theta[i], likelihood)
      }
      return(pairs)
    }

    ## START BOOTSTRAP ----
    msg <- sprintf(paste("\n [calc_MinDose] \n\nRecycled Bootstrap",
                         "\n\nParameters:",
                         "\n M = %d",
                         "\n N = %d",
                         "\n sigmab = %.2f \U00B1 %.2f",
                         "\n h = %.2f",
                         "\n\n Creating %d bootstrap replicates..."),
                   M, N, sigmab, sigmab.sd, h, N+M)
    message(msg)

    n <- length(data[ ,1])
    # Draw N+M samples of a normale distributed sigmab
    sigmab <- rnorm(N + M, sigmab, sigmab.sd)
    # Draw N+M random indices and their frequencies
    b2Pmatrix <- draw_Freq()
    # Finally draw N+M bootstrap replicates
    replicates <- create_Replicates(b2Pmatrix, sigmab)

    # MULTICORE: The call to 'Get_mle' is the bottleneck of the function.
    # Using multiple CPU cores can reduce the computation cost, but may
    # not work for all machines.
    if (multicore) {
      message(paste("\n Spawning", cores, "instances of R for parallel computation. This may take a few seconds..."))
      cl <- parallel::makeCluster(cores)
      message("\n Done! Applying the model to all replicates. This may take a while...")
      mle <- parallel::parLapply(cl, replicates, Get_mle)
      parallel::stopCluster(cl)
    } else {
      message("\n Applying the model to all replicates. This may take a while...")
      mle <- lapply(replicates, Get_mle)
    }

    # Final bootstrap calculations
    message("\n Calculating the likelihoods...")
    # Save 2nd- and 1st-level bootstrap results (i.e. estimates of gamma)
    b2mamvec <- as.matrix(sapply(mle[1:N], save_Gamma, simplify = TRUE))
    theta <- sapply(mle[c(N+1):c(N+M)], save_Gamma)
    # Calculate the probality/pseudo-likelihood
    Pmat <- lapply(replicates[c(N+1):c(N+M)], get_KDE)
    prodterm <- lapply(Pmat, get_ProductTerm, b2Pmatrix)
    # Save the bootstrap results as dose-likelihood pairs
    pairs <- make_Pairs(theta, b2mamvec, prodterm)

    ## --------- FIT POLYNOMIALS -------------- ##
    message("\n Fit curves to dose-likelihood pairs...")
    # polynomial fits of increasing degrees
    poly.three <- lm(pairs[ ,2] ~ poly(pairs[ ,1], degree = 3, raw = TRUE))
    poly.four <- lm(pairs[ ,2] ~ poly(pairs[ ,1], degree = 4, raw = TRUE))
    poly.five <- lm(pairs[ ,2] ~ poly(pairs[ ,1], degree = 5, raw = TRUE))
    poly.six <- lm(pairs[ ,2] ~ poly(pairs[ ,1], degree = 6, raw = TRUE))

    ## --------- FIT LOESS -------------- ##
    # Polynomials are probably not reasonable and often suffer badly from
    # overfitting, especially towards the margins of the fitted data. In this
    # particular use case polynomials may suggest a multimodal likelihood
    # distribution where actually none is given. The non-parametric
    # LOESS (LOcal polynomial regrESSion) often yields better results than
    # standard polynomials.
    loess <- loess(pairs[ ,2] ~ pairs[ ,1])

  }#EndOf::Bootstrap

  ##============================================================================##
  ## CONSOLE PRINT
  ##============================================================================##
  if (verbose) {
    if (!bootstrap) {
      cat("\n----------- meta data -----------\n")
      print(data.frame(n=length(data[ ,1]),
                       par=par,
                       sigmab=sigmab,
                       logged=log,
                       Lmax=-ests@min,
                       BIC=BIC,
                       row.names = ""))

      cat("\n--- final parameter estimates ---\n")
      print(round(data.frame(gamma=ifelse(!invert, bbmle::coef(ests)[["gamma"]], (bbmle::coef(ests)[["gamma"]]-x.offset)*-1),
                             sigma=bbmle::coef(ests)[["sigma"]],
                             p0=bbmle::coef(ests)[["p0"]],
                             mu=ifelse(par==4, ifelse(log,log(muend),muend),0),
                             row.names=""), 2))

      cat("\n------ confidence intervals -----\n")
      print(round(conf, 2))

      cat("\n------ De (asymmetric error) -----\n")
      print(round(data.frame(De=pal,
                             "lower"=ifelse(log, ifelse(!invert, exp(conf["gamma",1]), exp((conf["gamma",2]-x.offset)*-1)), conf["gamma",1]),
                             "upper"=ifelse(log, ifelse(!invert, exp(conf["gamma",2]), exp((conf["gamma",1]-x.offset)*-1)), conf["gamma",2]),
                             row.names=""), 2))

      cat("\n------ De (symmetric error) -----\n")
      print(round(data.frame(De=pal,
                             error=gamma_err,
                             row.names=""), 2))

    } else if (bootstrap) {
      message("\n Finished!")
    }
  }

  ##============================================================================##
  ## RETURN VALUES
  ##============================================================================##

  if (!bootstrap)
    pairs <- poly.three <- poly.four <- poly.five <- poly.six <- loess <- NULL

  newRLumResults.calc_MinDose <- set_RLum(
    class = "RLum.Results",
    data = list(summary = summary,
                data = data,
                args = args,
                call = call,
                mle = ests,
                BIC = BIC,
                confint = conf,
                profile = prof,
                bootstrap = list(
                  pairs = list(gamma=pairs),
                  poly.fits = list(poly.three = poly.three,
                                   poly.four = poly.four,
                                   poly.five = poly.five,
                                   poly.six = poly.six),
                  loess.fit = loess)))

  ##=========##
  ## PLOTTING
  if (plot)
    try(plot_RLum.Results(newRLumResults.calc_MinDose, ...))


  if (!debug)
    options(warn = 0)

  invisible(newRLumResults.calc_MinDose)

}
