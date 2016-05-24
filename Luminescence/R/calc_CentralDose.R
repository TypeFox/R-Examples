#' Apply the central age model (CAM) after Galbraith et al. (1999) to a given
#' De distribution
#'
#' This function calculates the central dose and dispersion of the De
#' distribution, their standard errors and the profile log likelihood function
#' for sigma.
#'
#' This function uses the equations of Galbraith & Roberts (2012). The
#' parameters \code{delta} and \code{sigma} are estimated by numerically solving
#' eq. 15 and 16. Their standard errors are approximated using eq. 17.
#' In addition, the profile log-likelihood function for \code{sigma} is
#' calculated using eq. 18 and presented as a plot. Numerical values of the 
#' maximum likelihood approach are \bold{only} presented in the plot and \bold{not}
#' in the console. A detailed explanation on maximum likelihood estimation can be found in the
#' appendix of Galbraith & Laslett (1993, 468-470) and Galbraith & Roberts
#' (2012, 15)
#'
#' @param data \code{\linkS4class{RLum.Results}} or \link{data.frame}
#' (\bold{required}): for \code{data.frame}: two columns with De
#' \code{(data[,1])} and De error \code{(values[,2])}
#' 
#' @param sigmab \code{\link{numeric}} (with default): spread in De values
#' given as a fraction (e.g. 0.2). This value represents the expected
#' overdispersion in the data should the sample be well-bleached (Cunningham &
#' Walling 2012, p. 100).
#' 
#' @param log \code{\link{logical}} (with default): fit the (un-)logged central
#' age model to De data
#' 
#' @param plot \code{\link{logical}} (with default): plot output
#' 
#' @param \dots further arguments (\code{trace, verbose}).
#' 
#' @return Returns a plot (optional) and terminal output. In addition an
#' \code{\linkS4class{RLum.Results}} object is returned containing the
#' following element:
#'
#' \item{summary}{\link{data.frame} summary of all relevant model results.}
#' \item{data}{\link{data.frame} original input data} \item{args}{\link{list}
#' used arguments} \item{call}{\link{call} the function call}
#' \item{profile}{\link{data.frame} the log likelihood profile for sigma}
#'
#' The output should be accessed using the function
#' \code{\link{get_RLum}}
#' @section Function version: 1.3.1
#' @author Christoph Burow, University of Cologne (Germany) \cr Based on a
#' rewritten S script of Rex Galbraith, 2010 \cr
#' @seealso \code{\link{plot}}, \code{\link{calc_CommonDose}},
#' \code{\link{calc_FiniteMixture}}, \code{\link{calc_FuchsLang2001}},
#' \code{\link{calc_MinDose}}
#' @references Galbraith, R.F. & Laslett, G.M., 1993. Statistical models for
#' mixed fission track ages. Nuclear Tracks Radiation Measurements 4, 459-470.
#' \cr \cr Galbraith, R.F., Roberts, R.G., Laslett, G.M., Yoshida, H. & Olley,
#' J.M., 1999. Optical dating of single grains of quartz from Jinmium rock
#' shelter, northern Australia. Part I: experimental design and statistical
#' models.  Archaeometry 41, 339-364. \cr \cr Galbraith, R.F. & Roberts, R.G.,
#' 2012. Statistical aspects of equivalent dose and error calculation and
#' display in OSL dating: An overview and some recommendations. Quaternary
#' Geochronology 11, 1-27. \cr \cr \bold{Further reading} \cr \cr Arnold, L.J.
#' & Roberts, R.G., 2009. Stochastic modelling of multi-grain equivalent dose
#' (De) distributions: Implications for OSL dating of sediment mixtures.
#' Quaternary Geochronology 4, 204-230. \cr \cr Bailey, R.M. & Arnold, L.J.,
#' 2006. Statistical modelling of single grain quartz De distributions and an
#' assessment of procedures for estimating burial dose. Quaternary Science
#' Reviews 25, 2475-2502. \cr \cr Cunningham, A.C. & Wallinga, J., 2012.
#' Realizing the potential of fluvial archives using robust OSL chronologies.
#' Quaternary Geochronology 12, 98-106. \cr \cr Rodnight, H., Duller, G.A.T.,
#' Wintle, A.G. & Tooth, S., 2006. Assessing the reproducibility and accuracy
#' of optical dating of fluvial deposits.  Quaternary Geochronology, 1 109-120.
#' \cr \cr Rodnight, H., 2008. How many equivalent dose values are needed to
#' obtain a reproducible distribution?. Ancient TL 26, 3-10.
#' @examples
#'
#' ##load example data
#' data(ExampleData.DeValues, envir = environment())
#'
#' ##apply the central dose model
#' calc_CentralDose(ExampleData.DeValues$CA1)
#'
#' @export
calc_CentralDose <- function(data, sigmab, log = TRUE, plot = TRUE, ...) {
  
  ## ============================================================================##
  ## CONSISTENCY CHECK OF INPUT DATA
  ## ============================================================================##
  
  if (!missing(data)) {
    if (!is(data, "data.frame") & !is(data, "RLum.Results")) {
      stop("[calc_CentralDose] Error: 'data' object has to be of type\n 
           'data.frame' or 'RLum.Results'!")
    } else {
      if (is(data, "RLum.Results")) {
        data <- get_RLum(data, signature(object = "De.values"))
      }
    }
  }
  try(colnames(data) <- c("ED", "ED_Error"), silent = TRUE)
  if (colnames(data[1]) != "ED" || colnames(data[2]) != "ED_Error") {
    cat(paste("Columns must be named 'ED' and 'ED_Error'"), fill = FALSE)
    stop(domain = NA)
  }
  if (!missing(sigmab)) {
    if (sigmab < 0 | sigmab > 1 & log) {
      cat(paste("sigmab needs to be given as a fraction between", "0 and 1 (e.g. 0.2)"), 
          fill = FALSE)
      stop(domain = NA)
    }
  }
  
  ## ============================================================================##
  ## ... ARGUMENTS
  ## ============================================================================##
  
  options <- list(verbose = TRUE,
                  trace = FALSE)
  
  options <- modifyList(options, list(...))
  
  
  ## ============================================================================##
  ## CALCULATIONS
  ## ============================================================================##
  
  # set default value of sigmab
  if (missing(sigmab)) 
    sigmab <- 0
  
  # calculate yu = log(ED) and su = se(logED)
  if (log) {
    yu <- log(data$ED)
    su <- sqrt((data$ED_Error / data$ED)^2 + sigmab^2)
  } else {
    yu <- data$ED
    su<- sqrt((data$ED_Error)^2 + sigmab^2)
  }

  
  # calculate starting values and weights
  sigma <- 0.15
  wu <- 1 / (sigma^2 + su^2)
  delta <- sum(wu * yu) / sum(wu)
  n <- length(yu)
  
  # compute mle's
  for (j in 1:200) {
    delta <- sum(wu * yu) / sum(wu)
    sigma <- sigma * sqrt(sum((wu^2) * (yu - delta)^2 / sum(wu)))
    wu <- 1 / (sigma^2 + su^2)
    
    # print iterations
    if (options$trace) 
      print(round(c(delta, sigma), 4))
  }
  
  # save parameters for terminal output
  out.delta <- ifelse(log, exp(delta), delta)
  out.sigma <- ifelse(log, sigma * 100, sigma / out.delta * 100)
  
  # log likelihood
  llik <- 0.5 * sum(log(wu)) - 0.5 * sum(wu * (yu - delta)^2)
  # save parameter for terminal output
  out.llik <- round(llik, 4)
  Lmax <- llik
  
  # standard errors
  sedelta <- 1 / sqrt(sum(wu))
  sesigma <- 1 / sqrt(2 * sigma^2 * sum(wu^2))
  
  # save parameters for terminal output
  if (log) {
    out.sedelta <- sedelta * 100
    out.sesigma <- sesigma
  } else {
    out.sedelta <- sedelta / out.delta * 100
    out.sesigma <- sqrt((sedelta / delta)^2 + 
                          (sesigma / out.delta * 100 / out.sigma)^2) * out.sigma / 100
    
  }
  
  # profile log likelihood
  sigmax <- sigma
  llik <- 0
  sig0 <- max(0, sigmax - 8 * sesigma)
  sig1 <- sigmax + 9.5 * sesigma
  sig <- try(seq(sig0, sig1, sig1 / 1000), silent = TRUE)
  
  if (class(sig) != "try-error") {
    # TODO: rewrite this loop as a function and maximise with mle2 ll is the actual
    # log likelihood, llik is a vector of all ll
    for (s in sig) {
      wu <- 1 / (s^2 + su^2)
      mu <- sum(wu * yu)/sum(wu)
      ll <- 0.5 * sum(log(wu)) - 0.5 * sum(wu * (yu - mu)^2)
      llik <- c(llik, ll)
    }
    llik <- llik[-1] - Lmax
  } #endif::try-error
  
  
  ## ============================================================================##
  ## TERMINAL OUTPUT
  ## ============================================================================##
  
  if (options$verbose) {
    cat("\n [calc_CentralDose]")
    cat(paste("\n\n----------- meta data ----------------"))
    cat(paste("\n n:                      ", n))
    cat(paste("\n log:                    ", log))
    cat(paste("\n----------- dose estimate ------------"))
    cat(paste("\n central dose [Gy]:      ", format(out.delta, digits = 2, nsmall = 2)))
    cat(paste("\n SE [Gy]:                ", format(out.delta * out.sedelta/100, 
                                                   digits = 2, nsmall = 2)))
    cat(paste("\n rel. SE [%]:            ", format(out.sedelta, digits = 2, nsmall = 2)))
    cat(paste("\n----------- overdispersion -----------"))
    cat(paste("\n OD [Gy]:                ", format(ifelse(log, sigma * out.delta, sigma), digits = 2, nsmall = 2)))
    cat(paste("\n SE [Gy]:                ", format(ifelse(log, sesigma * out.delta, sesigma), digits = 2, nsmall = 2)))
    cat(paste("\n OD [%]:                 ", format(out.sigma, digits = 2, nsmall = 2)))
    cat(paste("\n SE [%]:                 ", if (class(sig) != "try-error") {
      format(out.sesigma * 100, digits = 2, nsmall = 2)
    } else {
      "-"
    }))
    cat(paste("\n-------------------------------------\n\n"))
  }
  
  ## ============================================================================##
  ## RETURN VALUES
  ## ============================================================================##
  
  if (class(sig) == "try-error") {
    out.sigma <- 0
    out.sesigma <- NA
  }
  
  if(!log)
    sig <- sig / delta
    
  
  summary <- data.frame(de = out.delta, de_err = out.delta * out.sedelta / 100, 
                        OD = out.sigma, OD_err = out.sesigma * 100, Lmax = Lmax)
  
  call <- sys.call()
  args <- list(log = "TRUE", sigmab = sigmab)
  
  newRLumResults.calc_CentralDose <- set_RLum(class = "RLum.Results", 
                                              data = list(summary = summary, 
                                                          data = data, 
                                                          args = args,
                                                          call = call, 
                                                          profile = data.frame(sig = sig, 
                                                                               llik = llik)))
  
  ## =========## PLOTTING
  if (plot && class(sig) != "try-error") 
    try(plot_RLum.Results(newRLumResults.calc_CentralDose, ...))

  invisible(newRLumResults.calc_CentralDose)
}

