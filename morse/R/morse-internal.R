# Ugly hack to get rid of spurious notes in package check, caused by uses
# of dplyr::{rename, filter}. R is such a sad language.
if(getRversion() >= "2.15.1")  utils::globalVariables(c(
  "response", "Nreprocumul", "resp", "Mortality", "qinf95", "qsup95",
  "transf_conc", "obs", "pred", "..n..", "Points", "conc", "Line", "Nsurv",
  "time", "Conf.Int", "Cred.Lim", "Obs", "P50", "P2.5", "P97.5"
))


# Generates a character string vector from a data.frame using its replicate,
# conc and time columns. The result can be used as identifiers for the rows
# of the data.set.
#
#' @importFrom stringr str_c
#'
idCreate <- function(data) {
  str_c(data[, "replicate"],
        data[, "conc"],
        data[, "time"],
        sep = "_")
}

#' @importFrom dplyr filter
selectDataTT <- function(data, target.time) {
  # INPUT
  # - data: An object of class reproData or survData
  # - target.time: the time we want to consider as the last time for the analysis.
  # OUTPUT
  # - subset the dataframe for target.time used by function with target.time arg

  # target.time default
  if ("time" %in% colnames(data)) { # one time dataset check
    if (is.null(target.time)) {
      target.time <- max(data$time)
    }

    # correct target time
    if (!any(data$time == target.time))
      stop("target.time is not one of the possible time !")

    dataTT <- filter(data, time == target.time)
  } else {
    dataTT <- cbind(data, time = 1)
  }

  return(dataTT)
}

#' @import rjags
#' @importFrom coda raftery.diag
modelSamplingParameters <- function(model, parameters, n.chains, quiet = quiet) {
  # estimate the number of iteration required for the estimation
  # by using the raftery.diag
  # INPUTS:
  # - model: jags model from loading function
  # - parameters: parameters from loading function
  # - nchains: Number of chains desired
  # - quiet: silent option
  # OUTPUTS:
  # - niter: number of iteration (mcmc)
  # - thin: thining rate parameter
  # - burnin: number of iteration burned

  # number of iteration for the pilote run required by raftery.diag
  # default value: 3746
  niter.init <- 5000
  prog.b <- ifelse(quiet == TRUE, "none", "text") # plot progress bar option
  mcmc <- coda.samples(model, parameters, n.iter = niter.init, thin = 1,
                       progress.bar = prog.b)
  RL <- raftery.diag(mcmc)

  # check raftery.diag result (number of sample for the diagnostic procedure)
  if (n.chains < 2) stop('2 or more parallel chains required !')

  # extract raftery diagnostic results
  resmatrix <- RL[[1]]$resmatrix
  for (i in 2: length(RL)) {
    resmatrix <- rbind(resmatrix, RL[[i]]$resmatrix)
  }

  # creation of sampling parameters
  thin <- round(max(resmatrix[, "I"]) + 0.5) # autocorrelation
  niter <- max(resmatrix[, "Nmin"]) * thin # number of iteration
  burnin <- max(resmatrix[, "M"]) # burnin period

  return(list(niter = niter, thin = thin, burnin = burnin))
}

#' @import rjags
calcDIC <- function(m.M, sampling.parameters, quiet = quiet) {
  # calculate the dic for a jags model
  # INPUTS
  # - m.M:  jags model object
  # - niter: number of iterations for the sampling
  # - thin
  # OUTPUT:
  # - numeric value of the DIC

  prog.b <- ifelse(quiet == TRUE, "none", "text") # plot progress bar option

  # estimation of DIC
  dic <- dic.samples(m.M, n.iter = sampling.parameters$niter,
                     thin = sampling.parameters$thin, progress.bar = prog.b)

  # return penalised DIC
  return(round(sum(sapply(dic$deviance, mean) + sapply(dic$penalty, mean))))
}

#' @importFrom stats quantile
estimXCX <- function(mcmc, xcx, varx) {
  # create the table of estimated values of LCx or ECx
  # for the survival analyses

  # INPUT:
  # - mcmc:  list of estimated parameters for the model with each item representing
  # a chains
  # - xcx: vector of values of LCx or ECX
  # - varx: character string for lcx or ecx
  # OUTPUT:
  # - data frame with the estimated ECx and their CIs 95% (3 columns (values,
  # CIinf, CIsup) and length(x) rows)

  # Retrieving estimated parameters of the model
  mctot <- do.call("rbind", mcmc)
  b <- 10^mctot[, "log10b"]
  e <- 10^mctot[, "log10e"]

  # Calculation XCx median and quantiles
  XCx <- sapply(xcx, function(x) {e * ((100 / (100 - x)) - 1)^(1 / b)})

  q50 <- apply(XCx, 2, function(x) {quantile(x, probs = 0.5)})
  qinf95 <- apply(XCx, 2, function(x) {quantile(x, probs = 0.025)})
  qsup95 <- apply(XCx, 2, function(x) {quantile(x, probs = 0.975)})

  # defining names
  XCname <- sapply(xcx, function(x) {paste(varx, x, sep = '')})

  # create the dataframe with ECx median and quantiles
  res <- data.frame(median = q50, Q2.5 = qinf95, Q97.5 = qsup95,
                    row.names = XCname)

  return(res)
}

optLogTransform <- function(log.scale, x) {
  if(log.scale) log(x) else x
}

legendGgplotFit <- function(a.gplot) {
  # create multi legend for plot.reproFitTt and plot.survFitTt
  # to have differente legend for points and mean and CI curve
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

fCols <- function(data, fitcol, cicol, analyse) {
  
  if (analyse == "surv") {
    #points
    cols1 <- "black"
    names(cols1) <- unique(data$Points)
    # curve
    cols2 <- fitcol
    names(cols2) <- "loglogistic"
    # CI
    cols3 <- "black"
    names(cols3) <- "Confidence interval"
    # CI2
    cols4 <- cicol
    names(cols4) <- "Credible limits"
    
    return(list(cols1 = cols1,
                cols2 = cols2,
                cols3 = cols3,
                cols4 = cols4))
           
    } else if (analyse == "repro") {
    # fitted curve
    cols2 <- fitcol
    names(cols2) <- "loglogistic"
    # CI
    cols3 <- cicol
    names(cols3) <- "Credible limits"
    
    return(list(cols2 = cols2,
                cols3 = cols3))
    }
}

exclude_labels <- function(x) {
  x[-seq.int(1, length(x), 4)] <- ""
  return(x)
}
