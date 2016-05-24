#' @name Variability_Measures
#' @title Variability Measures
#' @rdname VarMeasures
#' @aliases by_id
#' @aliases sd_id
#' @aliases rmssd
#' @aliases rmssid_id
#' @aliases rolling_diff
#' @aliases rolling_diff_id
#'
#' @note These are a set of functions designed to calculate various
#' measures of variability either on a single data vector, or
#' calculate them by an ID.
#'
#' @param x A data vector to operate on.  Should be a numeric or
#'   integer vector, or coercible to such (e.g., logical).
#' @param ID an ID variable indicating how to split up the \code{x}
#'   vector.  Should be the same length as \code{x}.
#' @param fun The function to calculate by ID
#' @param long A logical indicating whether to return results in
#'   \dQuote{long} form (the default) or wide (if \code{FALSE}).
#' @param \dots Additional arguments passed on to \code{fun}
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
NULL


#' Variability Measures
#'
#' \code{by_id} - Internal function to allow a simple statistic (e.g., SD)
#' to be calculated individually by an ID variable and returned
#' either as per ID (i.e., wide form) or for every observation of an
#' ID (i.e., long form).
#' @return \code{by_id} - A vector the same length as \code{x}
#'   if \code{long=TRUE}, or the length of unique \code{ID}s if
#'   \code{long=FALSE}.
#' @rdname VarMeasures
by_id <- function(x, ID, fun, long=TRUE, ...) {
  if (long) {
    ave(x, ID, FUN = function(x) fun(x, ...))
  } else {
    tapply(x, ID, FUN = function(x) fun(x, ...))
  }
}

#' Variability Measures
#'
#' \code{sd_id} - Calculates the standard deviation of observations by \code{ID}.
#'
#' @return \code{sd_id} - A vector of the standard deviations by ID
#' @keywords utilities
#' @export
#' @rdname VarMeasures
#' @examples
#'
#' sd_id(mtcars$mpg, mtcars$cyl, long=TRUE)
#' sd_id(mtcars$mpg, mtcars$cyl, long=FALSE)
sd_id <- function(x, ID, long=TRUE) {
  by_id(x, ID, fun = sd, long = long, na.rm=TRUE)
}

#' Variability Measures
#'
#' \code{rmssd} - Calculates the root mean square of successive differences (RMSSD).
#'   Note that missing values are removed.
#'
#' @return \code{rmssd} - The RMSSD for the data.
#' @export
#' @rdname VarMeasures
#' @examples
#' rmssd(1:4)
#' rmssd(c(1, 3, 2, 4))
rmssd <- function(x) {
  x <- na.omit(diff(x))
  as.vector(sqrt(mean(x^2)))
}

#' Variability Measures
#'
#' \code{rmssd_id} - Calculates the RMSSD by ID.
#'
#' @return \code{rmssd_id} - A vector of the RMSSDs by ID
#' @export
#' @rdname VarMeasures
#' @examples
#' rmssd_id(mtcars$mpg, mtcars$cyl)
#' rmssd_id(mtcars$mpg, mtcars$cyl, long=FALSE)
rmssd_id <- function(x, ID, long=TRUE) {
  by_id(x, ID, fun = rmssd, long = long)
}

#' Variability Measures
#'
#' \code{rolling_diff} - Calculates the average rolling difference of the data.
#'   Within each window, the difference between the maximum and minimum value is
#'   computed and these are averaged across all windows.  The equation is:
#' \deqn{\frac{\sum_{t = 1}^{N - k} max(x_{t}, \ldots, x_{t + k}) - min(x_{t}, \ldots, x_{t + k})}{N - k}}
#'
#' @param window An integer indicating the size of the rolling window.
#'   Must be at least the length of \code{x}.
#' @return \code{rolling_diff} - The average of the rolling differences between maximum and minimum.
#' @export
#' @rdname VarMeasures
#' @examples
#' rolling_diff(1:7, window = 4)
#' rolling_diff(c(1, 4, 3, 4, 5))
rolling_diff <- function(x, window = 4) {
  stopifnot(length(x) >= window)

  index <- 1:(length(x) + 1 - window)

  mean(sapply(index, function(i) {
    x <- na.omit(x[i:(i + window - 1)])
    if (length(x) < 2) {
      NA
    } else {
      diff(range(x))
    }
  }), na.rm=TRUE)
}

#' Variability Measures
#'
#' \code{rolling_diff_id} - Calculates the average rolling difference by ID
#'
#' @return \code{rolling_diff_id} - A vector of the average rolling differences by ID
#' @export
#' @rdname VarMeasures
#' @examples
#' rolling_diff_id(mtcars$mpg, mtcars$cyl, window = 3)
rolling_diff_id <- function(x, ID, long=TRUE, window = 4) {
  by_id(x, ID, fun = rolling_diff, long = long, window = window)
}


#' Estimate the parameters for a Gamma distribution
#'
#' This is a simple function to estimate what the parameters for a Gamma
#' distribution would be from a data vector.  It is used internally to
#' generate start values.
#'
#' @param x a data vector to operate on
#' @return a list of the shape (alpha) and rate (beta) parameters
#'   and the mean and variance
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @keywords utilities
gamma_params <- function(x) {
  m <- mean(x, na.rm=TRUE)
  v <- var(x, na.rm=TRUE)
  beta <- m / v
  alpha <- m * beta
  list(alpha = alpha, beta = beta, mean = m, variance = v)
}

#' Estimates the parameters of a Gamma distribution from SDs
#'
#' This function calcualtes the parameters of a Gamma distribution
#' from the residuals from an individuals' own mean.
#' That is, the distribution of (standard) deviations from individuals'
#' own mean are calculated and then an estimate of the parameters of a
#' Gamma distribution are calculated.
#'
#' @param x A data vector to operate on
#' @param ID an ID variable of the same length as \code{x}
#' @return a list of the shape (alpha) and rate (beta) parameters
#'   and the mean and variance
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @export
#' @keywords utilities
#' @examples
#'
#' set.seed(1234)
#' y <- rgamma(100, 3, 2)
#' x <- rnorm(100 * 10, mean = 0, sd = rep(y, each = 10))
#' ID <- rep(1:100, each = 10)
#' res_gamma(x, ID)
res_gamma <- function(x, ID) {
  gamma_params(sd_id(x, ID, long = FALSE))
}

#' Calculates an empirical p-value based on the data
#'
#' This function takes a vector of statistics and calculates
#' the empirical p-value, that is, how many fall on the other
#' side of zero.  It calculates a two-tailed p-value.
#'
#' @param x a data vector to operate on
#' @param na.rm Logical whether to remove NA values. Defaults to \code{TRUE}
#' @return a named vector with the number of values falling at
#'   or below zero, above zero, and the empirical p-value.
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @export
#' @keywords utilities
#' @examples
#'
#' empirical_pvalue(rnorm(100))
empirical_pvalue <- function(x, na.rm = TRUE) {
  x <- as.integer(x <= 0)
  tmp <- table(factor(x, levels = 1:0, labels = c("<= 0", "> 0")))
  m <- mean(x, na.rm = na.rm)
  pval2 <- 2 * min(m, 1 - m)
  out <- c(as.vector(tmp), pval2)
  names(out) <- c(names(tmp), "p-value")

  out
}


#' nice formatting for p-values
#'
#' @param p a numeric pvalue
#' @param d the digits less than which should be displayed as less than
#' @param sd scientific digits for round
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @keywords utilities
#' @examples
#' varian:::pval_smartformat(c(1, .15346, .085463, .05673, .04837, .015353462,
#'   .0089, .00164, .0006589, .0000000053326), 3, 5)
pval_smartformat <- function(p, d = 3, sd = 5) {
  p.out <- ifelse(p < 1/(10^d),
         paste0("< .", paste(rep(0, d - 1), collapse = ""), "1"),
                  format(round(p, digits = d), digits = d, nsmall = d, scientific = sd))
  gsub("0\\.", ".", p.out)
}


#' Calculates summaries for a parameter
#'
#' This function takes a vector of statistics and calculates
#' several summaries: mean, median, 95% CI, and
#' the empirical p-value, that is, how many fall on the other
#' side of zero.
#'
#' @param x a data vector to operate on
#' @param digits Number of digits to round to for printing
#' @param pretty Logical value whether prettified values should be returned.
#'   Defaults to \code{FALSE}.
#' @param \dots Additional arguments passed to \code{pval_smartformat}
#'   to control p-value printing.
#' @param na.rm Logical whether to remove NA values. Defaults to \code{TRUE}
#' @return .
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @export
#' @keywords utilities
#' @examples
#'
#' param_summary(rnorm(100))
#' param_summary(rnorm(100), pretty = TRUE)
param_summary <- function(x, digits = 2, pretty = FALSE, ..., na.rm = TRUE) {
  res <- round(data.frame(Mean = mean(x, na.rm = na.rm),
    Median = median(x, na.rm = na.rm),
    SE = sd(x, na.rm = na.rm),
    LL2.5 = as.vector(quantile(x, probs = .025, na.rm = na.rm)),
    UL97.5 = as.vector(quantile(x, probs = .975, na.rm = na.rm))),
    digits = digits)

  p <- pval_smartformat(empirical_pvalue(x)[["p-value"]], ...)

  if (pretty) {
    out <- sprintf("%s [%s, %s], %s",
                   as.character(res$Mean),
                   as.character(res$LL2.5),
                   as.character(res$UL97.5),
                   ifelse(grepl("<", p), paste0("p ", p), paste0("p = ", p)))
  } else {
    res[, 'p-value'] <- p
    out <- res
  }

  return(out)
}


#' Simulate a Gamma Variability Model
#'
#' This function facilitates simulation of a Gamma Variability Model
#' and allows the number of units and repeated measures to be varied
#' as well as the degree of variability.
#'
#' @param n The number of repeated measures on each unit
#' @param k The number of units
#' @param mu The grand mean of the variable
#' @param mu.sigma The standard deviation of the random mean of the variable
#' @param sigma.shape the shape (alpha) parameter of the Gamma distribution
#'   controlling the residual variability
#' @param sigma.rate the rate (beta) parameter of the Gamma distribution
#'   controlling the residual variability
#' @param seed the random seed, used to make simulations reproductible.
#'   Defaults to 5346 (arbitrarily).
#' @return a list of the data, IDs, and the parameters used for the simulation
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @export
#' @import MASS
#' @keywords utilities
#' @examples
#' raw.sim <- simulate_gvm(12, 140, 0, 1, 4, .1, 94367)
#' sim.data <- with(raw.sim, {
#'   set.seed(265393)
#'   x2 <- MASS::mvrnorm(k, c(0, 0), matrix(c(1, .3, .3, 1), 2))
#'   y2 <- rnorm(k, cbind(Int = 1, x2) %*% matrix(c(3, .5, .7)) + sigma, sd = 3)
#'   data.frame(
#'     y = Data$y,
#'     y2 = y2[Data$ID2],
#'     x1 = x2[Data$ID2, 1],
#'     x2 = x2[Data$ID2, 2],
#'     ID = Data$ID2)
#' })
simulate_gvm <- function(n, k, mu, mu.sigma, sigma.shape, sigma.rate, seed = 5346) {
  set.seed(seed)
  m <- rnorm(k, mu, mu.sigma)
  sigma <- rgamma(k, sigma.shape, sigma.rate)
  y <- rnorm(n * k, mean = rep(m, each = n),
    sd = rep(sigma, each = n))
  y2 <- rnorm(k, sigma/sqrt(sigma.shape/(sigma.rate^2)), sd = 3)

  list(
    Data = data.frame(y = y, y2 = y2,
    ID1 = 1:(n * k), ID2 = rep(1:k, each = n)),
    n = n, k = k, mu = mu, mu.sigma,
    sigma.shape, sigma.rate, sigma = sigma, seed)
}

#' Wrapper for the stan function to parallelize chains
#'
#' This funcntion takes Stan model code, compiles the Stan model,
#' and then runs multiple chains in parallel.
#'
#' @param model_code A character string of Stan code
#' @param standata A data list suitable for Stan for the model given
#' @param totaliter The total number of iterations for inference.
#'   Note that the total number of iterations is automatically
#'   distributed across chains.
#' @param warmup How many warmup iterations should be used?  Note
#'   that every chain will use the same number of warmups and these
#'   will be \emph{added on top of the total iterations} for each chain.
#' @param thin The thin used, default to 1 indicating that all samples
#'   be saved.
#' @param chains The number of independent chains to run.
#' @param cl (optional) The name of a cluster to use to run the chains.
#'   If not specified, the function will make a new cluster.
#' @param cores (optional) If the \code{cl} argument is not used,
#'   this specifies the number of cores to make on the new cluster.
#'   If both \code{cl} and \code{cores} are missing, defaults to
#'   the minimum of the number of chains specified or the number of
#'   cores available on the machine.
#' @param seeds (optional) A vector of random seeds the same length as the number
#'   of independent chains being run, to make results replicable.
#'   If missing, random seeds will be generated and stored for reference
#'   in the output.
#' @param modelfit (optional) A compiled Stan model, if available, saves
#'   compiling \code{model_code}.
#' @param verbose A logical whether to print verbose output
#'   (defaults to \code{FALSE})
#' @param pars Parameter names from Stan to store
#' @param sample_file The sample file for Stan
#' @param diagnostic_file The diagnostic file for Stan
#' @param init A character string (\dQuote{random}) or a named list of starting values.
#' @param \dots Additional arguments, not currently used.
#' @return a named list with three elements, the \code{results},
#'   compiled Stan \code{model}, and the random \code{seeds}
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @export
#' @keywords utilities
#' @examples
#' # Make me!
parallel_stan <- function(model_code, standata, totaliter, warmup, thin = 1,
                  chains, cl, cores, seeds, modelfit, verbose = FALSE,
                  pars = NA, sample_file=NA, diagnostic_file = NA, init = "random", ...) {
  if (missing(cl)) {
    if (missing(cores)) {
      cores <- detectCores()
    }
    # make a cluster with nodes
    cl <- makeCluster(min(cores, chains))
  }
  on.exit(stopCluster(cl))

  if (!missing(modelfit)) {
    stanmodel <- modelfit
  } else {
    if (verbose) cat("Compiling Stan model\n")
    stanmodel <- stan_model(model_code = model_code, save_dso = TRUE)
  }

  if (verbose) cat("loading varian on workers\n")
  clusterEvalQ(cl, {
    require(varian)
  })

  if (missing(seeds)) {
    seeds <- sample(.Random.seed, chains)
  }

  eachiter <- ceiling(totaliter/chains)

  fooenv <- environment()

  clusterExport(cl, varlist = c(
      "stanmodel", "standata", "pars",
      "eachiter", "warmup", "thin", "seeds", "init",
      "sample_file", "diagnostic_file"),
    envir = fooenv)

  if (verbose) cat("Sampling from Stan\n")
  stanres <- parLapplyLB(cl, 1:chains, function(i) {
    stanres <- sampling(object = stanmodel, data = standata, pars = pars,
      chains = 1, iter = eachiter + warmup, warmup = warmup,
      thin = thin, seed = seeds[1], init = init, check_data = TRUE,
      sample_file = sample_file, diagnostic_file = diagnostic_file,
      chain_id = i)
  })

  if (chains > 1) {
    if (verbose) cat("Combining chains\n")

    stanres <- tryCatch(sflist2stanfit(stanres), error = function(e) return(results = e))
  } else {
    stanres <- stanres[[1]]
  }

  list(results = stanres, model = stanmodel, seeds = seeds)
}

