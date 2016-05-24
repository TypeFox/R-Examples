###########################################################################
##                                                                       ##
## mcarlo() - Monte Carlo simulation & p-values for MAT models           ##
##                                                                       ##
## Created       : 07-Jul-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 02-Nov-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
##                                                                       ##
###########################################################################
mcarlo <- function(object, ...) UseMethod("mcarlo")

mcarlo.default <- function(object, nsamp = 10000,
                           type = c("paired", "complete", "bootstrap", "permuted"),
                           replace = FALSE,
                           method = c("euclidean", "SQeuclidean",
                             "chord", "SQchord", "bray", "chi.square",
                             "SQchi.square", "information",
                             "chi.distance", "manhattan", "kendall", "gower",
                             "alt.gower", "mixed" ),
                           is.dcmat = FALSE, diag = FALSE,
                           ...) {
  if(missing(type))
    type <- "paired"
  TYPE <- match.arg(type)
  if(missing(method))
    method <- "euclidean"
  METHOD <- match.arg(method)
  if(is.dcmat)
    Dij <- object
  else
    Dij <- distance(object, method = METHOD)
  ## number of samples in modern set
  N <- nrow(object)
  mc.samp <- switch(TYPE,
                    paired = {
                      samp <- numeric(length = nsamp)
                      Dij <- Dij[lower.tri(Dij, diag = diag)]
                      for(i in 1:nsamp) {
                        choices <- sample(1:N, 1, replace = replace)
                        samp[i] <- Dij[choices]
                      }
                      samp
                    },
                    complete = stop("type = \"complete\" not yet implemented."),
                    bootstrap = {
                      if(diag)
                        num.dists <- ceiling((N*N) / 2)
                      else
                        num.dists <- (N*(N-1)) / 2
                      samp <- matrix(ncol = num.dists, nrow = nsamp)
                      for(i in 1:nsamp) {
                        choices <- sample(1:num.dists, num.dists, replace = TRUE)
                        samp[i, ] <- Dij[lower.tri(Dij, diag = diag)][choices]
                      }
                      as.vector(samp)
                    },
                    permuted = stop("type = \"permuted\" not yet implemented."))
  class(samp) <- c("mcarlo", type)
  attr(samp, "method") <- METHOD
  return(samp)
}

mcarlo.mat <- function(object, nsamp = 10000,
                       type = c("paired", "complete", "bootstrap", "permuted"),
                       replace = FALSE, diag = FALSE,
                       ...) {
  if(!inherits(object, "mat"))
    stop("'mat' method for 'mcarlo' only for objects of class 'mat'.")
  if(missing(type))
    type <- "paired"
  TYPE <- match.arg(type)
  METHOD <- object$method
  retval <- mcarlo.default(object$Dij, nsamp = nsamp,
                           type = TYPE, replace = replace, method = METHOD,
                           is.dcmat = TRUE, diag = diag,  ...)
  retval
}

mcarlo.analog <- function(object, nsamp = 10000,
                          type = c("paired", "complete", "bootstrap", "permuted"),
                          replace = FALSE, diag = FALSE,
                          ...) {
  if(!inherits(object, "analog"))
    stop("'object' not of class \"analog\".")
  if(missing(type))
    type <- "paired"
  TYPE <- match.arg(type)
  METHOD <- object$method
  if(is.null(object$train))
    stop("'object$train' missing. Refit 'object' with argument 'keep.train = TRUE'")
  retval <- mcarlo.default(object$train, nsamp = nsamp,
                           type = TYPE, replace = replace, method = METHOD,
                           is.dcmat = TRUE, diag = diag,  ...)
  retval
}

print.mcarlo <- function(x,
                         probs = c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99),
                         digits = min(3, getOption("digits") - 4),
                         ...) {
  summ <- fivenum(x)
  summ <- c(summ[1:3], mean(x), summ[4:5])
  names(summ) <- c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max")
  quant <- quantile(x, probs = probs)
  cat("\n")
  writeLines(strwrap("Simulated Dissimilarities", prefix = "\t"))
  cat("\n")
  cat(paste("Simulation type :", class(x)[2], "\n"))
  cat(paste("No. simulations :", length(x), "\n"))
  cat(paste("Coefficient     :", attr(x, "method"), "\n"))
  cat("\nSummary of simulated distribution:\n")
  print(summ, digits = digits)
  cat("\nPercentiles of simulated distribution:\n")
  print(quant, digits = digits)
  cat("\n")
  invisible(x)
}
