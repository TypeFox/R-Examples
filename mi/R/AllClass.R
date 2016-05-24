# Part of the mi package for multiple imputation of missing data
# Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015 Trustees of Columbia University
# Copyright (C) 1995-2012 The R Core Team
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## NOTE: If you change something here, also update the UML graph thingie

setClassUnion("MatrixTypeThing", c("matrix"))
setOldClass("family")
suppressWarnings(setClassUnion("WeAreFamily", c("family", "character"))) # arm + lme4 = warnings
setOldClass("mi_list")
setOldClass("mdf_list")

.known_imputation_methods <- c("ppd", "pmm", "mean", "median", "expectation", "mode", "mcar", NA_character_)
.known_families <- c("binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson",
                     "quasibinomial", "quasipoisson") # "quasi" is not supported at the moment (FIXME)
.known_links <- c("logit", "probit", "cauchit", "log", "cloglog", # for binomial()
                  "identity", "inverse",                          # for gaussian() plus "log", 
#                 "inverse", "identity", "log",                   # for Gamma()
                  "sqrt",                                         # for poisson() plus "log", and "identity",
                  "1/mu^2")                                       # for inverse.gaussian() plus "inverse", "identity" and "log"


# An important class in library(mi) is the missing_variable class, which is a virtual class
# for a variable that may (or may not) have missingness. The usual types of variables that
# we are interested in imputing all inherit (perhaps indirectly) from the missing_variable
# superclass, e.g. continuous, binary, etc. In principle, these class definitions should 
# provide ALL the necessary information for that variable, like the extent of its missingness
# and how the missing values will be (or have been) imputed. Thus, in principle, it should
# be possible to tweak the behavior of library(mi) simply by 1) creating a new class that
# inherits from the relevant existing class, 2) writing methods for the mi() and fit_model()
# generics and 3) perhaps a few other things that you will have to discover on your own.

## missing_variable is a virtual class for a variable that may (or may not) have missingness
setClass("missing_variable", 
         representation(
           variable_name = "character", # name of the variable but do not rely on for anything important 
           raw_data = "ANY",           ## DO NOT EVER CHANGE THE VALUES OF THIS SLOT
           data = "ANY",               ## Copy the raw_data into data and modify data as necessary
           n_total = "integer",         # total number of potential datapoints, i.e. length of raw_data
           all_obs = "logical",         # are ALL datapoints actually observed, i.e. not missing?
           n_obs = "integer",           # number of observed datapoints 
           which_obs = "integer",       # which datapoints are observed
           all_miss = "logical",        # are ALL datapoints missing, only true for latent variables
           n_miss = "integer",          # number of missing datapoints in the data slot (originally)
           which_miss = "integer",      # which datapoints are missing in the data slot
           n_extra = "integer",         # number of extra datapoints added (as missing)
           which_extra = "integer",     # which datapoints are extras
           n_unpossible = "integer",    # number of datapoints for which the variable could not be observed
           which_unpossible = "integer",# which datapoints could not be observed
           n_drawn = "integer",         # number of datapoints to impute
           which_drawn = "integer",     # which datapoints are imputed
           imputation_method = "character", # how to impute them
           family = "WeAreFamily",      # see help(family)
           known_families = "character",# families listed on help(family) plus multinomial()
           known_links = "character",   # see help(family)
           imputations = "MatrixTypeThing",      # iterations x n_drawn matrix of imputation history
           done = "logical",            # are we finished imputing?
           parameters = "MatrixTypeThing",       # history of estimated parameters in modeling this variable
           model = "ANY",               # last model fit
           fitted = "ANY",              # last fitted values
           "VIRTUAL"),
         prototype(
           variable_name = NA_character_,
           imputations = matrix(NA_real_, 0, 0),
           parameters = matrix(NA_real_, 0, 0),
           imputation_method = .known_imputation_methods,
           family = NA_character_,
           known_families = .known_families,
           known_links = .known_links
           ),
         validity = function(object) {
           out <- TRUE
           l <- length(object@raw_data)
           if(l == 0) return(out)
           
           if(sum(-object@n_total, object@n_obs, object@n_miss, object@n_extra, object@n_unpossible, na.rm = TRUE)) {
             out <- paste(object@variable_name, 
                          ": slots 'n_obs', 'n_miss', 'n_extra', and 'n_unpossible' must sum to 'n_total'")
           }
           else if(!(length(object@which_obs) %in%  c(0:1, object@n_obs))) {
             out <- paste(object@variable_name, ": 'n_obs' must equal the length of 'which_obs'")
           }
           else if(!(length(object@which_miss) %in% c(0:1, object@n_miss))) {
             out <- paste(object@variable_name, ": 'n_miss' must equal the length of 'which_miss'")
           }
           else if(!(length(object@which_extra) %in% c(0:1, object@n_extra))) {
             out <- paste(object@variable_name, ": 'n_extra' must equal the length of 'which_extra'")
           }
           else if(!(length(object@which_extra) %in% c(0:1, object@n_unpossible))) {
             out <- paste(object@variable_name, ": 'n_unpossible' must equal the length of 'which_unpossible'")
           }
           else if(sum(object@n_obs)) {
             temp <- sort(c(object@which_obs, object@which_miss, object@which_extra, object@which_unpossible))
             names(temp) <- NULL
             if(!identical(1:object@n_total, temp)) {
               out <- paste(object@variable_name, ": ''which_*' slots must be mutually exclusive and exhaustive")
             }
           }
           for(i in slotNames(object)) {
             if(i %in% c("raw_data", "data", "which_obs", "which_miss", "which_extra", 
                         "which_unpossible", "which_drawn", "imputation_method", "known_transformations",
                         "family", "known_families", "known_links", "levels", "cutpoints")) next
             if((l <- length(slot(object, i))) > 1) {
               out <- paste(object@variable_name, ": length of", i, "must be 0 or 1 but is", l)
               break
             }
           }
           return(out)
         }
         )

## this initialize() method gets called for everything that inherits from missing_variable
## but can be modified by a subsequently-called initialize() method
setMethod("initialize", "missing_variable", def = 
  function(.Object, NA.strings = c("", ".", "Na", "N/a", "N / a", "NaN",
                                   "Not Applicable", "Not applicable",
                                   "Not Available", "Not available", 
                                   "Not Ascertained", "Not ascertained",
                                   "Unavailable", "Unknown", "Missing",
                                   "Dk", "Don't Know", "Don't know", "Do Not Know", "Do not know"), ...) {
    .Object <- callNextMethod()
    if(length(.Object@raw_data) == 0) return(.Object)
    if(length(.Object@data) == 0) { # copy raw_data into data
      .Object@data <- .Object@raw_data
      names(.Object@data) <- .Object@variable_name
    }
    # bookkeeping
    infinites <- is.infinite(.Object@raw_data)
    if(any(infinites)) {
      warning(paste(.Object@variable_name, ": some observations are infinite, changing to NA"))
      .Object@data[infinites] <- NA
    }
    nans <- is.nan(.Object@raw_data)
    if(any(nans)) {
      warning(paste(.Object@variable_name, ": some observations are NaN, changing to NA"))
      .Object@data[nans] <- NA
    }
    NA.strings <- unique(c(NA.strings, toupper(NA.strings), tolower(NA.strings)))
    if(!is.numeric(.Object@raw_data)) for(i in seq_along(NA.strings)) {
      mark <- .Object@raw_data == NA.strings[i]
      if(any(mark, na.rm = TRUE)) {
        warning(paste(.Object@variable_name, ": some observations", NA.strings[i], "changing to NA"))
        .Object@data[mark] <- NA
      }
    }
    NAs <- which(is.na(.Object@data))
    if(length(NAs)) .Object@imputation_method <- "ppd" else .Object@imputation_method <- NA_character_
    .Object@n_miss <- length(NAs)
    .Object@which_miss <- NAs
    notNAs <- which(!is.na(.Object@data))
    .Object@n_obs <- length(notNAs)
    .Object@which_obs <- notNAs
    .Object@n_total <- length(NAs) + length(notNAs)
    .Object@all_miss <- length(notNAs) == 0
    .Object@all_obs  <- length(NAs) == 0
    if(!length(.Object@n_extra)) .Object@n_extra <- 0L
    if(!length(.Object@n_unpossible)) .Object@n_unpossible <- 0L
    .Object@n_drawn <- .Object@n_miss + .Object@n_extra
    .Object@which_drawn <- c(.Object@which_miss, .Object@which_extra)
    .Object@done <- FALSE
    return(.Object)
  })

setClass("irrelevant", representation("missing_variable"), 
         #   prototype(
         #     imputation_method = NA_character_,
         #     family = NA_character_)
         )

## a constant variable that has no missingness (and very few methods)
setClass("fixed", representation("irrelevant"), 
         validity = function(object) {
           out <- TRUE
           vals <- unique(object@raw_data)
           vals <- vals[!is.na(vals)]
           if(sum(object@n_miss)) {
             out <- paste(object@variable_name, ": fixed variables cannot have missingness")
           }
           else if(length(vals) > 1) {
             out <- paste(object@variable_name, ": purportedly 'fixed' variables cannot have multiple unique values")
           }
           return(out)
         }
         )

setClass("group", representation("irrelevant"))

## virtual class for categorical variables, which may be unordered, ordered, binary, or interval
setClass("categorical", 
         representation(
           "missing_variable", 
           levels = "character",
           "VIRTUAL"),
         prototype(
           known_families = c("multinomial", "binomial", "gaussian")
           )
         )

setMethod("initialize", "categorical", def = 
  function(.Object, ...) {
    .Object <- callNextMethod()
    l <- length(.Object@raw_data)
    if(l == 0) return(.Object)
    ## FIXME: check on the unused levels thing
    #   .Object@raw_data <- factor(.Object@raw_data)
    lev <- levels(factor(.Object@raw_data))
#     dummies <- t(sapply(.Object@raw_data, FUN = function(x) as.integer(x == lev)))[,-1, drop = FALSE]
#     if(ncol(dummies) == 1) colnames(dummies) <- .Object@variable_name
#     else colnames(dummies) <- lev[-1]
#     mark <- !apply(dummies, 2, FUN = function(x) all(x == 0, na.rm = TRUE))
#     dummies <- dummies[,mark, drop = FALSE]
#     lev <- c(lev[1], lev[-1][mark])
    .Object@levels <- lev
    .Object@data <- as.integer(factor(.Object@raw_data))
    return(.Object)
  })

## this is a hacked version of binomial()
multinomial <- function (link = "logit") 
{
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
    if (linktemp == "link") {
      warning("use of multinomial(link=link) is deprecated\n", 
              domain = NA)
      linktemp <- eval(link)
      if (!is.character(linktemp) || length(linktemp) != 
        1L) 
        stop("'link' is invalid", domain = NA)
    }
  }
  okLinks <- c("logit", "probit", "cloglog", "cauchit", "log")
  if (linktemp %in% okLinks) 
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name)) 
        linktemp <- stats$name
    }
    else {
      stop(gettextf("link \"%s\" not available for multinomial family; available links are %s", 
                    linktemp, paste(sQuote(okLinks), collapse = ", ")), 
           domain = NA)
    }
  }
  variance <- function(mu) mu * (1 - mu)
  validmu <- function(mu) all(mu > 0) && all(mu < 1)
  dev.resids <- binomial()$dev.resids
  aic <- function(y, n, mu, wt, dev) {
    m <- if (any(n > 1)) 
      n
    else wt
    -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * 
      y), round(m), mu, log = TRUE))
  }
  initialize <- expression({
    if (NCOL(y) == 1) {
      if (is.factor(y)) y <- y != levels(y)[1L]
      n <- rep.int(1, nobs)
      y[weights == 0] <- 0
      if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
      mustart <- (weights * y + 0.5)/(weights + 1)
      m <- weights * y
      if (any(abs(m - round(m)) > 0.001)) warning("non-integer #successes in a multinomial glm!")
    } else if (NCOL(y) == 2) {
      if (any(abs(y - round(y)) > 0.001)) warning("non-integer counts in a multinomial glm!")
      n <- y[, 1] + y[, 2]
      y <- ifelse(n == 0, 0, y[, 1]/n)
      weights <- weights * n
      mustart <- (n * y + 0.5)/(n + 1)
    } else stop("for the multinomial family, y must be a vector of 0 and 1's\n", 
                "or a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
  })
  simfun <- function(object, nsim) {
    ftd <- fitted(object)
    n <- length(ftd)
    ntot <- n * nsim
    wts <- object$prior.weights
    if (any(wts%%1 != 0)) 
      stop("cannot simulate from non-integer prior.weights")
    if (!is.null(m <- object$model)) {
      y <- model.response(m)
      if (is.factor(y)) {
        yy <- factor(1 + rbinom(ntot, size = 1, prob = ftd), 
                     labels = levels(y))
        split(yy, rep(seq_len(nsim), each = n))
      }
      else if (is.matrix(y) && ncol(y) == 2) {
        yy <- vector("list", nsim)
        for (i in seq_len(nsim)) {
          Y <- rbinom(n, size = wts, prob = ftd)
          YY <- cbind(Y, wts - Y)
          colnames(YY) <- colnames(y)
          yy[[i]] <- YY
        }
        yy
      }
      else rbinom(ntot, size = wts, prob = ftd)/wts
    }
    else rbinom(ntot, size = wts, prob = ftd)/wts
  }
  structure(list(family = "multinomial", link = linktemp, linkfun = stats$linkfun, 
                 linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simfun), 
            class = "family")
}

## unordered categorical, which corresponds to an unordered factor with more than 2 levels
setClass("unordered-categorical", representation("categorical", 
                                                 estimator = "character", 
                                                 use_NA = "logical",
                                                 rank = "integer"),
         prototype(
           estimator = "MNL",
           imputation_method = c("ppd", "pmm", "mode", "mcar", NA_character_),
           family = multinomial(link = "logit"),
           known_families = c("multinomial", "binomial"),
           known_links = c("logit", "probit", "cauchit", "log", "cloglog"),
           use_NA = FALSE,
           rank = NA_integer_
           ),
         validity = function(object) {
           out <- TRUE
           values <- unique(object@raw_data)
           values <- values[!is.na(values)]
           
           im <- getClass(class(object))@prototype@imputation_method
           if(length(values) > 0 && length(values) <= 2) {
             out <- paste(object@variable_name, "unordered-categoricals must have more than 2 levels; otherwise use binary")
           }
           else if(!all(object@imputation_method %in% im)) {
             out <- paste(object@variable_name, ": 'imputation_method' must be one of:\n", paste(im, collapse = ", "))
           }
#            else if(object@family$family != "multinomial") {
#              out <- "the 'family' slot of 'unordered-categorial' class must be 'multinomial(link = 'logit')'"
#            }
#            else if(object@family$link != "logit") {
#              out <- "the 'family' slot of 'unordered-categorial' class must be 'multinomial(link = 'logit')'"
#            }
           else if(!(object@estimator %in% c("MNL", "RNL"))) {
             out <- paste(object@variable_name, ": estimator not recognized")
           }
           else if(!(object@use_NA %in% c(TRUE, FALSE))) {
             out <- paste(object@variable_name, ": use_NA must be TRUE or FALSE")
           }
           return(out)
         }
)

## ordered categorical, which corresponds to an ordered factor
setClass("ordered-categorical", representation("categorical", 
                                               cutpoints = "numeric"),
         prototype(
           imputation_method = c("ppd", "pmm", "mode", "mcar", NA_character_),
           family = multinomial(link = "logit"),
           known_families = c("multinomial", "gaussian", "binomial", "quasibinomial"),
           known_links = "logit"
           ),
         validity = function(object) {
           out <- TRUE
           im <- getClass(class(object))@prototype@imputation_method
           if(!(object@family$family %in% getClass(class(object))@prototype@known_families)) { # interval and binary are validated separately
             out <- "the 'family' slot of 'ordered-categorial' class must be 'multinomial()'"
           }
           else if(object@family$family == "multinomial" && object@family$link != "logit") {
             out <- "the 'family' slot of 'ordered-categorial' class must be 'multinomial(link = 'logit')'"
           }
           else if(!all(object@imputation_method %in% im)) {
             out <- paste(object@variable_name, ": 'imputation_method' must be one of:\n", paste(im, collapse = ", "))
           }
           return(out)
         }
         )

## ordered categorical with known cutpoints that discretize a continuous variable (like income)
setClass("interval", representation("ordered-categorical"),
         prototype(
           imputation_method = c("ppd", NA_character_),
           family = gaussian(),
           known_families = "gaussian",
           known_links = c("identity", "inverse", "log")
           ),
         validity = function(object) {
           out <- TRUE
           if(!(object@imputation_method[1] == "ppd")) {
             out <- paste(object@variable_name, ": 'imputation_method' must be 'ppd'")
           }
           else if(object@family$family != "gaussian") {
             out <- "the 'family' slot of 'interval' class must be 'gaussian()'"
           }
           return(out)
         }
         )

## binary variable
# binary inherits from ordered-categorical because it often makes sense to think of
# those who are coded as 1 as having "more" of something than those who are coded as
# zero. Also, binary logit, probit, etc. are special cases of ordinal logit, probit,
# etc. with one cutpoint fixed at zero.
setClass("binary", representation("ordered-categorical"),
         prototype(
           family = binomial(link = "logit"),
           known_families = c("binomial", "quasibinomial"),
           known_links = c("logit", "probit", "cauchit", "log", "cloglog"),
           cutpoints = 0.0),
         validity = function(object) {
           out <- TRUE
           if(length(object@raw_data) == 0) return(out)
           vals <- unique(object@raw_data)
           vals <- vals[!is.na(vals)]
           kf <- getClass(class(object))@prototype@known_families
           kl <- getClass(class(object))@prototype@known_links
           if(length(vals) != 2) {
             out <- paste(object@variable_name, ": binary variables must have exactly two response categories")
           }
           else if(!identical(object@cutpoints, 0.0)) {
             out <- paste(object@variable_name, ": 'cutpoints' must be 0.0 for a binary variable")
           }
           else if(!(object@family$family %in% kf)) {
             out <- paste(object@variable_name, ": the 'family' slot of a object of class 'binary' must be one of", paste(kf, collapse = ", "))
           }
           else if(!(object@family$link %in% kl)) {
             out <- paste(object@variable_name, ": the 'link' slot of the 'family' slot of a object of class 'binary' must be one of", paste(kl, collapse = ", "))
           }
           return(out)
         }
         )

setMethod("initialize", "binary", def = 
  function(.Object, ...) {
    .Object <- callNextMethod()
    l <- length(.Object@raw_data)
    if(l == 0) return(.Object)
    .Object@data <- as.integer(.Object@data == max(.Object@data, na.rm = TRUE)) + 1L
    return(.Object)
  })

setClass("grouped-binary", representation("binary", strata = "character"),
         prototype(
           imputation_method = "pmm"
           ),
         validity = function(object) {
           out <- TRUE
           if(length(object@raw_data) == 0) return(out)
           if(!requireNamespace("survival")) {
             out <- "the 'survival' package must be installed to use 'grouped-binary' variables"
           }
           else if(length(object@strata) == 0) {
             warning(paste("you must specify the 'strata' slot for", object@variable_name, 
                           "see help('grouped-binary-class')"))
           }
           return(out)
         }
         )

setMethod("initialize", "grouped-binary", def = 
  function(.Object, ...) {
    .Object <- callNextMethod()
    l <- length(.Object@raw_data)
    if(l == 0) return(.Object)
    .Object@imputation_method <- "pmm"
    return(.Object)
  })

## count variables, which must be nonnegative integers
setClass("count", representation("missing_variable"),
         prototype(
           imputation_method = c("ppd", "pmm", "mean", "median", "expectation", "mcar", NA_character_),
           family = quasipoisson(),
           known_families = c("quasipoisson", "poisson"),
           known_links = c("log", "identity", "sqrt")
           ),
         validity = function(object) {
           out <- TRUE
           l <- length(object@raw_data)
           if(l == 0) return(out)
           im <- getClass(class(object))@prototype@imputation_method
           if(any(object@raw_data < 0, na.rm = TRUE)) {
             out <- paste(object@variable_name, ": counts must be nonnegative")
           }
           else if(any(object@raw_data != as.integer(object@raw_data), na.rm = TRUE)) {
             out <- paste(object@variable_name, ": must contain all nonnegative integers to use the 'count' class")
           }
           else if(!all(object@imputation_method %in% im)) {
             out <- paste(object@variable_name, ": 'imputation_method' must be one of:\n", paste(im, collapse = ", "))
           }
           else if(sum(object@n_unpossible)) {
             out <- paste(object@variable_name, ": unpossible observations not supported for count variables yet")
           }
           return(out)
         }
         )

.identity_transform <- function(y, ...) return(y)
.standardize_transform <- function(y, mean = stop("must supply mean"), sd = stop("must supply sd"), inverse = FALSE) {
  if(inverse) return(y * 2 * sd + mean)
  else return( (y - mean) / (2 * sd) )
}

## continuous variables, which may have inequality restrictions or transformation functions
setClass("continuous", 
         representation(
           "missing_variable",
           transformation = "function",
           inverse_transformation = "function",
           transformed = "logical", # TRUE -> in transformed state
           known_transformations = "character"
           ),
         prototype(
           imputation_method = c("ppd", "pmm", "mean", "median", "expectation", "mcar", NA_character_),
           transformed = TRUE,
           transformation = .standardize_transform,
           inverse_transformation = .standardize_transform,
           family = gaussian(), 
           known_families = c("gaussian", "Gamma", "inverse.gaussian", "binomial"), # binomial() is only for (SC_)proportions
           known_links = .known_links[.known_links != "sqrt"],
           known_transformations = c("standardize", "identity", "log", "logshift", "squeeze", "sqrt", "cuberoot", "qnorm")
           ),
         validity = function(object) {
           out <- TRUE
           im <- getClass(class(object))@prototype@imputation_method
           kf <- getClass(class(object))@prototype@known_families
           kl <- getClass(class(object))@prototype@known_links
           if(!all(object@imputation_method %in% im)) {
             out <- paste(object@variable_name, ": 'imputation_method' must be one of:\n", paste(im, collapse = ", "))
           }
           else if(sum(object@n_unpossible)) {
             out <- paste(object@variable_name, ": unpossible observations not supported for continuous variables yet")
           }
           else if(!(object@family$family %in% kf)) {
             out <- paste(object@variable_name, ": the 'family' slot of a object of class 'binary' must be one of", paste(kf, collapse = ", "))
           }
           else if(!(object@family$link %in% kl)) {
             out <- paste(object@variable_name, ": the 'link' slot of the 'family' slot of a object of class 'binary' must be one of", paste(kl, collapse = ", "))
           }
           return(out)
         }
         )

setMethod("initialize", "continuous", def = 
  function(.Object, ...) {
    .Object <- callNextMethod()
    l <- length(.Object@raw_data)
    if(l == 0) return(.Object)
    if(identical(.Object@transformation, .standardize_transform)) {
      mean <- mean(.Object@raw_data, na.rm = TRUE)
      sd <- sd(.Object@raw_data, na.rm = TRUE)
      formals(.Object@transformation)$mean <- formals(.Object@inverse_transformation)$mean <- mean
      formals(.Object@transformation)$sd <- formals(.Object@inverse_transformation)$sd <- sd
      formals(.Object@inverse_transformation)$inverse <- TRUE
    }
    else if(identical(.Object@transformation, .logshift)) {
      y <- .Object@raw_data
      if(any(y < 0, na.rm = TRUE)) a <- - min(y, na.rm = TRUE)
      else a <- 0
      a <- (a + min(y[y > 0], na.rm = TRUE)) / 2
      formals(.Object@transformation)$a <- formals(.Object@inverse_transformation)$a <- a
      formals(.Object@inverse_transformation)$inverse <- TRUE
    }
    .Object@data <- .Object@transformation(.Object@raw_data)
    .Object@data[.Object@which_miss] <- NA_real_
    return(.Object)
  })

setClass("bounded-continuous", representation("continuous", lower = "numeric", upper = "numeric"),
         prototype(
           imputation_method = "ppd",
           transformation = .identity_transform,
           inverse_transformation = .identity_transform
           ), 
         validity = function(object) {
           out <- TRUE
           #     if(any(object@raw_data <= object@lower, na.rm = TRUE)) {
           #       out <- paste(object@variable_name, ": all observed data must be strictly greater than 'lower'")
           #     }
           #     else if(any(object@raw_data >= object@upper, na.rm = TRUE)) {
           #       out <- paste(object@variable_name, ": all observed data must be strictly less than 'upper'")
           #     }
           if(any(object@lower > object@upper)) {
             out <- paste(object@variable_name, ": lower bounds must be less than or equal to upper bounds")
           }
           else if(object@imputation_method != "ppd") {
             out <- paste(object@variable_name, ": 'imputation_method' must be 'ppd' for 'bounded-continuous' variables ")
           }
           else if(!requireNamespace("truncnorm")) {
             out <- paste(object@variable_name, ": the 'truncnorm' package must be installed to use the 'bounded-continuous' class")
           }
           return(out)
         }
         )

setMethod("initialize", "bounded-continuous", def = 
  function(.Object, lower = -Inf, upper = Inf, ...) {
    .Object <- callNextMethod()
    l <- length(.Object@raw_data)
    if(l == 0) return(.Object)
    .Object@lower <- lower
    .Object@upper <- upper
    return(.Object)
  })

setClass("positive-continuous", representation("continuous"),
         prototype(
           transformation = log,
           inverse_transformation = exp,
           known_transformations = c("log", "sqrt", "squeeze", "qnorm") 
           ),
         validity = function(object) {
           out <- TRUE
           if(any(object@raw_data <= 0, na.rm = TRUE)) {
             out <- paste(object@variable_name, ": positive variables must be positive")
           }
           return(out)
         }
         )

## must be on the (0,1) interval
setClass("proportion", representation("positive-continuous", link.phi = "WeAreFamily"),
         prototype(
           transformed = FALSE,
           transformation = .identity_transform,
           inverse_transformation = .identity_transform,
           known_transformations = c("squeeze", "qnorm"),
           family = binomial(),
           known_families = c("binomial", "gaussian"),
           known_links = .known_links[.known_links != "sqrt"],
           link.phi = "log"),
         validity = function(object) {
           out <- TRUE
           kf <- getClass(class(object))@prototype@known_families
           kl <- getClass(class(object))@prototype@known_links
           if(any(object@raw_data > 1, na.rm = TRUE)) {
             out <- paste(object@variable_name, ": proportions must be on the unit interval")
           }
           else if(any(object@raw_data == 1, na.rm = TRUE)) {
             out <- paste(object@variable_name, ": some proportions are equal to 1.0 so use the SC_proportion class")
           }
           else if(!(object@family$family %in% kf)) {
             out <- paste(object@variable_name, ": the 'family' slot of a object of class 'proportion' must be one of", paste(kf, collapse = ", "))
           }
           else if(!(object@family$link %in% kl)) {
             out <- paste(object@variable_name, ": the 'link' slot of the 'family' slot of a object of class 'proportion' must be one of", paste(kl, collapse = ", "))
           }
           else if(object@family$family == "binomial" && !requireNamespace("betareg")) {
             out <- paste(object@variable_name, ": you must install the 'betareg' package to model proportions as proportions")
           }
           return(out)
         }
         )

# setClass("truncated-continuous", 
#   representation("continuous",
#                  lower = "ANY", 
#                  upper = "ANY", 
#                  n_lower = "integer",
#                  which_lower = "integer",
#                  n_upper = "integer",
#                  which_upper = "integer",
#                  n_both = "integer",
#                  which_both = "integer",
#                  n_truncated = "integer",
#                  which_truncated = "integer",
#                  "VIRTUAL")
# )
# 
# setClass("NN_truncated-continuous", representation("truncated-continuous", lower = "numeric", upper = "numeric"))
# 
# setMethod("initialize", "NN_truncated-continuous", def = 
# function(.Object, ...) {
#   .Object <- callNextMethod()
#   l <- length(.Object@raw_data)
#   if(l == 0) return(.Object)
#   if(identical(.Object@transformation, .standardize_transform)) {
#     mean <- mean(.Object@raw_data, na.rm = TRUE)
#     sd <- sd(.Object@raw_data, na.rm = TRUE)
#     formals(.Object@transformation)$mean <- formals(.Object@inverse_transformation)$mean <- mean
#     formals(.Object@transformation)$sd <- formals(.Object@inverse_transformation)$sd <- sd
#     formals(.Object@inverse_transformation)$inverse <- TRUE
#   }
#   .Object@data <- .Object@transformation(.Object@raw_data)
# 
#   if(length(.Object@lower) == 0 & length(.Object@upper) == 0) {
#     stop("at least one of 'lower' and 'upper' must be specified")
#   }
#   ## FIXME: Deal with interval censoring or force it to the interval class
#   .Object@n_both <- 0L
#   lowers <- .Object@raw_data <= .Object@lower
#   .Object@n_lower <- sum(lowers)
#   .Object@which_lower <- which(lowers)
#   uppers <- .Object@raw_data >= .Object@upper
#   .Object@n_uppers <- sum(uppers)
#   .Object@which_uppers <- which(uppers)
#   .Object@n_truncated <- .Object@n_lower + .Object@n_upper
#   .Object@which_truncated <- c(.Object@which_lower, .Object@which_upper)
#   return(.Object)
# })
# 
# setClass("FN_truncated-continuous", representation("truncated-continuous", lower = "function", upper = "numeric"))
# setClass("NF_truncated-continuous", representation("truncated-continuous", lower = "numeric", upper = "function"))
# setClass("FF_truncated-continuous", representation("truncated-continuous", lower = "function", upper = "function"))
# 
# setClass("censored-continuous", 
#   representation("continuous", 
#                  lower = "ANY", 
#                  upper = "ANY", 
#                  n_lower = "integer",
#                  which_lower = "integer",
#                  n_upper = "integer",
#                  which_upper = "integer",
#                  n_both = "integer",
#                  which_both = "integer",
#                  n_censored = "integer",
#                  which_censored = "integer",
#                  lower_indicator = "binary",
#                  upper_indicator = "binary",
#                  "VIRTUAL")
# )
# setClass("NN_censored-continuous", representation("censored-continuous", lower = "numeric", upper = "numeric"))
# setMethod("initialize", "NN_censored-continuous", def = 
# function(.Object, ...) {
#   .Object <- callNextMethod()
#   l <- length(.Object@raw_data)
#   if(l == 0) return(.Object)
#   if(identical(.Object@transformation, .standardize_transform)) {
#     mean <- mean(.Object@raw_data, na.rm = TRUE)
#     sd <- sd(.Object@raw_data, na.rm = TRUE)
#     formals(.Object@transformation)$mean <- formals(.Object@inverse_transformation)$mean <- mean
#     formals(.Object@transformation)$sd <- formals(.Object@inverse_transformation)$sd <- sd
#     formals(.Object@inverse_transformation)$inverse <- TRUE
#   }
#   .Object@data <- .Object@transformation(.Object@raw_data)
# 
#   if(length(.Object@lower) == 0 & length(.Object@upper) == 0) {
#     stop("at least one of 'lower' and 'upper' must be specified")
#   }
#   ## FIXME: Deal with interval censoring or force it to the interval class
#   .Object@n_both <- 0L
#   lowers <- .Object@raw_data <= .Object@lower
#   .Object@n_lower <- sum(lowers, na.rm = TRUE)
#   .Object@which_lower <- which(lowers)
#   if(.Object@n_lower > 0) {
#     .Object@lower_indicator <- missing_variable(as.ordered(lowers), type = "binary", 
#                                  variable_name = paste(.Object@variable_name, "lower", sep = ""))
#   }
#   uppers <- .Object@raw_data >= .Object@upper
#   .Object@n_upper <- sum(uppers, na.rm = TRUE)
#   .Object@which_upper <- which(uppers)
#   if(.Object@n_upper > 0) {
#     .Object@lower_indicator <- missing_variable(as.ordered(uppers), type = "binary", 
#                                  variable_name = paste(.Object@variable_name, "upper", sep = ""))
#   }
#   .Object@n_censored <- .Object@n_lower + .Object@n_upper
#   .Object@which_censored <- c(.Object@which_lower, .Object@which_upper)
#   return(.Object)
# })
# 
# setClass("FN_censored-continuous", representation("censored-continuous", lower = "function", upper = "numeric"))
# setClass("NF_censored-continuous", representation("censored-continuous", lower = "numeric", upper = "function"))
# setClass("FF_censored-continuous", representation("censored-continuous", lower = "function", upper = "function"))

setClass("semi-continuous", representation("continuous", indicator = "ordered-categorical"),
         prototype(
           transformation = .identity_transform,
           inverse_transformation = .identity_transform)
         )

.logshift <- function(y, a, inverse = FALSE) {
  if(inverse) exp(y) - a
  else log(y + a)
}

setClass("nonnegative-continuous", representation("semi-continuous"),
         prototype(transformation = .logshift,
                   inverse_transformation = .logshift,
                   known_transformations = c("logshift", "squeeze", "identity")),
         validity = function(object) {
           out <- TRUE
           if(any(object@raw_data < 0, na.rm = TRUE)) {
             out <- paste(object@variable_name, ": nonnegative variables must be nonnegative")
           }
           return(out)
         }
         )

setMethod("initialize", "nonnegative-continuous", def = 
  function(.Object, ...) {
    .Object <- callNextMethod()
    l <- length(.Object@raw_data)
    if(l == 0) return(.Object)
    is_zero <- as.integer(.Object@raw_data == 0)
    if(any(is_zero, na.rm = TRUE)) {
      .Object@indicator <- missing_variable(is_zero, type = "binary", 
                                            variable_name = paste(.Object@variable_name, ":is_zero", sep = ""))
    }
    .Object@data <- .Object@transformation(.Object@raw_data)
    if(!all(is.finite(.Object@data[!is.na(.Object)]))) {
      stop(paste(.Object@variable_name, ": some transformed values are infinite or undefined"))
    }
    return(.Object)
  })

.squeeze_transform <- function(y, inverse = FALSE) {
  n <- length(y)
  if(inverse) (y * n - .5) / (n - 1)
  else (y * (n - 1) + .5) / n
}

## some values are zero and / or one
setClass("SC_proportion", representation("nonnegative-continuous", link.phi = "WeAreFamily"),
         prototype(
           transformation = .squeeze_transform,
           inverse_transformation = .squeeze_transform,
           known_transformations = c("squeeze", "qnorm"),
           family = binomial(),
           known_families = "binomial",
           known_links = getClass("binary")@prototype@known_links,
           link.phi = "log"
           ),
         validity = function(object) {
           out <- TRUE
           if(any(object@data > 1, na.rm = TRUE)) {
             out <- paste(object@variable_name, ": proportions must be less than or equal to 1")
           }
           else if(object@family$family != "binomial") {
             out <- paste(object@variable_name, ": 'family' must be 'binomial'")
           }
           else if(!identical(body(object@transformation), body(.squeeze_transform))) {
             out <- paste(object@variable_name, ": 'transformation' must be 'squeeze'")
           }
           else if(!requireNamespace("betareg")) {
             out <- paste(object@variable_name, ": you must install the 'betareg' package to model proportions")
           }
           return(out)
         }
         )

setMethod("initialize", "SC_proportion", def = 
  function(.Object, ...) {
    .Object <- callNextMethod()
    l <- length(.Object@raw_data)
    if(l == 0) return(.Object)
    
    if(any(.Object@raw_data == 0, na.rm = TRUE)) {
      if(any(.Object@raw_data == 1, na.rm = TRUE)) {
        is_bound <- ifelse(.Object@raw_data == 0, -1, ifelse(.Object@raw_data == 1, 1, 0))
        .Object@indicator <- missing_variable(is_bound, type = "ordered-categorical", 
                                              variable_name = paste(.Object@variable_name, ":is_bound", sep = ""))
      }
      else {
        is_zero <- as.integer(.Object@raw_data == 0)
        .Object@indicator <- missing_variable(is_zero, type = "binary", 
                                              variable_name = paste(.Object@variable_name, ":is_zero", sep = ""))
      }
    }
    else {
      is_one <- as.integer(.Object@raw_data == 1)
      .Object@indicator <- missing_variable(is_one, type = "binary", 
                                            variable_name = paste(.Object@variable_name, ":is_one", sep = ""))
    }
    return(.Object)
  })

# A missing_data.frame is a another important S4 class that is not unlike a data.frame, except
# that its "columns" (actually list elements) are objects that inherit from the missing_variable
# class. The missing_data.frame class should, in principle, contain ALL the necessary information
# regarding how the missing_variables relate to each other. Together, the missing_variable class(es)
# and the missing_data.frame class supplant the mi.info S4 class in previous versions of library(mi).

.get_slot <- 
  function(object, name, simplify = TRUE) {
    if(isS4(object)) return(slot(object, name))
    else if(is.list(object)) sapply(object, FUN = slot, name = name, simplify = simplify)
    else stop("'object' not supported")
  }

setOldClass("data.frame")
setClass("missing_data.frame", 
         representation(       
           variables = "list",      # of missing_variables
           no_missing = "logical",  # basically a collection of the all_obs slots of the missing_variables
           patterns = "factor",     # indicates which missingness_pattern an observation belongs to
           DIM = "integer",         # observations x variables
           DIMNAMES = "list",       # list of rownames and colnames
           postprocess = "function",# makes additional variables from existing variables (interactions, etc.)
           index = "list",          # this indicate which variables to exclude when modeling a given variable
           X = "MatrixTypeThing",   # ALL variables (categorical variables are in dummy-variable form)
           weights = "list",        # this gets passed to bayesglm() and similar modeling functions
           priors = "list",         # the elements of this get passed to bayesglm() and other modeling functions in arm
           correlations = "matrix", # has SMCs and Spearman correlations
           done = "logical",        # are we done?
           workpath = "character"),
         contains = "data.frame", 
         prototype(postprocess = function() stop("postprocess does not work yet"), X = matrix(NA_real_, 0, 0), done = FALSE),
         validity = function(object) {
           out <- TRUE
           l <- length(object@variables)
           if(l == 0) return(out)
           
           if(!all(sapply(object@variables, FUN = is, class2 = "missing_variable"))) {
             out <- "all of the list elements in 'variables' must inherit from the 'missing_variable' class"
           }
           else if(length(unique(.get_slot(object@variables, "n_total"))) > 1) {
             out <- "all missing_variables must have the same 'n_total'"
           }
           else if(!is.numeric(object@X)) {
             out <- "'X' must be a numeric matrix"
           }
           missingness <- .get_slot(object@variables, "which_miss", simplify = FALSE)
           varnames <- .get_slot(object@variables, "variable_name")
           names(missingness) <- varnames
           missingness <- missingness[sapply(missingness, length) > 0]
           if(length(missingness) > 1) { ## FIXME: Very slow
             combos <- combn(length(missingness), 2) 
             dupes <- apply(combos, 2, FUN = function(x) {
               mx1 <- missingness[[x[1]]]
               mx2 <- missingness[[x[2]]]
               if(length(mx1) == length(mx2)) {
                 if(identical(mx1, mx2)) return(1L)
               }
               else if(length(mx1) > length(mx2)) {
                 if(all(mx2 %in% mx1)) return(2L)
               }
               else if(all(mx1 %in% mx2)) return(3L)
               return(0L)
             })
             if(any(dupes == 1L)) {
               temp <- matrix(names(missingness)[combos[,which(dupes == 1L)]], ncol = 2, byrow = TRUE)
               cat("NOTE: The following pairs of variables appear to have the same missingness pattern.\n",
                   "Please verify whether they are in fact logically distinct variables.\n")
               print(temp)
               # 	warning("Potentially duplicated variables detected by duplicated variable detector")
             }
             else if(any(dupes == 2L)) {
               temp <- matrix(names(missingness)[combos[,which(dupes == 2L)]], ncol = 2, byrow = TRUE)
               cat("NOTE: In the following pairs of variables, the missingness pattern of the second is a subset of the first.\n",
                   "Please verify whether they are in fact logically distinct variables.\n")
               print(temp)
             }
             else if(any(dupes == 3L)) {
               temp <- matrix(names(missingness)[combos[,which(dupes == 3L)]], ncol = 2, byrow = TRUE)
               cat("NOTE: In the following pairs of variables, the missingness pattern of the first is a subset of the second.\n",
                   "Please verify whether they are in fact logically distinct variables.\n")
               print(temp)
             }
           }
           return(out)
         }
         )

.set_priors <- function(variables, mu = 0) { ## FIXME: maybe add an option to draw from such a t distribution?
  foo <- function(y) {
    out <- list(prior.mean = mu, prior.scale  = 2.5, prior.df = 1, 
                prior.mean.for.intercept = mu, prior.scale.for.intercept = 10, prior.df.for.intercept = 1)
    if(is(y, "irrelevant") | y@all_obs) return(NULL)
    else if(is(y, "binary")) {
      if(y@family$link == "probit") {
        out[[2]] <- out[[2]] * dnorm(0) / dlogis(0)
        out[[4]] <- out[[4]] * dnorm(0) / dlogis(0)
      }
    }
    else if(is(y, "categorical")) {
      out <- list(prior.mean = mu, prior.scale = 2.5, prior.df = 1, prior.counts.for.bins = 1/(1 + length(y@levels)))
    }
    return(out)
  }
  
  out <- lapply(variables, FUN = function(y) foo(y))
  for(i in seq_along(variables)) if(is(y <- variables[[i]], "semi-continuous")) out[[y@indicator@variable_name]] <- foo(y@indicator)
  return(out)
}

setMethod("initialize", "missing_data.frame", def =
  function(.Object, include_missingness = TRUE, skip_correlation_check = FALSE, ...) {
    .Object <- callNextMethod()
    l <- length(.Object@variables)
    if(l == 0) return(.Object)
    varnames <- names(.Object@variables)
    if(is.null(varnames)) {
      if(is.null(.Object@DIMNAMES[[2]])) names(.Object@variables) <- sapply(.Object@variables, FUN = .get_slot, name = "variable_name")
      else names(.Object@variables) <- .Object@DIMNAMES[[2]]
    }
    else for(i in 1:l) .Object@variables[[i]]@variable_name <- varnames[i]
    .Object@DIM <- c(.Object@variables[[1]]@n_total, l)
    .Object@no_missing <- sapply(.Object@variables, FUN = .get_slot, name = "all_obs")
    if(length(.Object@DIMNAMES) == 0) .Object@DIMNAMES <- list(NULL, names(.Object@variables))

    Z <- lapply(.Object@variables, FUN = function(y) {
      if(is(y, "irrelevant")) return(NULL) else return(is.na(y))
    })
    Z <- as.matrix(as.data.frame(Z[!sapply(Z, is.null)]))
    if(any(apply(Z, 1, all))) {
      warning("Some observations are missing on all included variables.\n",
              "Often, this indicates a more complicated model is needed for this missingness mechanism")
    }
    uZ <- unique(Z)
    if(nrow(uZ) == 1) {
      if(all(uZ[1,] == 0)) patterns <- factor(rep("nothing", nrow(Z)))
      else patterns <- factor(colnames(uZ)[which(uZ[1,] == 1)], nrow(Z))
    }
    else {
      uZ <- uZ[order(rowSums(uZ)),,drop = FALSE]
      patterns <- apply(Z, 1, FUN = function(x) which(apply(uZ, 1, FUN = function(u) all(u == x))))
      pattern_labels <- apply(uZ, 1, FUN = function(x) paste(names(x)[x], collapse = ", "))
      if(length(pattern_labels)) {
        if(pattern_labels[1] == "") pattern_labels[1] <- "nothing"
        pattern_lables <- paste("missing:", pattern_labels)
        patterns <- factor(patterns, labels = pattern_labels, ordered = FALSE)
      }
      else patterns <- factor(patterns)
    }
    .Object@patterns <- patterns

    if(!length(.Object@workpath)) {
      .Object@workpath <- file.path(tempdir(), paste("mi", as.integer(Sys.time()), sep = ""))
    }
    dir.create(.Object@workpath, showWarnings = FALSE)

    if(is(.Object, "allcategorical_missing_data.frame")) return(.Object)
    
    Z <- Z[,!duplicated(t(Z)), drop = FALSE]
    Z <- Z[,apply(Z, 2, FUN = function(x) length(unique(x))) > 1, drop = FALSE]
    ## FIXME: What to do if two columns of Z are collinear?
    if(ncol(Z) > 0) colnames(Z) <- paste("missing", colnames(Z), sep = "_")
    else include_missingness <- FALSE

    X <- lapply(.Object@variables, FUN = function(x) {
      if(is(x, "irrelevant"))  return(NULL)
      else if(is(x, "categorical")) return(.cat2dummies(x))
      else if(is(x, "semi-continuous")) {  
        out <- cbind(x@data, .cat2dummies(x@indicator))
        colnames(out) <- c(x@variable_name, 
                           paste(x@variable_name, 2:ncol(out) - 1, sep = "_"))
        return(out)
      }
      else if(is(x, "censored-continuous")) {
        temp <- x@data
        if(x@n_lower) temp <- cbind(temp, lower = x@lower_indicator@data)
        if(x@n_upper) temp <- cbind(temp, upper = x@upper_indicator@data)
        if(x@n_both)  stop("FIXME: censoring on both sides not supported yet")
        return(temp)
      }
      else if(is(x, "truncated-continuous")) {
        temp <- x@data
        n <- length(temp)
        if(x@n_lower) temp <- cbind(lower = x@lower_indicator@data, temp)
        if(x@n_upper) temp <- cbind(upper = x@upper_indicator@data, temp)
        if(x@n_both)  stop("FIXME: censoring on both sides not supported yet")
        return(temp)
      }
      else return(x@data)
    }) ## NOTE: Might need to make this more complicated in the future
    X <- X[!sapply(X, is.null)]
    
    index <- vector("list", length = length(X))
    names(index) <- names(X)
    start <- 2L
    end <- 0L
    for(i in seq_along(index)) {
      end <- start + NCOL(X[[i]]) - 1L
      index[[i]] <- start:end
      start <- end + 1L
    }
    if(include_missingness) for(i in seq_along(index)) {
      nas <- is.na(.Object@variables[[i]])
      check <- apply(Z, 2, FUN = function(x) all(x == nas))
      index[[i]] <- c(index[[i]], which(check) + start - 1)
    }
    else for(i in seq_along(index)) index[[i]] <- c(index[[i]], start:(start + ncol(Z) - 1))
    
    grouped <- names(which(sapply(.Object@variables, is, class2 = "grouped-binary")))
    for(i in grouped) index[[i]] <- c(index[[i]], index[[.Object@variables[[i]]@strata]], 1)
    
    .Object@index <- index
    .Object@X <- cbind("(Intercept)" = 1, as.matrix(as.data.frame(X)), Z)
    
    correlations <- matrix(NA_real_, l,l)
    if(!skip_correlation_check) for(i in 1:(l - 1)) { ## FIXME: Put SMCs in the lower triangle
      if(is(.Object@variables[[i]], "irrelevant")) next
      x <- try(rank(xtfrm(.Object@variables[[i]]@raw_data)), silent = TRUE)
      if(!is.numeric(x)) next
      for(j in (i + 1):l) {
        if(is(.Object@variables[[j]], "irrelevant")) next
        y <- try(rank(xtfrm(.Object@variables[[j]]@raw_data)))
        if(!is.numeric(y)) next
        rho <- cor(x, y, use = "pair", method = "pearson") # on ranks
        if(is.finite(rho) && abs(rho) == 1) {
          warning(paste(names(.Object@variables)[i], "and", names(.Object@variables)[j], 
                        "have the same rank ordering.\n",
                        "Please verify whether they are in fact distinct variables.\n"))
        }
        if(is.finite(rho)) correlations[i,j] <- rho
      }
    }
    
    .Object@correlations <- correlations
    .Object@priors <- .set_priors(.Object@variables)
    .Object
  })

setClass("allcategorical_missing_data.frame", 
         representation("missing_data.frame", "Hstar" = "integer",
                        "parameters" = "list","latents" = "unordered-categorical"),
         
         prototype = prototype(Hstar = 20L),
         validity = function(object) {
           out <- TRUE
           types <- sapply(object@variables, 
                           FUN = function(y) is(y, "irrelevant") | is(y, "categorical"))
           if(!all(types)) {
             out <- "all variable classes must be 'irrelevant' or 'categorical'"
           }
           else if(length(object@Hstar) && object@Hstar < 1) {
             out <- "'Hstar' must be >= 1"
           }
           return(out)
         })

setMethod("initialize", "allcategorical_missing_data.frame", def =
  function(.Object, include_missingness = TRUE, ...) {
    .Object <- callNextMethod()
    l <- length(.Object@variables)
    n <- nrow(.Object)
    uc <- factor(rep(NA_integer_, n))
    .Object@latents <- new("unordered-categorical", raw_data = rep(NA_integer_, n))
    .Object@priors <- list(a = rep(1, ncol(.Object)), a_alpha = 1, b_alpha = 1)
    names(.Object@priors$a) <- colnames(.Object)
    return(.Object)
  })

setClass("experiment_missing_data.frame", representation("missing_data.frame", 
                                                         concept = "factor",
                                                         case = "character"),
         validity = function(object) {
           out <- TRUE
           l <- length(object@concept)
           if(l != length(object@variables)) {
             out <- "length of 'concept' must equal the number of variables"
           }
           else if(!all(levels(object@concept) %in% c("outcome", "covariate", "treatment"))) {
             out <- "all elements of 'concept' must be exactly one of 'outcome', 'covariate', or 'treatment'"
           }
           else if(sum(object@concept == "treatment") != 1) {
             out <- "there must be exactly one variable designated 'treatment'"
           }
           else if(!is(object@variables[[which(object@concept == "treatment")]], "binary")) {
             out <- "the 'treatment' variable must be of class 'binary'"
           }
           else if(object@variables[[which(object@concept == "treatment")]]@n_miss) {
             out <- "'treatment' variable cannot have any missingness"
           }
           else if(length(object@case) > 1) {
             out <- "'case' must be exactly one of 'outcomes', 'covariates', or 'both'"
           }
           else if(length(object@case) && !(object@case %in% c("outcomes", "covariates", "both"))) {
             out <- "'case' must be exactly one of 'outcomes', 'covariates', or 'both'"
           }
           return(out)
         })

setMethod("initialize", "experiment_missing_data.frame", def =
  function(.Object, include_missingness = TRUE, ...) {
    .Object <- callNextMethod()
    l <- 1  ## FIXME
    if(l == 0) return(.Object)
    names(.Object@concept) <- .Object@DIMNAMES[[2]]
    outcomes   <- any(!.Object@no_missing[.Object@concept == "outcomes"])
    covariates <- any(!.Object@no_missing[.Object@concept == "covariates"])
    
    .Object@case <- if(outcomes & covariates) "both" else if(outcomes) "covariates" else "outcomes"
    return(.Object)
  })

.empty_mdf_list <- list()
class(.empty_mdf_list) <- "mdf_list"
setClass("multilevel_missing_data.frame", representation("missing_data.frame", 
                                                         groups = "character", 
                                                         mdf_list = "mdf_list"),
         prototype(
           mdf_list = .empty_mdf_list
           ),
         
         validity = function(object) {
           out <- TRUE
           return(out)
         }
         )

setMethod("initialize", "multilevel_missing_data.frame", def =
  function(.Object, include_missingness = TRUE, ...) {
    .Object <- callNextMethod()
    classes <- sapply(.Object@variables, class)
    for(i in .Object@groups) classes[names(classes) == i] <- "fixed"
    df <- complete(.Object, m = 0L)
    mdf_list <- missing_data.frame(df, by = .Object@groups, types = classes)
    .Object@mdf_list <- mdf_list
    return(.Object)
  })

## an object of class mi merely holds the results of a call to mi(), primary the list of missing_data.frames
setClass("mi",
         representation(
           call       = "call",
           data = "list",             # of missing_data.frames
           total_iters  = "integer"), # how many iterations were conducted (can be a vector)
         )

## an object of class pooled has regression results using the Rubin rules
setClass("pooled", 
         representation(
           formula = "formula",
           fit = "character",
           models = "list",
           coefficients = "numeric",
           ses = "numeric",
           pooled_summary  = "ANY",
           call = "language"),
         )
