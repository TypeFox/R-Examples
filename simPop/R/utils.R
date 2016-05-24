# calculate parameters for parallel processing
# returns a list with following elements
# - have_win: TRUE if windows system, FALSE else
# - nr_cores: number of cpus to use
# - parallel: TRUE if parallel computing should be applied
parallelParameters <- function(nr_cpus=NULL, nr_strata) {
  control <- list()
  control$have_win <- Sys.info()["sysname"] == "Windows"

  cpus_available <- detectCores()
  if ( !is.null(nr_cpus) && nr_cpus > cpus_available) {
    stop("there are only",cpus_available,"cpus available in your system!\n")
  }
  if ( !is.null(nr_cpus) && nr_cpus < 1) {
    stop("we must use at least one cpu!\n")
  }

  if ( !is.null(nr_cpus) ) {
    control$nr_cores <- min(nr_cpus, nr_strata)
  } else {
    if ( cpus_available > 1 & nr_strata > 1 ) {
      control$nr_cores <- min(cpus_available-1, nr_strata) # keep at least one core available
    } else {
      control$nr_cores <- 1
    }
  }
  control$parallel <- ifelse(control$nr_cores>1, TRUE, FALSE)
  control
}

## check break points for categorization
checkBreaks <- function(x) {
  if(!is.numeric(x) || length(x) < 2) {
    stop("'breaks' must be a numeric vector with length >= 2")
  }
  invisible()
}

## check whether selected columns of a data frame are factors
# TODO: add argument 'num': a logical indicating whether numeric variables other than integers should be converted to factors
checkFactor <- function(x, select) {
  if ( "data.table" %in% class(x) ) {
    convert <- select[sapply(x[, select, with=F], function(x) !is.factor(x))]
    if ( length(convert) > 0 ) {
      x[,(convert):=lapply(.SD, as.factor),.SDcols=convert]
    }
  } else {
    convert <- select[sapply(x[, select], function(x) !is.factor(x))]
    n <- length(convert)
    if ( n == 1 ) {
      x[[convert]] <- as.factor(x[[convert]])
    } else {
      x[, convert] <- as.data.frame(lapply(x[, convert], as.factor))
    }
  }
  invisible(x)
}

# drop unused factor levels for factors in a data.frame
dropLevels <- function(x, select = names(x)) {
  check <- select[sapply(x[, select], function(x) is.factor(x))]
  n <- length(check)
  if(n == 1) x[, check] <- x[, check][, drop=TRUE]
  else if(n > 1) {
    x[, check] <- as.data.frame(lapply(x[, check], "[", drop=TRUE))
  }
  x
}

# # create factor that includes NA
# factorNA <- function(x, always = FALSE) {
#   always <- isTRUE(always)
#   if(is.factor(x)) {
#     l <- levels(x)
#     if(NA %in% l || !(always || any(is.na(x)))) x
#     else {
#       l <- c(l, NA)
#       factor(x, levels=c(levels(x), NA), exclude=c())
#     }
#   } else {
#     if(always) {
#       factor(c(NA, x), exclude=c())[-1]  # little trick
#     } else factor(x, exclude=c())
#   }
# }


## get which observations contain NAs (and need to be excluded)
getExclude <- function(x, ...) UseMethod("getExclude")
getExclude.default <- function(x, ...) which(is.na(x))
getExclude.data.frame <- function(x, ...) {
  unique(which(is.na(x), arr.ind=TRUE)[, 1])
}
getExclude.data.table <- function(x, ...) {
  unique(which(is.na(x), arr.ind=TRUE)[, 1])
}


### exclude observations
#excludeData <- function(x, exclude = NULL, ...) UseMethod("excludeData")
#excludeData.default <- function(x, exclude = NULL, ...) {
#	if(length(exclude) == 0) x else x[-exclude]
#}
#excludeData.data.frame <- function(x, exclude = NULL, ...) {
#	if(length(exclude) == 0) x else x[-exclude, , drop=FALSE]
#}

## adjust probabilities estimated with multinomial model to account for
## structural zeros
adjustProbs <- function(probs, grid, pNames, limit = NULL, censor = NULL) {
  set0 <- matrix(FALSE, nrow(probs), ncol(probs),
                 dimnames=dimnames(probs))
  target <- colnames(probs)
  # account for structural zeros via argument 'limit'
  if(is.list(limit)) {
    limit <- limit[sapply(limit, inherits, "list")]
    limit <- limit[names(limit) %in% names(grid)]
    # loop over predictors for which to censor probabilities
    for(i in seq_along(limit)) {
      predNameI <- names(limit)[i]
      predI <- grid[, predNameI]
      limitI <- limit[[i]]
      # loop over supplied outcomes of current predictor to
      # find probabilities to be set to zero
      for(j in seq_along(limitI)) {
        cat <- names(limitI)[j]
        ok <- intersect(target, limitI[[j]])
        if(length(ok) > 0) {
          set0[predI == cat, setdiff(target, ok)] <- TRUE
        }
      }
    }
  }
  # account for structural zeros via argument 'censor'
  if(is.list(censor)) {
    censor <- censor[sapply(censor, inherits, c("list", "data.frame"))]
    censor <- censor[names(censor) %in% target]
    if(length(censor)) {
      for(i in seq_along(censor)) {
        censorNameI <- names(censor)[i]
        censorI <- censor[[i]]
        if(is.list(censorI)) {
          # loop over supplied predictors to find
          # probabilities to be set to zero
          for(j in seq_along(censorI)) {
            predNameJ <- names(censorI)[j]
            predJ <- grid[, predNameJ]
            set0[predJ %in% censorI[[j]], censorNameI] <- TRUE
          }
        } else {
          # set probabilities to zero for supplied
          # combinations of predictors
          cNames <- apply(censorI, 1, paste, collapse=".")
          set0[pNames %in% cNames, censorNameI] <- TRUE
        }
      }
    }
  }
  # set indicated probabilities to zero
  probs[set0] <- 0
  # set the non-censored probabilities to non-zero if all
  # probabilities of a row are zero
  adjust <- which(apply(probs == 0, 1, all))
  for(i in adjust) {
    ok <- which(!set0[i,])
    probs[i, ok] <- 1/length(ok)
  }
  probs
}

## get breakpoints for categorizing continuous or semi-continuous variables







## categorize continuous or semi-continuous variables







## get name of categorized variable
getCatName <- function(name) paste(name, "Cat", sep="")


## get name of categorical variable for household head
getHeadName <- function(name) paste(name, "Head", sep="")


## truncated Pareto distribution
truncPareto <- function(n, loc, scale, shape, lower, upper) {
  # initializations
  x <- numeric(n)
  left <- n
  # fit generalized pareto distribution with lower and upper bound
  while(left > 0) {
    xTmp <- rgpd(left, loc=loc, scale=scale, shape=shape)
    ind  <- which(xTmp > lower & xTmp <= upper)
    lind <- length(ind)
    if(lind > 0) {
      x[(n - left + 1):(n - left + lind)]  <- xTmp[ind]
      left <- left - lind
    }
  }
  return(x)
}


## logit regression (designed for internal use, hence no error handling)
#logitreg <- function(x, y, weights = rep(1, length(y)),
#    intercept = TRUE, start = rep(0, p), ...) {
#    # function to be minimized (log-likelihood)
#    fmin <- function(beta, X, y, w) {
#        p <- plogis(X %*% beta)
#        -sum(2 * w * ifelse(y, log(p), log(1-p)))
#    }
#    # gradient
#    gmin <- function(beta, X, y, w) {
#        eta <- as.numeric(X %*% beta)
#        p <- plogis(eta)
#        -2 * (w * dlogis(eta) * ifelse(y, 1/p, -1/(1-p))) %*% X
#    }
#    # some preparations
#    if(is.null(dim(x))) dim(x) <- c(length(x), 1)
#    dn <- dimnames(x)[[2]]
#    if(!length(dn)) dn <- paste("Var", 1:ncol(x), sep="")
#    p <- ncol(x) + intercept
#    if(intercept) {
#        x <- cbind(1,x)
#        dn <- c("(Intercept)", dn)
#    }
#    if(is.factor(y)) y <- (unclass(y) != 1)
#    # optimize and return result
#    fit <- optim(start, fmin, gmin, X = x, y = y, w = weights, method = "BFGS", ...)
#    names(fit$par) <- dn
#    return(fit)
#}
# intercept is assumed to be included in 'x'
logitreg <- function(x, y, weights = rep(1, length(y)), start = rep(0, p), ...) {
  # function to be minimized (log-likelihood)
  fmin <- function(beta, X, y, w) {
    p <- plogis(X %*% beta)
    -sum(2 * w * ifelse(y, log(p), log(1-p)))
  }
  # gradient
  gmin <- function(beta, X, y, w) {
    eta <- as.numeric(X %*% beta)
    p <- plogis(eta)
    -2 * (w * dlogis(eta) * ifelse(y, 1/p, -1/(1-p))) %*% X
  }
  # some preparations
  if(is.null(dim(x))){
    dim(x) <- c(length(x), 1)
  }
  dn <- dimnames(x)[[2]]
  if(!length(dn)) dn <- paste("Var", 1:ncol(x), sep="")
  p <- ncol(x)
  if(is.factor(y)) y <- (unclass(y) != 1)
  # optimize and return result
  fit <- optim(start, fmin, gmin, X = x, y = y, w = weights, method = "BFGS", ...)
  names(fit$par) <- dn
  return(fit)
}


## search for donor cells for splitting continuous variables
# i ........... indices corresponding to a zero of a table
# indDonors ... matrix containing indices of all non-zero elements
# donors ...... vector indices of the non-zero elements

# sequential minimum distance
seqMinDist <- function(i, indDonors, donors) {
  # initializations
  donor <- NA
  k <- 1
  # loop over dimensions to find donor cells with matching category
  while(is.na(donor) && k <= length(i)) {
    sel <- which(indDonors[, k] == i[k])  # donor cells in same category
    if(length(sel)) {
      tmpIndDonors <- indDonors[sel, , drop=FALSE]
      tmpDonors <- donors[sel]
      donor <- minDist(i[-k], tmpIndDonors[, -k, drop=FALSE], tmpDonors)
    }
    k <- k + 1
  }
  # if not successful, use donor with minimum overall distance
  if(is.na(donor)) donor <- minDist(i, indDonors, donors)
  donor
}

# minimum distance
minDist <- function(i, indDonors, donors) {
  # calculate distance from element defined by i for each element of m
  d <- colSums(abs(t(indDonors) - i))
  # return index of selected donor (minimum distance)
  donors[which.min(d)]
}








# prepare input for regression-based estimation methods
# returns a list with required independent variables and the model formula
# regModel can be either 'basic', 'availabe' or a formula
regressionInput <- function(obj, additional, regModel) {
  cn_pop <- colnames(popData(obj))
  cn_samp <- colnames(sampleData(obj))

  if ( length(additional) == 1 & class(regModel)=="formula" ) {
    regModel <- list(regModel)
  }

  if ( length(additional) != length(regModel) ) {
    stop("makeRegInput:: dimensions do not match!\n")
  }

  if ( length(additional) == 1 & class(regModel)=="formula" ) {
    regModel <- list(regModel)
  }

  regInput <- list(); length(regInput) <- length(additional)
  names(regInput) <- additional
  for ( i in seq_along(additional) ) {
    if ( class(regModel[[i]]) == "formula" ) {
      # check the formula
      fmt <- format(regModel[[i]])
      if ( substr(fmt,1,1)=="~" ) {
        regModel[[i]] <- as.formula(paste0(additional[i], fmt))
      }
      regvars <- all.vars(regModel[[i]])
      if ( regvars[1] != additional[i] ) {
        stop("dependent variable in model-formula does not match with variable name specified in 'additional'!\n")
      }
      regvars <- regvars[-1]
      ii <- regvars %in% cn_pop
      if ( !all(ii) ) {
        stop(paste0("Some variables (", paste0(regvars[!ii], collapse=", "), ") required for the regression model are not available in the synthetic population!\n"))
      }
      ii <- regvars %in% cn_samp
      if ( !all(ii) ) {
        stop(paste0("Some variables (", paste0(regvars[!ii], collapse=", "), ") required for the regression model are not available in the sample!\n"))
      }
      regInput[[i]]$predNames <- regvars
      regInput[[i]]$formula <- format(regModel[[i]])
    } else {
      if ( !regModel[i] %in% c("basic","available") ) {
        stop("regModel[",i,"] is neither a formula nor 'basic' or 'available'!\n")
      }
      if ( regModel[[i]] == "basic" ) {
        regInput[[i]] <- list()
        # get basic household variables and hh sizes
        regInput[[i]]$predNames <- c(obj@basicHHvars, sampleObj(obj)@hhsize)
        regInput[[i]]$formula <- paste0(additional[i],"~", paste(regInput[[i]]$predNames, collapse="+"))
      }
      if ( regModel[[i]] == "available" ) {
        # get all available variables from sample/pop (excluding ids,...)
        regInput[[i]]$predNames <- setdiff(cn_pop, c(popObj(obj)@hhid, popObj(obj)@pid, popObj(obj)@weight))
        # we also use variables previously generated
        if ( i > 1 ) {
          regInput[[i]]$predNames <- c(regInput[[i]]$predNames, additional[1:(i-1)])
        }
        regInput[[i]]$formula <- paste0(additional[i],"~", paste(regInput[[i]]$predNames, collapse="+"))
      }
    }
  }
  return(regInput)
}


# recodes a factor:
# NA -> _tmpmiss
# non-available codes are removed (droplevels)
# required when estimating a model on a strata and
# not all levels are present in the response
cleanFactor <- function(x) {
  missstr <- "_tmpmiss"
  if ( !is.factor(x) ) {
    return(x)
  }
  y <- as.character(x)
  indna <- is.na(y)
  levN <- levels(droplevels(x))
  if ( any(indna) ) {
    y[indna] <- missstr
    levN <- c(levN,missstr)
  }
  factor(y, levels=levN)
}
