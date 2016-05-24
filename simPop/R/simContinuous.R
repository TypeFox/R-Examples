generateValues_multinom <- function(dataSample, dataPop, params) {
  excludeLevels <- params$excludeLevels
  maxit <- params$maxit
  MaxNWts <- params$MaxNWts
  command <- params$command
  eps <- params$eps
  limit <- params$limit
  censor <- params$censor
  hasNewLevels  <- params$hasNewLevels
  newLevels <- params$newLevels
  name <- params$name
  response <- params$response

  # unique combinations in the stratum of the population need
  # to be computed for prediction
  indGrid <- split(1:nrow(dataPop), dataPop, drop=TRUE)
  grid <- dataPop[sapply(indGrid, function(i) i[1]), , drop=FALSE]
  grid <- as.data.frame(grid)

  # in sample, observations with NAs have been removed to fit the
  # model, hence population can have additional levels
  # these need to be removed since those probabilities cannot
  # be predicted from the model
  if ( excludeLevels ) {
    exclude <- mapply(function(pop, new) pop %in% new,
      pop=grid[, hasNewLevels, drop=FALSE], new=newLevels[hasNewLevels]
    )
    if ( is.null(dim(exclude)) ) {
      exclude <- which(any(exclude))
    } else {
      exclude <- which(apply(exclude, 1, any))
    }
  } else {
    exclude <- integer()
  }
  # fit multinomial model
  mod <- eval(parse(text=command))  # fitted model
  # predict probabilities
  if ( length(exclude) == 0 ) {
    probs <- predict(mod, newdata=grid, type="probs")
  } else {
    probs <- predict(mod, newdata=grid[-exclude, , drop=FALSE], type="probs")
  }
  # set too small probabilities to exactly 0
  if ( !is.null(eps) ) {
    probs[probs < eps] <- 0
  }
  # ensure it works for missing levels of response
  ind <- as.integer(which(table(dataSample[[name]]) > 0))
  if ( length(ind) > 2 && (nrow(grid)-length(exclude)) == 1 ) {
    probs <- t(probs)
  }
  # account for structural zeros
  if ( (!is.null(limit) || !is.null(censor)) && !is.null(dim(probs)) ) {
    if ( length(exclude) == 0 ) {
      probs <- adjustProbs(probs, grid, names(indGrid), limit, censor)
    } else {
      probs <- adjustProbs(probs, grid[-exclude, , drop=FALSE], names(indGrid)[-exclude], limit, censor)
    }
  }
  # local function for sampling from probabilities
  if ( length(ind) == 1 ) {
    resample <- function(k, n, p) rep.int(1, n[k])
  } else if ( length(ind) == 2 ) {
    resample <- function(k, n, p) spSample(n[k], c(1-p[k],p[k]))
  } else {
    resample <- function(k, n, p) spSample(n[k], p[k,])
  }
  # generate realizations for each combination
  sim <- as.list(rep.int(NA, length(indGrid)))
  if ( length(exclude) == 0 ) {
    ncomb <- as.integer(sapply(indGrid, length))
    sim <- lapply(1:length(ncomb), resample, ncomb, probs)
  } else {
    ncomb <- as.integer(sapply(indGrid[-exclude], length))
    sim[-exclude] <- lapply(1:length(ncomb), resample, ncomb, probs)
  }
  sim <- unsplit(sim, dataPop, drop=TRUE)
  # return realizations
  levels(response)[ind][sim]
}

generateValues_lm <- function(dataSample, dataPop, params) {
  if ( !nrow(dataSample) ) {
    return(numeric())
  }
  coef <- params$coef
  excludeLevels <- params$excludeLevels
  hasNewLevels <- params$hasNewLevels
  newLevels <- params$newLevels
  command <- params$command
  predNames <- params$predNames
  additional <- params$additional
  const <- params$const
  formula <- params$formula
  levels <- params$levels
  residuals <- params$residuals
  log <- params$log

  # fix: for each predictor, the level set must be equal in dataSample and dataPop
  for ( i in predNames ) {
    both <- intersect(levels(dataSample[[i]]), levels(dataPop[[i]]))
    a <- as.character(dataSample[[i]])
    a[!a%in%both] <- NA
    b <- as.character(dataPop[[i]])
    b[!b %in%both] <- NA
    dataSample[[i]] <- factor(a, levels=both)
    dataPop[[i]] <- factor(b, levels=both)
  }
  # unique combinations in the stratum of the population need to be computed for prediction
  indGrid <- split(1:nrow(dataPop), dataPop, drop=TRUE)
  grid <- dataPop[sapply(indGrid, function(i) i[1]), , drop=FALSE]
  grid <- as.data.frame(grid)

  # in sample, observations with NAs have been removed to fit the
  # model, hence population can have additional levels
  # these need to be removed since those probabilities cannot
  # be predicted from the model
  if ( excludeLevels ) {
    exclude <- mapply(function(pop, new) pop %in% new,
        pop=grid[, hasNewLevels, drop=FALSE], new=newLevels[hasNewLevels]
    )
    if ( is.null(dim(exclude)) ) {
      exclude <- which(any(exclude))
    } else {
      exclude <- which(apply(exclude, 1, any))
    }
    if ( length(exclude) > 0 ) {
      grid <- grid[-exclude, , drop=FALSE]
    }
    for ( j in predNames[hasNewLevels] ) {
      # drop new factor levels
      grid[, j] <- factor(as.character(grid[, j]), levels=levels(dataSample[[j]]))
    }
  } else {
    exclude <- integer()
  }
  # fit linear model
  mod <- eval(parse(text=command))
  # add coefficients from auxiliary model if necessary
  #tmp <- coef
  #coef[names(coef(mod))] <- coef(mod)
  #mod$coefficients <- coef
  # prediction
  # add 0 variable to combinations for use of 'model.matrix'
  newdata <- cbind(grid, 0)
  names(newdata) <- c(predNames, additional[1])
  newdata <- model.matrix(formula, data=newdata)

  if ( length(exclude) == 0 ) {
    pred <- spPredict(mod, newdata)
  } else {
    pred <- as.list(rep.int(NA, length(indGrid)))
    pred[-exclude] <- spPredict(mod, newdata)
  }
  pred <- unsplit(pred, dataPop, drop=TRUE)
  # add error terms
  if ( residuals ) {
    error <- sample(residuals(mod), size=nrow(dataPop), replace=TRUE)
  } else {
    mu <- median(residuals(mod))
    sigma <- mad(residuals(mod))
    error <- rnorm(nrow(dataPop), mean=mu, sd=sigma)
  }
  # return realizations
  sim <- pred + error
  if ( log ) {
    res <- exp(sim)  # transform back
    if ( !is.null(const) ) {
      res <- res - const  # subtract constant
    }
    return(res)
  } else {
    return(sim)
  }
}

generateValues_poisson <- function(dataSample, dataPop, params) {
  if ( !nrow(dataSample) ) {
    return(rep(0,nrow(dataPop)))
  }
  coef <- params$coef
  excludeLevels <- params$excludeLevels
  hasNewLevels <- params$hasNewLevels
  newLevels <- params$newLevels
  command <- params$command
  predNames <- params$predNames
  additional <- params$additional
  const <- params$const
  formula <- params$formula
  levels <- params$levels
  residuals <- params$residuals
  log <- params$log

  # fix: for each predictor, the level set must be equal in dataSample and dataPop
  for ( i in predNames ) {
    both <- intersect(levels(dataSample[[i]]), levels(dataPop[[i]]))
    a <- as.character(dataSample[[i]])
    a[!a%in%both] <- NA
    b <- as.character(dataPop[[i]])
    b[!b %in%both] <- NA
    dataSample[[i]] <- factor(a, levels=both)
    dataPop[[i]] <- factor(b, levels=both)
  }
  # unique combinations in the stratum of the population need to be computed for prediction
  indGrid <- split(1:nrow(dataPop), dataPop, drop=TRUE)
  grid <- dataPop[sapply(indGrid, function(i) i[1]), , drop=FALSE]
  grid <- as.data.frame(grid)

  # in sample, observations with NAs have been removed to fit the
  # model, hence population can have additional levels
  # these need to be removed since those probabilities cannot
  # be predicted from the model
  if ( excludeLevels ) {
    exclude <- mapply(function(pop, new) pop %in% new,
        pop=grid[, hasNewLevels, drop=FALSE], new=newLevels[hasNewLevels]
    )
    if ( is.null(dim(exclude)) ) {
      exclude <- which(any(exclude))
    } else {
      exclude <- which(apply(exclude, 1, any))
    }
    if ( length(exclude) > 0 ) {
      grid <- grid[-exclude, , drop=FALSE]
    }
    for ( j in predNames[hasNewLevels] ) {
      # drop new factor levels
      grid[, j] <- factor(as.character(grid[, j]), levels=levels(dataSample[[j]]))
    }
  } else {
    exclude <- integer()
  }
  # fit linear model
  mod <- eval(parse(text=command))
  # add coefficients from auxiliary model if necessary
  #tmp <- coef
  #coef[names(coef(mod))] <- coef(mod)
  #mod$coefficients <- coef
  # prediction
  # add 0 variable to combinations for use of 'model.matrix'
  newdata <- cbind(grid, 0)
  names(newdata) <- c(predNames, additional[1])

  if ( length(exclude) == 0 ) {
    pred <- round(predict(mod, newdata=newdata,type="response"))
  } else {
    pred <- as.list(rep.int(NA, length(indGrid)))
    pred[-exclude] <- round(predict(mod, newdata=newdata,type="response"))
  }
  pred <- unsplit(pred, dataPop, drop=TRUE)
  # add error terms
# addition of an error term, not implemented for Poisson Regression yet
#  if ( residuals ) {
#    error <- sample(residuals(mod), size=nrow(dataPop), replace=TRUE)
#  } else {
#    mu <- median(residuals(mod))
#    sigma <- mad(residuals(mod))
#    error <- rnorm(nrow(dataPop), mean=mu, sd=sigma)
#  }
  # return realizations
  return(pred)
  sim <- pred #+ error
  #return(sim)
}

generateValues_binary <- function(dataSample, dataPop, params) {
  excludeLevels <- params$excludeLevels
  hasNewLevels <- params$hasNewLevels
  newLevels <- params$newLevels
  predNames <- params$predNames
  name <- params$name
  weight <- params$weight
  useAux <- params$useAux
  tol <- params$tol
  eps <- params$eps
  if ( !nrow(dataSample) ) {
    return(numeric())
  }
  # if all y values are the same return the same value for everybody
  if(length(unique(dataSample[[name]]))==1){
    return(rep(dataSample[[name]][1],nrow(dataPop)))
  }
  # unique combinations in the stratum of the population need to be computed for prediction
  indGrid <- split(1:nrow(dataPop), dataPop, drop=TRUE)
  grid <- dataPop[sapply(indGrid, function(i) i[1]), , drop=FALSE]
  grid <- as.data.frame(grid)

  # in sample, observations with NAs have been removed to fit the
  # model, hence population can have additional levels
  # these need to be removed since those probabilities cannot
  # be predicted from the model
  if ( excludeLevels ) {
    exclude <- mapply(function(pop, new) pop %in% new,
      pop=grid[, hasNewLevels, drop=FALSE],
      new=newLevels[hasNewLevels]
    )
    if ( is.null(dim(exclude)) ) {
      exclude <- which(any(exclude))
    } else {
      exclude <- which(apply(exclude, 1, any))
    }
    if ( length(exclude) > 0 ) {
      grid <- grid[exclude, , drop=FALSE]
    }
    for ( j in predNames[hasNewLevels] ) {
      # drop new factor levels
      grid[, j] <- factor(as.character(grid[, j]), levels=levels(dataSample[[j]]))
    }
  } else {
    exclude <- integer()
  }
  # add 0 variable to combinations for use of 'model.matrix'
  Xnew <- cbind(grid, 0)
  names(Xnew) <- c(predNames, name)
  Xnew <- model.matrix(params$command, data=Xnew)

  # fit logit model
  X <- model.matrix(params$command, data=dataSample)
  y <- dataSample[[name]]
  weights <- dataSample[[weight]]
  mod <- logitreg(X, y, weights=weights)
  # add parameters from auxiliary model if necessary
  if ( useAux ) {
    indPar <- abs(mod$par) < tol
    mod$par[indPar] <- params$par[indPar]

    # remove non-existing combinations from mod$par
    # reason: auxiliary model is estimated on total population
    ii <- setdiff(names(mod$par), colnames(Xnew))
    if ( length(ii) >0 ) {
      mod$par <- mod$par[!names(mod$par)%in%ii]
    }
  }

  # predict probabilities
  tmp <- exp(Xnew %*% mod$par)
  # avoid integer overflow
  p <- ifelse(is.infinite(tmp), 1, as.numeric(tmp / (1 + tmp)))
  # set too small probabilities to exactly 0
  if ( !is.null(eps) ) {
    p[p < eps] <- 0
  }
  # generate realizations for each combination
  if ( length(exclude) == 0 ) {
    ncomb <- as.integer(sapply(indGrid, length))
    sim <- lapply(1:length(ncomb), function(k) {
      spSample(ncomb[k], c(1-p[k], p[k])) - 1
    })
  } else {
    ncomb <- as.integer(sapply(indGrid[-exclude], length))
    sim <- as.list(rep.int(NA, length(indGrid)))
    sim[-exclude] <- lapply(1:length(ncomb), function(k) {
      spSample(ncomb[k], c(1-p[k], p[k])) - 1
    })
  }
  # return realizations
  unsplit(sim, dataPop, drop=TRUE)
}

genVals <- function(dataSample, dataPop, params, typ) {
  # unify level-set of predictors
  for ( i in params$predNames ) {
    dataSample[[i]] <- cleanFactor(dataSample[[i]])
    dataPop[[i]] <- cleanFactor(dataPop[[i]])
  }

  if ( !typ %in% c("multinom","lm","binary","poisson") ) {
    stop("unsupported value for argument 'type' in genVals()\n")
  }
  if ( typ=="binary") {
    res <- generateValues_binary(dataSample, dataPop, params)
  }else if ( typ=="lm") {
    res <- generateValues_lm(dataSample, dataPop, params)
  }else if ( typ=="multinom" ) {
    res <- generateValues_multinom(dataSample, dataPop, params)
  }else if ( typ=="poisson") {
    res <- generateValues_poisson(dataSample, dataPop, params)
  }
  res
}

runModel <- function(dataS, dataP, params, typ) {
  x <- NULL
  strata <- params$strata
  pp <- parallelParameters(nr_cpus=params$nr_cpus, nr_strata=length(levels(dataS[[strata]])))
  indStrata <- params$indStrata
  predNames <- params$predNames
  additional <- unique(c(params$additional, params$name, params$weight))
  if ( pp$parallel ) {
    # windows
    if ( pp$have_win ) {
      cl <- makePSOCKcluster(pp$nr_cores)
      registerDoParallel(cl,cores=pp$nr_cores)
      valuesCat <- foreach(x=levels(dataS[[strata]]), .options.snow=list(preschedule=TRUE)) %dopar% {
        genVals(
          dataSample=dataS[dataS[[strata]] == x,],
          dataPop=dataP[indStrata[[x]], predNames, with=F],
          params,
          typ=typ)
      }
      stopCluster(cl)
    }
    # linux/mac
    if ( !pp$have_win ) {
      valuesCat <- mclapply(levels(dataS[[strata]]), function(x) {
        genVals(
          dataSample=dataS[dataS[[strata]] == x,],
          dataPop=dataP[indStrata[[x]], predNames, with=F],
          params=params,
          typ=typ)
      },mc.cores=pp$nr_cores)
    }
  } else {
    valuesCat <- lapply(levels(dataS[[strata]]), function(x) {
       if(params$verbose){
         cat("Current by group for the binary model:",x,"\n")
       }
      genVals(
        dataSample=dataS[dataS[[strata]] == x,c(predNames, additional), with=F],
        dataPop=dataP[indStrata[[x]], predNames, with=F],
        params=params,
        typ=typ)
    })
  }

  # check for errors
  res <- sapply(valuesCat, class)
  if ( any(res=="try-error") ) {
    stop(paste0("Error in estimating the linear model. Try to specify a more simple model!\n"))
  }

  if ( typ=="multinom" ) {
    response <- dataS[[params$name]]
    valuesCat <- factor(unsplit(valuesCat, dataP[[strata]]), levels=levels(response))
  }
  if ( typ=="binary" ) {
    valuesCat <- unsplit(valuesCat, dataP[[strata]], drop=FALSE)
  }
  if ( typ%in%c("poisson","lm") ) {
    valuesCat <- unsplit(valuesCat, dataP[[strata]], drop=FALSE)
  }
  return(valuesCat)
}


#' Simulate continuous variables of population data
#'
#' Simulate continuous variables of population data using multinomial
#' log-linear models combined with random draws from the resulting categories
#' or (two-step) regression models combined with random error terms. The
#' household structure of the population data and any other categorical
#' predictors need to be simulated beforehand.
#'
#' If \code{method} is \code{"lm"}, the behavior for two-step models is
#' described in the following.
#'
#' If \code{zeros} is \code{TRUE} and \code{log} is not \code{TRUE} or the
#' variable specified by \code{additional} does not contain negative values, a
#' log-linear model is used to predict whether an observation is zero or not.
#' Then a linear model is used to predict the non-zero values.
#'
#' If \code{zeros} is \code{TRUE}, \code{log} is \code{TRUE} and \code{const}
#' is specified, again a log-linear model is used to predict whether an
#' observation is zero or not. In the linear model to predict the non-zero
#' values, \code{const} is added to the variable specified by \code{additional}
#' before the logarithms are taken.
#'
#' If \code{zeros} is \code{TRUE}, \code{log} is \code{TRUE}, \code{const} is
#' \code{NULL} and there are negative values, a multinomial log-linear model is
#' used to predict negative, zero and positive observations. Categories for the
#' negative values are thereby defined by \code{breaks}. In the second step, a
#' linear model is used to predict the positive values and negative values are
#' drawn from uniform distributions in the respective classes.
#'
#' If \code{zeros} is \code{FALSE}, \code{log} is \code{TRUE} and \code{const}
#' is \code{NULL}, a two-step model is used if there are non-positive values in
#' the variable specified by \code{additional}. Whether a log-linear or a
#' multinomial log-linear model is used depends on the number of categories to
#' be used for the non-positive values, as defined by \code{breaks}. Again,
#' positive values are then predicted with a linear model and non-positive
#' values are drawn from uniform distributions.
#'
#' The number of cpus are selected automatically in the following manner. The
#' number of cpus is equal the number of strata. However, if the number of cpus
#' is less than the number of strata, the number of cpus - 1 is used by
#' default. This should be the best strategy, but the user can also overwrite
#' this decision.
#'
#' @name simContinuous
#' @param simPopObj a \code{\linkS4class{simPopObj}} holding household survey
#' data, population data and optionally some margins.
#' @param additional a character string specifying the additional continuous
#' variable of \code{dataS} that should be simulated for the population data.
#' Currently, only one additional variable can be simulated at a time.
#' @param method a character string specifying the method to be used for
#' simulating the continuous variable. Accepted values are \code{"multinom"},
#' for using multinomial log-linear models combined with random draws from the
#' resulting categories, \code{"lm"}, for using (two-step) regression
#' models combined with random error terms and \code{"poisson"} for using Poisson regression for count variables.
#' @param zeros a logical indicating whether the variable specified by
#' \code{additional} is semi-continuous, i.e., contains a considerable amount
#' of zeros. If \code{TRUE} and \code{method} is \code{"multinom"}, a separate
#' factor level for zeros in the response is used. If \code{TRUE} and
#' \code{method} is \code{"lm"}, a two-step model is applied. The first step
#' thereby uses a log-linear or multinomial log-linear model (see
#' \dQuote{Details}).
#' @param breaks an optional numeric vector; if multinomial models are
#' computed, this can be used to supply two or more break points for
#' categorizing the variable specified by \code{additional}. If \code{NULL},
#' break points are computed using weighted quantiles.
#' @param lower,upper optional numeric values; if multinomial models are
#' computed and \code{breaks} is \code{NULL}, these can be used to specify
#' lower and upper bounds other than minimum and maximum, respectively. Note
#' that if \code{method} is \code{"multinom"} and \code{gpd} is \code{TRUE}
#' (see below), \code{upper} defaults to \code{Inf}.
#' @param equidist logical; if \code{method} is \code{"multinom"} and
#' \code{breaks} is \code{NULL}, this indicates whether the (positive) default
#' break points should be equidistant or whether there should be refinements in
#' the lower and upper tail (see \code{\link{getBreaks}}).
#' @param probs numeric vector with values in \eqn{[0, 1]}; if \code{method} is
#' \code{"multinom"} and \code{breaks} is \code{NULL}, this gives probabilities
#' for quantiles to be used as (positive) break points. If supplied, this is
#' preferred over \code{equidist}.
#' @param gpd logical; if \code{method} is \code{"multinom"}, this indicates
#' whether the upper tail of the variable specified by \code{additional} should
#' be simulated by random draws from a (truncated) generalized Pareto
#' distribution rather than a uniform distribution.
#' @param threshold a numeric value; if \code{method} is \code{"multinom"},
#' values for categories above \code{threshold} are drawn from a (truncated)
#' generalized Pareto distribution.
#' @param est a character string; if \code{method} is \code{"multinom"}, the
#' estimator to be used to fit the generalized Pareto distribution.
#' @param limit an optional named list of lists; if multinomial models are
#' computed, this can be used to account for structural zeros. The names of the
#' list components specify the predictor variables for which to limit the
#' possible outcomes of the response. For each predictor, a list containing the
#' possible outcomes of the response for each category of the predictor can be
#' supplied. The probabilities of other outcomes conditional on combinations
#' that contain the specified categories of the supplied predictors are set to
#' 0. Currently, this is only implemented for more than two categories in the
#' response.
#' @param censor an optional named list of lists or \code{data.frame}s; if
#' multinomial models are computed, this can be used to account for structural
#' zeros. The names of the list components specify the categories that should
#' be censored. For each of these categories, a list or \code{data.frame}
#' containing levels of the predictor variables can be supplied. The
#' probability of the specified categories is set to 0 for the respective
#' predictor levels. Currently, this is only implemented for more than two
#' categories in the response.
#' @param log logical; if \code{method} is \code{"lm"}, this indicates whether
#' the linear model should be fitted to the logarithms of the variable
#' specified by \code{additional}. The predicted values are then
#' back-transformed with the exponential function. See \dQuote{Details} for
#' more information.
#' @param const numeric; if \code{method} is \code{"lm"} and \code{log} is
#' \code{TRUE}, this gives a constant to be added before log transformation.
#' @param alpha numeric; if \code{method} is \code{"lm"}, this gives trimming
#' parameters for the sample data. Trimming is thereby done with respect to the
#' variable specified by \code{additional}. If a numeric vector of length two
#' is supplied, the first element gives the trimming proportion for the lower
#' part and the second element the trimming proportion for the upper part. If a
#' single numeric is supplied, it is used for both. With \code{NULL}, trimming
#' is suppressed.
#' @param residuals logical; if \code{method} is \code{"lm"}, this indicates
#' whether the random error terms should be obtained by draws from the
#' residuals. If \code{FALSE}, they are drawn from a normal distribution
#' (median and MAD of the residuals are used as parameters).
#' @param keep logical; if multinomial models are computed, this indicates
#' whether the simulated categories should be stored as a variable in the
#' resulting population data. If \code{TRUE}, the corresponding column name is
#' given by \code{additional} with postfix \code{"Cat"}.
#' @param maxit,MaxNWts control parameters to be passed to
#' \code{\link[nnet]{multinom}} and \code{\link[nnet]{nnet}}. See the help file
#' for \code{\link[nnet]{nnet}}.
#' @param tol if \code{method} is \code{"lm"} and \code{zeros} is \code{TRUE},
#' a small positive numeric value or \code{NULL}. When fitting a log-linear
#' model within a stratum, factor levels may not exist in the sample but are
#' likely to exist in the population. However, the coefficient for such factor
#' levels will be 0. Therefore, coefficients smaller than \code{tol} in
#' absolute value are replaced by coefficients from an auxiliary model that is
#' fit to the whole sample. If \code{NULL}, no auxiliary log-linear model is
#' computed and no coefficients are replaced.
#' @param nr_cpus if specified, an integer number defining the number of cpus
#' that should be used for parallel processing.
#' @param eps a small positive numeric value, or \code{NULL} (the default). In
#' the former case and if (multinomial) log-linear models are computed,
#' estimated probabilities smaller than this are assumed to result from
#' structural zeros and are set to exactly 0.
#' @param regModel allows to specify the model that should be for the
#' simulation of the additional continuous variable. The following choices are
#' possible: \itemize{ \item'basic'only the basic household-variables
#' (generated with \code{\link{simStructure}}) are used.  \item'available'all
#' available variables (that are common in the sample and the syntetic
#' population (e.g. previously generated variables) are used for the modeling.
#' Should be used with care because all variables are automatically used as
#' factors!  \item formula-object: Users may also specify a specific formula
#' (class 'formula') that will be used. Checks are performed that all required
#' variables are available.}
#' @param byHousehold if NULL, simulated values are used as is. If either \code{'sum'},
#' \code{'mean'} or \code{'random'} is specified, the values are aggregated and each member
#' of the household gets the same value (mean, sum or a random value) assigned.
#' @param imputeMissings if TRUE, missing values in variables that are used for
#' the underlying model are imputed using hock-deck.
#' @param seed optional; an integer value to be used as the seed of the random
#' number generator, or an integer vector containing the state of the random
#' number generator to be restored.
#' @param verbose (logical) if \code{TRUE}, additional output is written to the promt
#' @param by defining which variable to use as split up variable of the estimation. Defaults to the strata variable.
#' @return An object of class \code{\linkS4class{simPopObj}} containing survey
#' data as well as the simulated population data including the continuous
#' variable specified by \code{additional} and possibly simulated categories
#' for the desired continous variable.
#' @note The basic household structure and any other categorical predictors
#' need to be simulated beforehand with the functions
#' \code{\link{simStructure}} and \code{\link{simCategorical}}, respectively.
#' @author Bernhard Meindl and Andreas Alfons (based on code by Stefan Kraft)
#' @seealso \code{\link{simStructure}}, \code{\link{simCategorical}},
#' \code{\link{simComponents}}, \code{\link{simEUSILC}}
#' @keywords datagen
#' @export
#' @examples
#'
#' data(eusilcS)
#' \dontrun{
#' ## approx. 20 seconds computation time
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' simPop <- simStructure(data=inp, method="direct",
#'   basicHHvars=c("age", "rb090", "hsize", "pl030", "pb220a"))
#'
#' regModel = ~rb090+hsize+pl030+pb220a
#'
#' # multinomial model with random draws
#' eusilcM <- simContinuous(simPop, additional="netIncome",
#'               regModel = regModel,
#'               upper=200000, equidist=FALSE, nr_cpus=1)
#' class(eusilcM)
#' }
#'
#' \dontrun{
#' # two-step regression
#' eusilcT <- simContinuous(simPop, additional="netIncome",
#'               regModel = "basic",
#'               method = "lm")
#' class(eusilcT)
#' }
#'
simContinuous <- function(simPopObj, additional = "netIncome",
  method = c("multinom", "lm","poisson"), zeros = TRUE,
  breaks = NULL, lower = NULL, upper = NULL,
  equidist = TRUE, probs = NULL, gpd = TRUE,
  threshold = NULL, est = "moments", limit = NULL,
  censor = NULL, log = TRUE, const = NULL,
  alpha = 0.01, residuals = TRUE, keep = TRUE,
  maxit = 500, MaxNWts = 1500,
  tol = .Machine$double.eps^0.5,
  nr_cpus=NULL, eps = NULL, regModel="basic", byHousehold=NULL,
  imputeMissings=FALSE, seed, verbose=FALSE,by="strata") {

  x <- hhid <- vals <- id <- V1 <- randId <- NULL

  if ( !is.null(byHousehold) ) {
    if ( !byHousehold %in% c("mean","sum","random") ) {
      stop("invalid value for argument 'byHousehold'. Allowed values are 'mean', 'sum' or 'random'!\n")
    }
  }

  samp <- simPopObj@sample
  pop <- simPopObj@pop
  basic <- simPopObj@basicHHvars
  if(by=="strata"){
    strata <- samp@strata
  }else if(!is.null(by)){
    strata <- by
  }
  weight <- samp@weight

  dataS <- samp@data
  dataP <- pop@data

  if ( additional %in% names(dataP)) {
    stop(paste0("Variable '",additional,"' already available in the synthetic population!\n"))
  }

  ## initializations
  if ( !missing(seed) ) {
    set.seed(seed)
  }

  if ( length(additional) != 1 ) {
    stop("currently only one additional variable can be generated at a time")
  }
  if ( !additional %in% colnames(samp@data) ) {
    stop("variable 'additional' must be included in the sample of input 'simPopObj'!\n")
  }

  regInput <- regressionInput(simPopObj, additional=additional, regModel=regModel)
  predNames <- regInput[[1]]$predNames
  estimationModel <- regInput[[1]]$formula

  varNames <- unique(c(predNames, weight, additional, strata))
  if(!strata%in%colnames(dataS)){
    stop(strata," is defined as by variable, but not in the sample data set.")
  }
  dataS <- dataS[,varNames, with=F]

  method <- match.arg(method)
  zeros <- isTRUE(zeros)
  log <- isTRUE(log)
  if(log&&method=="poisson"){
    log <- FALSE
    warning("For Poisson regression the log=TRUE parameter is ignored and the numeric variable is not transformed.")
  }
  if(!is.null(alpha)&&method=="poisson"){
    alpha <- NULL
    warning("For Poisson regression the alpha!=NULL is not yet implemented and therefore set to NULL.")
  }
  if ( is.numeric(alpha) && length(alpha) > 0 ) {
    alpha <- rep(alpha, length.out=2)
    if ( !all(is.finite(alpha)) || any(alpha < 0) || sum(alpha) >= 1 ) {
      alpha <- NULL
      warning("invalid parameter 'alpha': trimming is not applied\n")
    }
  } else {
    alpha <- NULL
  }

  # observations with missings are excluded from simulation
  #exclude <- getExclude(dataS[,c(additional,predNames),with=F]) # fixes #31?
  #if ( length(exclude) ) {
  #  dataS <- dataS[-exclude,]
  #}

  # temporarily impute (using hotdeck) / or check (if imputeMissings=FALSE)
  # missing values in additional variables in the sample
  if ( is.null(strata) ) {
    modelVars <- setdiff(predNames, c(weight,basic,pop@hhsize))
  } else {
    modelVars <- setdiff(predNames, c(strata,weight,basic,pop@hhsize))
  }

  if ( length(modelVars) > 0 & imputeMissings ) {
    dataS_orig <- dataS[,modelVars,with=F]
    dataS <- hotdeck(dataS, variable=modelVars, domain_var=strata, imp_var=FALSE)
  }

  # check for NAs and warn user
  if ( !imputeMissings) {
    naTab <- dataS[,lapply(.SD, is.na), .SDcols=c(additional,predNames)]
    perc.miss <- sum(rowSums(naTab)!=0) / nrow(dataS)
    if ( perc.miss > 0 ) {
      wm <- paste0("There are ~",formatC(100*perc.miss,format="f", digits=1),"% ")
      wm <- paste0(wm, "observations in the response/predictors with at least one missing variable.\n")
      wm <- paste0(wm, "If you get errors in the estimation procedure, consider to recode these missing ")
      wm <- paste0(wm, "values (e.g. by assigning an additional category) or try to specify a different model.\n\n")
      for ( z in 1:ncol(naTab) ) {
        vv <- colnames(naTab)[z]
        missv <- sum(naTab[[z]])
        missp <- formatC(100*missv/nrow(dataS),format="f", digits=1)
        if ( vv == additional ) {
          wm <- paste0(wm, "Variable '",vv,"' (response): ",missv," missing values (~",missp,"%).\n")
        } else {
          wm <- paste0(wm, "Variable '",vv,"' (predictor): ",missv," missing values (~",missp,"%).\n")
        }
      }
      if(verbose) warning(wm)
    }
  }

  # variables are coerced to factors
  select <- unique(c(predNames, strata)) # strata always included
  dataS <- checkFactor(dataS, select)
  if(!strata%in%colnames(dataP)){
    stop(strata," is defined as by variable, but not in the population data set.")
  }
  dataP <- checkFactor(dataP, select)

  # sample data of variable to be simulated
  additionalS <- dataS[[additional]]

  ## determine which models to fit and do further initializations
  haveBreaks <- !is.null(breaks)
  if ( method == "multinom" ) {
    useMultinom <- TRUE
    useLogit <- FALSE
    useLm <- FALSE
    usePoisson <- FALSE
    # define break points (if missing)
    if ( haveBreaks ) {
      checkBreaks(breaks)
      breaks <- if(zeros) union(breaks, 0) else unique(breaks)
      breaks <- sort(breaks)
    } else {
      if ( is.null(upper) && gpd ) {
        upper <- Inf
      }
      breaks <- getBreaks(additionalS, dataS[[weight]], zeros, lower, upper, equidist, probs)
    }
  } else {
    if(method=="lm"){
      useLm <- TRUE
      usePoisson <- FALSE
    }else if(method=="poisson"){
      useLm <- FALSE
      usePoisson <- TRUE
    }

    if ( log ) {
      if ( is.null(const) ) {
        ## use log-transformation
        # check for negative values
        neg <- which(additionalS < 0)
        haveNeg <- length(neg) > 0
        if ( haveNeg ) {
          # define break points for negative values
          if ( haveBreaks ) {
            checkBreaks(breaks)
            breaks <- c(unique(breaks[breaks < 0]), 0)
          } else {
            breaks <- getBreaks(additionalS[neg], dataS[[weight]][neg], zeros=TRUE, lower, upper)
          }
          if ( zeros || length(breaks) > 2 ) {
            useMultinom <- TRUE
            breaks <- c(breaks, Inf)  # add Inf to breakpoints
          } else {
            useMultinom <- FALSE
          }
          useLogit <- !useMultinom
        } else {
          useLogit <- zeros || any(additionalS == 0)
          useMultinom <- FALSE
        }
      } else {
        # check constant
        if ( !is.numeric(const) || length(const) == 0 ) {
          stop("'const' must be numeric\n")
        } else {
          const <- const[1]
        }
        # set control parameters
        useLogit <- zeros || any(additionalS == 0)
        useMultinom <- FALSE
      }
    } else {
      # logistic model is used in case of semi-continuous variable
      useLogit <- zeros
      # multinomial model is not needed
      useMultinom <- FALSE
    }
  }

  ## some general preparations for the simulation
  # list indStrata contains the indices of dataP split by strata
  N <- nrow(dataP)
  indP <- 1:N
  indStrata <- split(indP, dataP[[strata]])
  #fpred <- paste(predNames, collapse = " + ")  # for formula
  # check if population data contains factor levels that do not exist
  # in the sample
  newLevels <- lapply(predNames, function(nam) {
    levelsS <- levels(dataS[[nam]])
    levelsP <- levels(dataP[[nam]])
    levelsP[!(levelsP %in% levelsS)]
  })
  hasNewLevels <- sapply(newLevels, length) > 0
  excludeLevels <- any(hasNewLevels)

  ## preparations for multinomial or binomial logit model
  if ( useMultinom || useLogit ) {
    name <- getCatName(additional)
    estimationModel <- gsub(additional, name, estimationModel)
    # remove strata variable from estimation model if we are computing on multiple cores
    # else multinom fails because the variable has only one factor!
    if ( !is.null(dataS[[strata]]) ) {
      if ( parallelParameters(nr_cpus, length(levels(dataS[[strata]])))$nr_cores > 1 ) {
        estimationModel <- gsub(paste0("[+]",strata),"",estimationModel)
      }
    }
  }

  if ( useMultinom ) {
    ## some preparations
    dataS[[name]] <- getCat(additionalS, breaks, zeros, right=TRUE)
    response <- dataS[[name]]  # response variable
    # check threshold for GPD (if supplied)
    if ( !useLm && gpd && !is.null(threshold) && length(threshold) != 1 ) {
      stop("'threshold' must be a single numeric value")
    }

    ## simulate categories
    # TODO: share code with 'simCategorical'
    params <- list()
    params$excludeLevels <- excludeLevels
    # command needs to be constructed as string
    # this is actually a pretty ugly way of fitting the model
    params$command <- paste("suppressWarnings(multinom(", estimationModel,
      ", weights=", weight, ", data=dataSample, trace=FALSE",
      ", maxit=maxit, MaxNWts=MaxNWts))", sep="")
    params$maxit <- maxit
    params$MaxNWts <- MaxNWts
    params$eps <- eps
    params$limit <- limit
    params$censor <- censor
    params$hasNewLevels  <- hasNewLevels
    params$newLevels <- newLevels
    params$name <- name
    params$response <- response
    params$strata <- strata
    params$nr_cpus <- nr_cpus
    params$indStrata <- indStrata
    params$predNames <- predNames
    params$additional <- c(additional, weight)
    params$verbose <- verbose
    if(verbose) cat("running multinom with the following model:\n")
    if(verbose) cat(gsub("))",")",gsub("suppressWarnings[(]","",params$command)),"\n")

    # run in parallel if possible
    valuesCat <- runModel(dataS, dataP, params, typ="multinom")

    ## simulate (semi-)continuous values
    tcat <- table(valuesCat)
    ncat <- length(tcat)

    icat <- 1:ncat
    values <- as.list(rep.int(NA, ncat))
    # zeros

    if ( zeros ) {
      # bug: missing 0 even though zeros is not null?
      izero <- which(breaks == 0)
      values[izero] <- 0
      tcat <- tcat[-izero]
      ncat <- length(tcat)
      icat <- icat[-izero]
    }

    # values to be simulated with linear model or draws from Pareto
    # distribution
    if ( useLm ) {
      # last breakpoint is Inf, the one before is 0
      nunif <- ncat - 1  # leave category of positive values out
    } else {
      nbreaks <- length(breaks)
      if ( gpd ) {
        if ( is.null(threshold) ) {
          if ( !haveBreaks && (!isTRUE(equidist) || !is.null(probs)) ) {
            ngpd <- nbreaks-2
          } else {
            ngpd <- nbreaks-1
          }
        } else if ( any(tmp <- breaks >= threshold) ) {
          ngpd <- min(which(tmp))
        } else {
          ngpd <- nbreaks
        }
      } else {
        ngpd <- nbreaks
      }
      if ( gpd && ngpd <= ncat ) {
        # adjust threshold and fit GPD
        threshold <- breaks[ngpd]  # adjust threshold
        estPar <- fitgpd(additionalS, threshold, est)  # fit GPD
        estPar <- estPar[["fitted.values"]]  # parameters of GPD
        # generalized pareto distribution
        igpd <- ngpd:ncat
        values[icat[igpd]] <- lapply(igpd, function(i) {
          truncPareto(tcat[i], loc=threshold, scale=estPar["scale"], shape=estPar["shape"], breaks[i], breaks[i+1])
        })
      }
      nunif <- ngpd - 1
    }
    # uniform distribution
    if ( nunif > 0 ) {
      iunif <- 1:nunif
      values[icat[iunif]] <- lapply(iunif, function(i) {
        runif(tcat[i], breaks[i], breaks[i+1])
      })
    }
    # turn list into vector of values
    values <- unsplit(values, valuesCat)
  }

  if ( useLogit ) {
    ## some preparations
    if ( log && is.null(const) && haveNeg ) {
      indS <- additionalS > 0
    } else {
      indS <- additionalS != 0
    }
    dataS[[name]] <- as.integer(indS)
    estimationModel <- as.formula(estimationModel)  # formula for model
    # auxiliary model for all strata (used in case of empty combinations)
    useAux <- !is.null(tol)
    if(method=="poisson"){
      useAux <- FALSE
    }
    if ( useAux ) {
      if ( length(tol) != 1 || tol <= 0 ) {
        stop("'tol' must be a single small positive value!\n")
      }
      #nas <- sum(!complete.cases(dataS))
      #if ( length(nas) > 0 ) {
      #  warning("\nwe Hotdeck-imputation of missing values sample data required!\n")
      #  dataS <- hotdeck(dataS, variable=predNames, domain_var=samp@strata, imp_var=FALSE)
      #}
      X <- model.matrix(estimationModel, data=dataS)
      y <- dataS[[name]]
      weights <- dataS[[weight]]
      mod <- logitreg(X, y, weights=weights)
      par <- mod$par
    } else {
      par <- NULL
    }

    ## simulate binary vector
    params <- list()
    params$excludeLevels <- excludeLevels
    params$hasNewLevels <- hasNewLevels
    params$newLevels <- newLevels
    params$predNames <- predNames
    params$tol <- tol
    params$eps <- eps
    params$weight <- weight
    params$useAux <- useAux
    params$name <- name
    params$strata <- strata
    params$nr_cpus <- nr_cpus
    params$indStrata <- indStrata
    params$predNames <- predNames
    params$additional <- additional
    params$par <- par
    params$command <- estimationModel
    params$verbose <- verbose
    # run in parallel if possible

    valuesCat <- runModel(dataS, dataP, params, typ="binary")
  }

  if ( useLm || usePoisson ) {
    ## some preparations
    if ( useMultinom ) {
      catLm <- names(tcat)[ncat]  # category for positive values
      dataS <- dataS[response == catLm, , drop=FALSE]
      indP <- valuesCat == catLm
    } else if( useLogit ) {
      dataS <- dataS[indS]  # select only non-zeros
      indP <- valuesCat == 1  # indicates non-zeros in population
    }
    if ( useMultinom || useLogit ) {
      # adjust population data
      if ( any(indP) ) {
        dataPop <- dataP[indP]
      } else {
        dataPop <- dataP
      }
      # list indStrata is adjusted so that it only contains
      # indices of persons in population with non-zero value
      indStrata <- split(1:nrow(dataPop), dataPop[[strata]])
    } else {
      dataPop <- dataP
    }

    ## trim data (if specified)
    if ( !is.null(alpha) ) {
      additional <- additional[1]
      additionalS <- dataS[[additional]]
      p <- c(alpha[1], 1-alpha[2])
      bounds <- quantileWt(additionalS, dataS[[weight]], p)
      select <- additionalS > bounds[1] & additionalS < bounds[2]
      dataSample <- dataS[select, , drop=FALSE]
      # check if all relevant levels of predictor variables are still
      # contained in sample after trimming
      # if not, trimming is not applied and a warning message is generated
      check <- unlist(sapply(predNames, function(i) {
        table(dataS[[i]]) > 0 & table(dataSample[[i]]) == 0
      }))
      if ( any(check) ) {
        dataSample <- dataS
        warning("trimming could not be applied\n")
      }
    } else {
      dataSample <- dataS
    }

    ## fit linear model
    # formula for linear model
    if ( log ) {
      fname <- paste("log(", additional, if(!is.null(const)) " + const", ")", sep = "")
    } else {
      fname <- additional
    }
    fstring <- paste0(fname, " ~ ", tail(unlist(strsplit(as.character(estimationModel),"~")),1))
    formula <- as.formula(fstring)
    # auxiliary model for all strata (used in case of empty combinations)
    weights <- dataSample[[weight]]
    if(useLm){
      mod <- lm(formula, weights=weights, data=dataSample)
      coef <- coef(mod)
    }else if(usePoisson){
      mod <- glm(formula, weights=weights, data=dataSample,family=poisson())
      coef <- coef(mod)
    }


    # simulate values
    params <- list()
    params$coef <- coef
    if(useLm){
      params$command <- paste("lm(", fstring,", weights=", weight, ", data=dataSample)", sep="")
    }else if(usePoisson){
      params$command <- paste("glm(", fstring,", weights=", weight, ", data=dataSample,family=poisson())", sep="")
    }
    #params$name <- fname
    params$name <- additional
    params$excludeLevels <- excludeLevels
    params$hasNewLevels <- hasNewLevels
    params$newLevels <- newLevels
    params$predNames <- predNames
    params$additional <- c(additional, weight)
    params$const <- const
    params$formula <- formula
    params$residuals <- residuals
    params$log <- log
    params$strata <- strata
    params$nr_cpus <- nr_cpus
    params$indStrata <- indStrata
    params$predNames <- predNames
    params$verbose <- verbose
    if(useLm){
      valuesTmp <- runModel(dataSample, dataPop, params, typ="lm")
    }else if(usePoisson){
      valuesTmp <- runModel(dataSample, dataPop, params, typ="poisson")
    }

    ## put simulated values together
    if ( useMultinom ) {
      values[which(indP == 1)] <- valuesTmp
    } else {
      if ( useLogit ) {
        if ( log && is.null(const) && haveNeg ) {
          # only one category for non-positive values (two breakpoints, one of them is 0)
          values <- rep.int(NA, N)
          values[indP] <- runif(sum(indP), breaks[1], breaks[2])
        } else {
          values <- ifelse(is.na(indP), NA, 0) # only zeros
        }
        values[indP] <- valuesTmp
      } else {
        values <- valuesTmp
      }
    }
  }

  # reset imputed variables in sample
  if ( imputeMissings ) {
    for ( i in 1:ncol(dataS_orig)) {
      cmd <- paste0("dataS[,",colnames(dataS_orig)[i],":=dataS_orig$",colnames(dataS_orig)[i],"]")
      eval(parse(text=cmd))
    }
  }

  # attach new variable(s) to population data
  if ( useMultinom && keep ) {
    dataP[[name]] <- valuesCat
  }

  # calculate mean of new variable by household
  if ( !is.null(byHousehold) ) {
    xx <- data.table(id=1:length(values), hhid=dataP[[pop@hhid]], vals=values)
    setkey(xx, hhid)

    if ( byHousehold=="mean" ) {
      yy <- xx[,mean(vals, na.rm=TRUE), by=key(xx)]
    }
    if ( byHousehold=="sum" ) {
      yy <- xx[,sum(vals, na.rm=TRUE), by=key(xx)]
    }
    if ( byHousehold=="random" ) {
      zz <- xx[,.N,by=key(xx)]
      ids1 <- zz[N==1]
      yy <- xx[hhid%in%ids1$hhid]
      yy[,id:=NULL]
      ids2 <- zz[N>1]
      if ( nrow(ids2) > 0 ) {
        xx2 <- xx[hhid %in% ids2$hhid]
        xx2[,randId:=sample(1:nrow(xx2))]
        setkey(xx2, hhid, randId)
        setkey(xx2, hhid)
        yy2 <- unique(xx2)
        yy2[,c("id","randId"):=NULL]
        yy <- rbind(yy, yy2)
      }
      setkey(yy, hhid)
      yy[,V1:=vals]
      yy[,vals:=NULL]
    }
    xx <- merge(xx, yy, all.x=TRUE)
    setkey(xx, id)
    xx[is.nan(V1), V1:=NA]
    values <- xx$V1
  }

  # return simulated data
  dataP[[additional]] <- values
  simPopObj@pop@data <- dataP
  invisible(simPopObj)
}
