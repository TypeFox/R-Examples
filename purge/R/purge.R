#' Purge training data from a model
#'
#' Most R model implementations store the training data
#' within the fitted object, often many times. It can be useful
#' to remove the embedded data for portability, especially
#' if the only required functionality is to predict on new data.
#'
#' @importFrom methods as
#' @param model A fitted R model object
#' @return A fitted R model object, purged of its training data,
#'         but retaining its predict functionality on new data
#' @examples
#' x <- rnorm(1000)
#' y <- x + rnorm(1000)
#' unpurged.model <- lm(y ~ x)
#' purged.model <- purge(unpurged.model)
#' object.size(unpurged.model)
#' object.size(purged.model)
#'
#' @export
purge <- function(model) {
  UseMethod('purge')
}

#' @describeIn purge Default purge returns a copy of the model
#' @export
purge.default <- function(model) {
  model.copy <- model
  return(model.copy)
}

#' @describeIn purge Purges a glm model
#' @export
purge.glm <- function(model) {
  model.copy <- model
  model.copy$y <- c()
  model.copy$model <- c()

  model.copy$residuals <- c()
  model.copy$fitted.values <- c()
  model.copy$effects <- c()
  model.copy$qr$qr <- c()
  model.copy$linear.predictors <- c()
  model.copy$weights <- c()
  model.copy$prior.weights <- c()
  model.copy$data <- c()

  model.copy$family$variance <- c()
  model.copy$family$dev.resids <- c()
  model.copy$family$aic <- c()
  model.copy$family$validmu <- c()
  model.copy$family$simulate <- c()
  attr(model.copy$terms,".Environment") <- c()
  attr(model.copy$formula,".Environment") <- c()

  return(model.copy)
}

#' @describeIn purge Purges an lm model
#' @export
purge.lm <- purge.glm

#' @describeIn purge Purges a merMod, linear mixed-effects model
#' @export
purge.merMod <- function(model) {
  model.copy <- model
  model.copy@resp <- lme4::lmerResp$new(y=numeric())
  model.copy@frame <- model.copy@frame[0, ]
  attr(attr(model.copy@frame, "terms"), ".Environment") <- c()
  attr(attr(model.copy@frame, "formula"), ".Environment") <- c()
  model.copy@call <- call('round', 0.0)
  model.copy@flist <- list()
  for (named.key in names(model@flist)) {
    model.copy@flist[[named.key]] <- unique(model@flist[[named.key]][0])
  }
  attr(model.copy@flist, 'assign') <- attr(model@flist, 'assign')

  # set internal lme4::merPredD object
  ss <- 10
  X <- t(matrix(lme4::getME(model, 'X')[1:ss,], ncol=ss))
  attr(X, 'assign') <- attr(lme4::getME(model, 'X'), 'assign')
  X.dimnames <- list()
  X.dimnames[[1]] <- attr(model@pp$X, 'dimnames')[[1]][1:ss]
  X.dimnames[[2]] <- attr(model@pp$X, 'dimnames')[[2]]
  attr(X, 'dimnames') <- X.dimnames

  Zt <- as(matrix(lme4::getME(model, 'Zt')[,1:ss], ncol=ss), 'dgCMatrix')
  Zt.Dimnames <- vector('list', 2)
  Zt.Dimnames[[1]] <-  attr(model@pp$Zt, 'Dimnames')[[1]]
  # don't set, truncates list:
  # Zt.Dimnames[[2]] <-  attr(model@pp$Zt, 'Dimnames')[[2]]
   attr(Zt, 'Dimnames') <- Zt.Dimnames

  Lambdat <- lme4::getME(model, 'Lambdat')
  Lind <- lme4::getME(model, 'Lind')
  theta <- lme4::getME(model, 'theta')
  n <- ss

  model.copy@pp <- lme4::merPredD(X, Zt, Lambdat, Lind, theta, n)
  model.copy@pp$setDelb(lme4::getME(model, 'beta'))
  model.copy@pp$setDelu(lme4::getME(model, 'u'))

  return(model.copy)
}

#' @describeIn purge Purges a glmerMod, generalized linear mixed-effects model
#' @export
purge.glmerMod <- function(model) {
  model.copy <- model
  model.copy@resp <- lme4::glmResp$new(family=model@resp$family,
                                       y=numeric())
  model.copy@frame <- model.copy@frame[0, ]
  attr(attr(model.copy@frame, "terms"), ".Environment") <- c()
  attr(attr(model.copy@frame, "formula"), ".Environment") <- c()
  model.copy@call <- call('round', 0.0)
  model.copy@flist <- list()
  for (named.key in names(model@flist)) {
    model.copy@flist[[named.key]] <- unique(model@flist[[named.key]][0])
  }
  attr(model.copy@flist, 'assign') <- attr(model@flist, 'assign')

  # set internal lme4::merPredD object copied from a tiny merMod model
  ss <- 10
  X <- t(matrix(lme4::getME(model, 'X')[1:ss,], ncol=ss))
  attr(X, 'assign') <- attr(lme4::getME(model, 'X'), 'assign')
  X.dimnames <- list()
  X.dimnames[[1]] <- attr(model@pp$X, 'dimnames')[[1]][1:ss]
  X.dimnames[[2]] <- attr(model@pp$X, 'dimnames')[[2]]
  attr(X, 'dimnames') <- X.dimnames

  Zt <- as(matrix(lme4::getME(model, 'Zt')[,1:ss], ncol=ss), 'dgCMatrix')
  Zt.Dimnames <- vector('list', 2)
  Zt.Dimnames[[1]] <-  attr(model@pp$Zt, 'Dimnames')[[1]]
  # don't set, truncates list:
  # Zt.Dimnames[[2]] <-  attr(model@pp$Zt, 'Dimnames')[[2]]
  attr(Zt, 'Dimnames') <- Zt.Dimnames

  Lambdat <- lme4::getME(model, 'Lambdat')
  Lind <- lme4::getME(model, 'Lind')
  theta <- lme4::getME(model, 'theta')
  n <- ss

  model.copy@pp <- lme4::merPredD(X, Zt, Lambdat, Lind, theta, n)
  model.copy@pp$setDelb(lme4::getME(model, 'beta'))
  model.copy@pp$setDelu(lme4::getME(model, 'u'))

  return(model.copy)
}

#' @describeIn purge Purges an rpart model
#' @export
purge.rpart <- function(model) {
  model.copy <- model
  model.copy$where <- c()
  model.copy$y <- c()
  attr(model.copy$terms,".Environment") <- c()
  return(model.copy)
}

#' @describeIn purge Purges a random forest model
#' @export
purge.randomForest <- function(model) {
  model.copy <- model
  model.copy$y <- c()
  model.copy$predicted <- c()
  model.copy$oob.times <- c()
  attr(model.copy$terms,".Environment") <- c()
  return(model.copy)
}

#' @describeIn purge Purges a ranger model for classification, regression, or survival
#' @export
purge.ranger <- function(model) {
  model.copy <- model
  if (model.copy$treetype %in% c('Classification', 'Regression')) {
    model.copy$predictions <- c()
  }
  if (model.copy$treetype == 'Survival') {
    model.copy$chf <- c()
    model.copy$survival <- c()
  }
  return(model.copy)
}

#' @describeIn purge Purges a coxph model
#' @export
purge.coxph <- function(model) {
  model.copy <- model
  model.copy$y <- c()
  model.copy$residuals <- c()
  model.copy$linear.predictors <- c()
  attr(model.copy$terms,".Environment") <- c()
  attr(model.copy$formula,".Environment") <- c()
  return(model.copy)
}
