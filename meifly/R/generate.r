#' Fit all combinations of x variables ($2^p$).
#' 
#' This technique generalises \code{\link{fitbest}}.  While it is much
#' slower it will work for any type of model.
#'
#' @param y vector y values
#' @param x matrix of x values
#' @param method name of method used to fit the model, e.g
#'    \code{\link{lm}},\code{\link[MASS]{rlm}}
#' @param ... other arguments passed on to \code{method}
#' @keywords regression
#' @export
#' @examples
#' y <- swiss$Fertility
#' x <- swiss[, -1]
#' mods <- fitall(y, x, "lm")
fitall <- function(y, x, method = "lm", ...) {
  data <- cbind(y=y, x)

  combs <- do.call(expand.grid, rep(list(c(FALSE, TRUE)), ncol(x)))[-1, ]

  vars <- apply(combs, 1, function(i) names(x)[i])
  form <- paste("y ~ ", lapply(vars, paste, collapse=" + "), sep = "")
  form <- lapply(form, as.formula, env = baseenv())
  
  message("Fitting ", length(form), " models...")

  method <- as.name(method)
  fitmodel <- function(f) {
     eval(substitute(method(f, data = data, model = FALSE, ...), 
       list(f = f, method = method)))
  }
  
  models <- llply(form, fitmodel, .progress = "text")
  names(models) <- seq_along(models)
  
  new_ensemble(models, data)
}

#' Use the leaps package to generate the best subsets.
#' 
#' @param formula model formula
#' @param data data frame
#' @param nbest number of subsets of each size to record
#' @param ... other arguments passed to \code{\link[leaps]{regsubsets}}
#' @keywords regression
#' @export
#' @importFrom leaps regsubsets
#' @examples
#' y <- swiss$Fertility
#' mods <- fitbest(Fertility ~ ., swiss)
fitbest <- function(formula, data, nbest=10, ...) {
  b <- regsubsets(formula, data=data, nbest=nbest, ...)
  mat <- summary(b, matrix.logical = TRUE)$which

  intercept <- c("", "-1")[as.numeric(mat[,1])]
  vars <- apply(mat[,-1], 1, function(x) colnames(mat[, -1])[x])
  form <- paste(formula[[2]], " ~ ", lapply(vars, paste, collapse=" + "), sep = "")
  form <- lapply(form, as.formula)

  models <- lapply(form, function(f) eval(substitute(lm(f, data=data), list(f=f, data=data))))
  names(models) <- seq_along(models)
  
  new_ensemble(models, data)
}

#' General ensemble of models from models in global workspace'
#' 
#' @param modeltype model class
#' @param dataset if specified, all models must use this dataset
#' @param pattern pattern of model object names to match
#' @keywords regression
findmodels <- function(modeltype = "lm", dataset, pattern) {
  ls <- ls(".GlobalEnv", pattern=pattern)
  mods <- ls[sapply(ls, function(x) inherits(get(x), modeltype))]
  if (!missing(dataset)) {
    data.name <- function(x) as.character(x$call[["data"]])
    mods <- mods[sapply(mods, function(x) data.name == dataset)]
  }
  
  models <- lapply(mods, get)
  class(models) <- c("ensemble", class(models))
  models
}

#' Generate linear models by bootstrapping observations
#' 
#' @param formula model formula
#' @param data data set
#' @param n number of bootstrapped data sets to generate
#' @keywords regression
#' @export
lmboot <- function(formula, data, n=100) {
  models <- replicate(n, 
    lm(formula, data = data[sample(nrow(data), replace=TRUE), ]), 
    simplify = FALSE)
  names(models) <- seq_along(models)
  
  new_ensemble(models, data)
}

