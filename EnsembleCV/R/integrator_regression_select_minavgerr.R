# extending the base classes of EnsembleBase package
Regression.Select.MinAvgErr.Config <- setClass("Regression.Select.MinAvgErr.Config", slots=c(instance.list="Instance.List"), contains = "Regression.Select.Config")
Regression.Select.MinAvgErr.FitObj <- setClass("Regression.Select.MinAvgErr.FitObj", contains = "Regression.Select.FitObj")

# extending the base method of EnsembleBase package
setMethod("Regression.Select.Fit", "Regression.Select.MinAvgErr.Config",
  function(object, X, y, print.level=1) {
    config.list <- lapply(object@instance.list@instances, function(instance) instance@config)
    config.list.unique <- unique(config.list)
    errors <- apply(X, 2, object@errfun, y)
    errors.avg <- sapply(config.list.unique, function(config) {
      index.select <- sapply(config.list, function(conf) identical(conf, config))
      mean(errors[index.select])
    })
    
    errors.minavg <- min(errors.avg)
    config.opt <- config.list.unique[[which(errors.avg==errors.minavg)[1]]]
    index.select <- sapply(config.list, function(conf) identical(conf, config.opt))
    pred <- rowMeans(X[,index.select,drop=F]) # this won't be used in the ecv code, instead base learners trained on full training set are used to generate pred
    est <- list(config.opt=config.opt, error.opt=errors.minavg, errors=errors.avg)
    ret <- Regression.Select.MinAvgErr.FitObj(config=object, est=est, pred=pred)

    return (ret)
  }
)

predict.Regression.Select.MinAvgErr.FitObj <- function(object, Xnew=NULL, config.list, ...) {
  if (is.null(Xnew)) stop("this method requires Xnew argument")
  if (ncol(Xnew)!=length(config.list)) stop("dimension mismatch between Xnew and config.list")
  
  index.opt <- which(sapply(config.list, function(config) identical(config, object@est$config.opt)))
  if (length(index.opt)>1) stop("multiple matches found in config.list")
  
  return (Xnew[,index.opt])
}

summary.Regression.Select.MinAvgErr.FitObj <- function(object, ...) {
  ret <- c(object@est, list(n=length(object@est$errors)))
  class(ret) <- "summary.Regression.Select.MinAvgErr.FitObj"
  return (ret)
}

print.summary.Regression.Select.MinAvgErr.FitObj <- function(x, ...) {
  cat("number of base learner configurations considered:", x$n, "\n")
  cat("best configuration:\n")
  print(x$config.opt)
  cat("minimum error:", x$error.opt, "\n")
}

## the following integrator classes and methods are not used in the ecv.regression code; instead, the above select classes and methods are directly
## whether using the operators directly is a better software design or committing to integrator-level function calls in ecv.regression code remains
## to be determined as they each offer their own pros and cons

# extending the abstract Regression.Integrator class from EnsembleBase package
# this is currently the standard CV-based integrator
Regression.Integrator.SelectMinAvgErr.Config <- setClass("Regression.Integrator.SelectMinAvgErr.Config"
  , slots=c(instance.list="Instance.List"), contains = "Regression.Integrator.Config")
Regression.Integrator.SelectMinAvgErr.FitObj <- setClass("Regression.Integrator.SelectMinAvgErr.FitObj", contains = "Regression.Integrator.FitObj")

setMethod("Regression.Integrator.Fit", "Regression.Integrator.SelectMinAvgErr.Config",
  function(object, X, y, print.level=1) {
    my.select.config <- Regression.Select.MinAvgErr.Config(errfun=object@errfun, instance.list=object@instance.list)
    est.select <- Regression.Select.Fit(my.select.config, X=X, y=y)
    
    est <- list(select=est.select)
    ret <- Regression.Integrator.SelectMinAvgErr.FitObj(config=object, est=est, pred=est.select@pred)
    return (ret)
  }
)

predict.Regression.Integrator.SelectMinAvgErr.FitObj <- function(object, Xnew=NULL, ...) {
  if (is.null(Xnew)) return (object@pred)
  newpred <- predict(object@est$select, Xnew)
  return (newpred)
}



