# extending the abstract Regression.Select operation from EnsembleBase package

Regression.Select.MinErr.Config <- setClass("Regression.Select.MinErr.Config", contains = "Regression.Select.Config")

Regression.Select.MinErr.FitObj <- setClass("Regression.Select.MinErr.FitObj", contains = "Regression.Select.FitObj")

setMethod("Regression.Select.Fit", "Regression.Select.MinErr.Config",
  function(object, X, y, print.level=1) {
    error <- apply(X, 2, object@errfun, y)
    error.min <- min(error)
    index.min <- which(error==error.min)[1] # in case we find multiple optimal indexes
    est <- list(index.min=index.min, error.min=error.min, error=error)
    ret <- Regression.Select.MinErr.FitObj(config=object, est=est, pred=X[,index.min])
    return (ret)
  }
)

predict.Regression.Select.MinErr.FitObj <- function(object, Xnew=NULL, ...) {
  if (is.null(Xnew)) return (object@pred)
  return (Xnew[,object@est$index.min])
}

summary.Regression.Select.MinErr.FitObj <- function(object, ...) {
  ret <- c(object@est, list(n=length(object@est$error)))
  class(ret) <- "summary.Regression.Select.MinErr.FitObj"
  return (ret)
}

print.summary.Regression.Select.MinErr.FitObj <- function(x, ...) {
  cat("number of predictors:", x$n, "\n")
  cat("index of best predictor:", x$index.min, "\n")
  cat("minimum error:", x$error.min, "\n")
}

# extending the abstract Regression.Integrator class from EnsembleBase package
# this integrator is superseded by SelectMinAvgErr as the standard CV-based selection of model with minimum error (it averages errors across repeated CV partitions)
Regression.Integrator.SelectMinErr.Config <- setClass("Regression.Integrator.SelectMinErr.Config", contains = "Regression.Integrator.Config")
Regression.Integrator.SelectMinErr.FitObj <- setClass("Regression.Integrator.SelectMinErr.FitObj", contains = "Regression.Integrator.FitObj")

setMethod("Regression.Integrator.Fit", "Regression.Integrator.SelectMinErr.Config",
  function(object, X, y, print.level=1) {
    my.select.config <- Regression.Select.MinErr.Config(errfun=object@errfun)
    est.select <- Regression.Select.Fit(my.select.config, X=X, y=y)
    
    est <- list(select=est.select)
    ret <- Regression.Integrator.SelectMinErr.FitObj(config=object, est=est, pred=est.select@pred)
    return (ret)
  }
)

predict.Regression.Integrator.SelectMinErr.FitObj <- function(object, Xnew=NULL, ...) {
  if (is.null(Xnew)) return (object@pred)
  newpred <- predict(object@est$select, Xnew)
  return (newpred)
}




