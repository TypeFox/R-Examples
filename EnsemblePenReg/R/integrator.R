# defining base classes and generic methods to support sweep operation
# TODO: reconcile with same code in EnsemblePCReg; move these base classes to EnsembleBase

# base class for sweep configurations
setClass("Regression.Sweep.Config", slots=c(n="OptionalNumeric"), contains="VIRTUAL")
# base class for output of sweep training
setClass("Regression.Sweep.FitObj", slots = c(config="Regression.Sweep.Config", est="ANY", pred="matrix"), contains = "VIRTUAL")
# generic method for training sweep operations
setGeneric("Regression.Sweep.Fit", function(object, X, y, print.level=1) standardGeneric("Regression.Sweep.Fit"))
# class for output of sweep cv training
Regression.Sweep.CV.FitObj <- setClass("Regression.Sweep.CV.FitObj"
  , slots = c(sweep.list="list", pred="matrix", partition="OptionalInteger")
)

# sweep cv training function
Regression.Sweep.CV.Fit <- function(config, X, y, partition, print.level=1) {
  nfolds <- max(partition)
  if (length(partition)!=nrow(X)) stop("length of fold parameter does not match number of rows in data")
  #pred <- array(NA, dim=dim(X))
  pred <- array(NA, dim=c(nrow(X),config@n)) # config@n plays the rols of maximum possible size, dim(X) was appropriate only for PCR
  sweep.list <- list()
  n.min <- +Inf
  for (i in 1:nfolds) {
    if (print.level>=1) cat("processing fold", i, "of", nfolds, "\n")
    index_predict <- which(partition==i)
    X_train <- X[-index_predict,]
    X_predict <- X[index_predict,]
    y_train <- y[-index_predict]
    y_predict <- y[index_predict]
    regtmp <- Regression.Sweep.Fit(object=config, X=X_train, y=y_train, print.level=print.level-1)
    pred[index_predict,1:regtmp@config@n] <- predict(regtmp, Xnew=X_predict)
    sweep.list[[i]] <- regtmp
    n.min <- min(n.min, regtmp@config@n) # TODO: loop-carried, dependency, cannot be parallelized in current form
    config <- regtmp@config # for penreg, we want to use the same sequence of lambda's for all folds, so the first fold dictates the config for the rest
  }
  # this is needed specifically for pcr sweep since we are allowing for n in each fold to be determined using the try-error method
  # alternatively, we can keep this function general, and expect the PCR sweep to define its own CV method
  for (i in 1:nfolds) {
    sweep.list[[i]]@config@n <- n.min
  }
  pred <- pred[,1:n.min,drop=F]
  
  ret <- Regression.Sweep.CV.FitObj(sweep.list=sweep.list, pred=pred, partition=partition)
  return (ret)
}
# sweep cv predict
predict.Regression.Sweep.CV.FitObj <- function(object, Xnew=NULL, ...) {
  if (is.null(Xnew)) return (object@pred)
  nfolds <- length(object@sweep.list)
  ret <- 0
  for (n in 1:nfolds) {
    ret <- ret + predict(object@sweep.list[[n]], Xnew=Xnew)
  }
  ret <- ret/nfolds
  return (ret)
}

## end of sweep operation code that must be moved to EnsembleBase

### PenReg sweep methods and classes ###
# alpha=1.0 -> lasso, alpha=0.0 -> ridge
# lambda: shrinkage parameter
Regression.Sweep.PenReg.Config <- setClass("Regression.Sweep.PenReg.Config", slots=c(alpha="numeric", lambda="numeric"), contains = "Regression.Sweep.Config")
Regression.Sweep.PenReg.FitObj <- setClass("Regression.Sweep.PenReg.FitObj", contains = "Regression.Sweep.FitObj")

# training
setMethod("Regression.Sweep.Fit", "Regression.Sweep.PenReg.Config",
  function(object, X, y) {
    if (length(object@lambda)==0) {
      # if lambda is not provided, we let glmnet automatically determine lambda
      penreg <- glmnet(x=X, y=y, family="gaussian", alpha=object@alpha, nlambda=object@n)
    } else {
      penreg <- glmnet(x=X, y=y, family="gaussian", alpha=object@alpha, lambda=object@lambda)
    }
    object@lambda <- penreg$lambda
    object@n <- length(object@lambda)
    pred <- predict(penreg, newx=X)
    ret <- Regression.Sweep.PenReg.FitObj(config=object, est=penreg, pred=pred)
  }
)
# prediction
predict.Regression.Sweep.PenReg.FitObj <- function(object, Xnew=NULL, ...) {
  if (is.null(Xnew)) return (object@pred)
  return (predict(object@est, Xnew))
}

### PCR integrator methods and classes ###

## Select.MinErr operation, TODO: this must be exported in EnsembleCV and imported in this package, or moved to EnsembleBase
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
## end of Select.MinErr operation

# can the Sweep operation be abstracted out and turned into an argument to a more generic function, thus uniting PCR and PenReg?
Regression.Integrator.PenReg.SelMin.Config <- setClass("Regression.Integrator.PenReg.SelMin.Config", slots = c(partition="integer", n="OptionalNumeric", alpha="numeric"), contains = "Regression.Integrator.Config")
Regression.Integrator.PenReg.SelMin.FitObj <- setClass("Regression.Integrator.PenReg.SelMin.FitObj", contains = "Regression.Integrator.FitObj")
setMethod("Regression.Integrator.Fit", "Regression.Integrator.PenReg.SelMin.Config",
  function(object, X, y, print.level=1) {
    # 1) cv sweep - penreg method
    my.penreg.config <- Regression.Sweep.PenReg.Config(n=object@n, alpha=object@alpha)
    est.penreg <- Regression.Sweep.CV.Fit(my.penreg.config, X=X, y=y, partition=object@partition, print.level=print.level)
    # 2) select - min method
    my.select.config <- Regression.Select.MinErr.Config(errfun=object@errfun)
    est.select <- Regression.Select.Fit(my.select.config, X=est.penreg@pred, y=y)
    
    est <- list(penreg=est.penreg, select=est.select)
    ret <- Regression.Integrator.PenReg.SelMin.FitObj(config=object, est=est, pred=est.select@pred)
    return (ret)
  }
)
predict.Regression.Integrator.PenReg.SelMin.FitObj <- function(object, Xnew=NULL, ...) {
  if (is.null(Xnew)) return (object@pred)
  newpred.penreg <- predict(object@est$penreg, Xnew)
  newpred <- predict(object@est$select, newpred.penreg)
  attr(newpred, "newpred.penreg") <- newpred.penreg
  return (newpred)
}




