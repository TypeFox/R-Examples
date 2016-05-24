# defining base classes and generic methods to support sweep operation; TODO: move these base classes to EnsembleBase?

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
  pred <- array(NA, dim=dim(X))
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

### PCR sweep methods and classes ###
Regression.Sweep.PCR.Config <- setClass("Regression.Sweep.PCR.Config", contains = "Regression.Sweep.Config")
Regression.Sweep.PCR.FitObj <- setClass("Regression.Sweep.PCR.FitObj", contains = "Regression.Sweep.FitObj")
# training
setMethod("Regression.Sweep.Fit", "Regression.Sweep.PCR.Config",
  function(object, X, y) {
    maxpc <- min(dim(X))
    pcr <- prcomp(X)
    Xpc <- cbind(1,predict(pcr, X))
    XtX <- t(Xpc)%*%Xpc
    Xty <- t(Xpc)%*%y
    predmat <- array(NA, dim=c(nrow(X),maxpc))
    beta.list <- list()
    for (i in 1:maxpc) {
      # add a try wrapper, trap error and exit loop; this allows us to go to higher pc's and not arbitrarily choose a cutoff
      betatmp <- try(solve(XtX[1:(i+1),1:(i+1)], Xty[1:(i+1)]),T) # using normal equations, faster but less stable
      if (inherits(betatmp, "try-error")) {
        #cat("singularity encountered at i=", i, "\n")
        i <- i-1
        break
      }
      predmat[,i] <- Xpc[,1:(i+1)]%*%betatmp
      beta.list[[i]] <- betatmp
    }
    object@n <- i
    est <- list(pcr=pcr, beta.list=beta.list, dimX=dim(X))
    ret <- Regression.Sweep.PCR.FitObj(config=object, est=est, pred=predmat)
    return (ret)
  }
)
# prediction
predict.Regression.Sweep.PCR.FitObj <- function(object, Xnew=NULL, ...) {
  if (is.null(Xnew)) return (object@pred)
  maxpc <- object@config@n
  Xpc.new <- cbind(1,predict(object@est$pcr, Xnew)[,1:maxpc,drop=F])
  newpred <- array(NA, dim=c(nrow(Xnew),maxpc))
  for (i in 1:maxpc) {
    newpred[,i] <- Xpc.new[,1:(i+1)]%*%object@est$beta.list[[i]]
  }
  return (newpred)
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

Regression.Integrator.PCR.SelMin.Config <- setClass("Regression.Integrator.PCR.SelMin.Config", slots = c(partition="integer"), contains = "Regression.Integrator.Config")
Regression.Integrator.PCR.SelMin.FitObj <- setClass("Regression.Integrator.PCR.SelMin.FitObj", contains = "Regression.Integrator.FitObj")
setMethod("Regression.Integrator.Fit", "Regression.Integrator.PCR.SelMin.Config",
  function(object, X, y, print.level=1) {
    # 1) cv sweep - pcr method
    my.pcr.config <- Regression.Sweep.PCR.Config(n=NULL)
    est.pcr <- Regression.Sweep.CV.Fit(my.pcr.config, X=X, y=y, partition=object@partition, print.level=print.level)
    # 2) select - min method
    my.select.config <- Regression.Select.MinErr.Config(errfun=object@errfun)
    est.select <- Regression.Select.Fit(my.select.config, X=est.pcr@pred, y=y)
    
    est <- list(pcr=est.pcr, select=est.select)
    ret <- Regression.Integrator.PCR.SelMin.FitObj(config=object, est=est, pred=est.select@pred)
    return (ret)
  }
)
predict.Regression.Integrator.PCR.SelMin.FitObj <- function(object, Xnew=NULL, ...) {
  if (is.null(Xnew)) return (object@pred)
  newpred.pcr <- predict(object@est$pcr, Xnew)
  newpred <- predict(object@est$select, newpred.pcr)
  attr(newpred, "newpred.pcr") <- newpred.pcr
  return (newpred)
}




