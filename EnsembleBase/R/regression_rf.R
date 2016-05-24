RF.Regression.Config <- setClass("RF.Regression.Config", slots = c(ntree="numeric", nodesize="numeric", mtry.mult="numeric")
  , validity = function(object) {
    if (object@ntree==round(object@ntree) && object@ntree>1 && 
          object@nodesize==round(object@nodesize) && object@nodesize>0 &&
          object@mtry.mult>0) TRUE
    else "invalid parameters"
  }
  , contains = "Regression.Config"
)

RF.Regression.FitObj <- setClass("RF.Regression.FitObj", contains = "Regression.FitObj")

make.configs.rf.regression <- function(df=expand.grid(ntree=c(100,500),mtry.mult=c(1,2),nodesize=c(2,5,25,100))) {
  ret <- lapply(1:nrow(df), function(i) {
    RF.Regression.Config(ntree=df$ntree[i], nodesize=df$nodesize[i], mtry.mult=df$mtry.mult[i])
  })
}

setMethod("BaseLearner.Fit", "RF.Regression.Config",
  function(object, formula, data, tmpfile=NULL, print.level=1) {
    y <- data[,all.vars(formula)[1]]
    varnames <- labels(terms(formula))
    est <- randomForest(formula, data, ntree=object@ntree
      , nodesize=object@nodesize
      , mtry=max(floor(object@mtry.mult*length(varnames)/3), 1), do.trace=print.level>=1, keep.forest=T)
    pred <- as.vector(est$predicted)
    if (!is.null(tmpfile)) {
      save(est, file=tmpfile, compress=FALSE)
      rm(est); gc()
    }
    ret <- RF.Regression.FitObj(config=object
      , est=if (is.null(tmpfile)) est else tmpfile
      , pred=pred)
    return (ret)
  }
)

predict.RF.Regression.FitObj <- function(object, newdata=NULL, ...) {
  if (is.null(newdata)) return (object@pred)
  if (is.character(object@est)) object@est <- load.object(object@est)
  newpred <- as.vector(predict(object@est, newdata=newdata))
  #rm(object); gc()
  return (newpred)
}






