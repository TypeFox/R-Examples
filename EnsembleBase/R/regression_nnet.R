NNET.Regression.Config <- setClass("NNET.Regression.Config", slots = c(decay="numeric", size="numeric", maxit="numeric")
  , validity = function(object) {
    if (object@decay>=0 && round(object@size)==object@size && object@size>0 &&
      object@maxit>0 && round(object@maxit)==object@maxit) TRUE
    else "invalid parameters"
  }
  , contains = "Regression.Config"
)

NNET.Regression.FitObj <- setClass("NNET.Regression.FitObj"
  , slots = c(y.range="numeric", y.min="numeric"), contains = "Regression.FitObj")

make.configs.nnet.regression <- function(df=expand.grid(decay=c(1e-4,1e-2,1,100), size=c(5,10,20,40), maxit=2000)) {
  ret <- lapply(1:nrow(df), function(i) {
    NNET.Regression.Config(decay=df$decay[i], size=df$size[i], maxit=df$maxit[i])
  })
}

setMethod("BaseLearner.Fit", "NNET.Regression.Config",
  function(object, formula, data, tmpfile=NULL, print.level=1) {
    respVar <- all.vars(formula)[1]
    y.original <- data[,respVar]
    y.min <- min(data[,respVar])
    y.range <- diff(range(data[,respVar]))
    if (y.range>0) {
      data[,respVar] <- (data[,respVar]-y.min)/y.range
    } else { # shouldn't we do something more drastic if y.range=0?!
      data[,respVar] <- (data[,respVar]-y.min)
    }
    est <- nnet(formula, data, size=object@size, decay=object@decay, maxit=object@maxit, trace=print.level>=1)
    pred <- predict(est)
    if (y.range>0) pred <- as.vector(pred*y.range+y.min)
    else pred <- as.vector(pred+y.min)
    if (!is.null(tmpfile)) {
      save(est, file=tmpfile, compress=FALSE)
      rm(est); gc()
    }
    ret <- NNET.Regression.FitObj(config=object
      , est=if (is.null(tmpfile)) est else tmpfile
      , pred=pred, y.range=y.range, y.min=y.min)
    return (ret)
  }
)

predict.NNET.Regression.FitObj <- function(object, newdata=NULL, ...) {
  if (is.null(newdata)) return (object@pred)
  if (is.character(object@est)) object@est <- load.object(object@est)
  newpred <- as.vector(predict(object@est, newdata=newdata))
  if (object@y.range>0) newpred <- newpred*object@y.range+object@y.min
  else newpred <- newpred+object@y.min
  #rm(object); gc()
  return (newpred)
}






