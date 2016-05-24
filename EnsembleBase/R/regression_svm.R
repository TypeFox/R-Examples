SVM.Regression.Config <- setClass("SVM.Regression.Config", slots = c(cost="numeric", epsilon="numeric", kernel="character")
  , validity = function(object) {
   if (object@cost>=0 && object@epsilon>0 && object@kernel %in% c("linear","polynomial","radial","sigmoid")) TRUE
   else "invalid parameters"
  }
  , contains = "Regression.Config"
)

SVM.Regression.FitObj <- setClass("SVM.Regression.FitObj", contains = "Regression.FitObj")

make.configs.svm.regression <- function(df=expand.grid(cost=c(0.1,0.5,1.0,5.0,10,50,75,100), epsilon=c(0.1,0.25), kernel="radial")) {
  ret <- lapply(1:nrow(df), function(i) {
    SVM.Regression.Config(cost=df$cost[i], epsilon=df$epsilon[i], kernel=as.character(df$kernel[i]))
  })
}

setMethod("BaseLearner.Fit", "SVM.Regression.Config",
  function(object, formula, data, tmpfile=NULL, print.level=1) {
    respVar <- all.vars(formula)[1]
    y <- data[,respVar]
    est <- svm(formula, data, kernel=object@kernel, cost=object@cost, epsilon=object@epsilon)
    pred <- as.vector(predict(est))
    if (!is.null(tmpfile)) {
      save(est, file=tmpfile, compress=FALSE)
      rm(est); gc()
    }
    ret <- SVM.Regression.FitObj(config=object
      , est=if (is.null(tmpfile)) est else tmpfile
      , pred=pred)
    return (ret)
  }
)

predict.SVM.Regression.FitObj <- function(object, newdata=NULL, ...) {
  if (is.null(newdata)) return (object@pred)
  if (is.character(object@est)) object@est <- load.object(object@est)
  newpred <- as.vector(predict(object@est, newdata=newdata, na.action=na.pass))
  #rm(object); gc()
  return (newpred)
}



