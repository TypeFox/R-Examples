KNN.Regression.Config <- setClass("KNN.Regression.Config", slots = c(kernel="character", k="numeric")
  , validity = function(object) {
    if (object@kernel %in% c("rectangular","epanechnikov","triweight","gaussian","tiangular","biweight","cos","inv","rank","optimal") && 
          object@k>0 && object@k==round(object@k)) TRUE
    else "invalid parameters"
  }
  , contains = "Regression.Config"
)

KNN.Regression.FitObj <- setClass("KNN.Regression.FitObj"
  , slots = c(formula="formula", data="data.frame"), contains = "Regression.FitObj")

make.configs.knn.regression <- function(df=expand.grid(kernel=c("rectangular","epanechnikov","triweight","gaussian"),k=c(5,10,20,40))) {
  ret <- lapply(1:nrow(df), function(i) {
    KNN.Regression.Config(kernel=as.character(df$kernel[i]), k=df$k[i])
  })
}

setMethod("BaseLearner.Fit", "KNN.Regression.Config",
  function(object, formula, data, tmpfile=NULL, print.level=1) {
    y <- data[,all.vars(formula)[1]]
    est <- kknn(formula, data, data, k=object@k, kernel=object@kernel)
    pred <- est$fitted.values
    if (!is.null(tmpfile)) {
      save(est, file=tmpfile, compress=FALSE)
      rm(est); gc()
    }
    ret <- KNN.Regression.FitObj(config=object
      , est=if (is.null(tmpfile)) est else tmpfile
      , pred=pred, formula=formula, data=data)
    return (ret)
  }
)

predict.KNN.Regression.FitObj <- function(object, newdata=NULL, ...) {
  if (is.null(newdata)) return (object@pred)
  if (is.character(object@est)) object@est <- load.object(object@est)
  newreg <- kknn(object@formula, object@data, newdata, k=object@config@k, kernel=object@config@kernel)
  #rm(object); gc()
  return (newreg$fitted.values)
}


