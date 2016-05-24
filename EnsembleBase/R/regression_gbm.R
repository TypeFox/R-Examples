GBM.Regression.Config <- setClass("GBM.Regression.Config", slots = c(n.trees="numeric", interaction.depth="numeric", shrinkage="numeric", bag.fraction="numeric")
  , validity = function(object) {
    if (object@n.trees==round(object@n.trees) && object@n.trees>0 && 
          object@interaction.depth==round(object@interaction.depth) && object@interaction.depth>0 && 
          object@shrinkage>0 && object@bag.fraction>=0 && object@bag.fraction<=1.0) TRUE
    else "invalid parameters"
  }
  , contains = "Regression.Config"
)

GBM.Regression.FitObj <- setClass("GBM.Regression.FitObj", contains = "Regression.FitObj")

make.configs.gbm.regression <- function(df=expand.grid(n.trees=c(1000,2000),interaction.depth=c(3,4),shrinkage=c(0.001,0.01,0.1,0.5), bag.fraction=0.5)) {
  ret <- lapply(1:nrow(df), function(i) {
    GBM.Regression.Config(n.trees=df$n.trees[i], interaction.depth=df$interaction.depth[i], shrinkage=df$shrinkage[i], bag.fraction=df$bag.fraction[i])
  })
}

setMethod("BaseLearner.Fit", "GBM.Regression.Config",
  function(object, formula, data, tmpfile=NULL, print.level=1) {
    y <- data[,all.vars(formula)[1]]
    est <- gbm(formula, distribution="gaussian", data=data, n.trees=object@n.trees, interaction.depth=object@interaction.depth
               , bag.fraction=object@bag.fraction, shrinkage=object@shrinkage, verbose=print.level>=1)
    pred <- predict(est, newdata=data, n.trees=object@n.trees)
    if (!is.null(tmpfile)) {
      save(est, file=tmpfile, compress=FALSE)
      rm(est); gc()
    }
    ret <- GBM.Regression.FitObj(config=object
      , est=if (is.null(tmpfile)) est else tmpfile
      , pred=pred)
    return (ret)
  }
)

predict.GBM.Regression.FitObj <- function(object, newdata=NULL, ...) {
  if (is.null(newdata)) return (object@pred)
  if (is.character(object@est)) object@est <- load.object(object@est)
  newpred <- predict(object@est, newdata=newdata, n.trees=object@config@n.trees)
  #rm(object); gc()
  return (newpred)
}

