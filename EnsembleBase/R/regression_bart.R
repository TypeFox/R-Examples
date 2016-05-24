BART.Regression.Config <- setClass("BART.Regression.Config"
  , slots = c(num_trees = "numeric", k = "numeric", q = "numeric", nu = "numeric")
  , validity = function(object) {
    if (object@num_trees == round(object@num_trees) && object@num_trees > 0 && 
        object@k > 0 && # seems like k can be fractional, so no validation check included for being integer 
        object@q >= 0 && object@q <= 1.0 &&
        object@nu == round(object@nu) && object@nu > 0
        ) TRUE
    else "invalid parameters"
  }
  , contains = "Regression.Config"
)

BART.Regression.FitObj <- setClass("BART.Regression.FitObj", slots = c(mm = "list"), contains = "Regression.FitObj")

make.configs.bart.regression <- function(df = rbind(cbind(expand.grid(num_trees = c(50, 100), k = c(2,3,4,5)), q = 0.9, nu = 3)
  , cbind(expand.grid(num_trees = c(50, 100), k = c(2,3,4,5)), q = 0.75, nu = 10))) {
  ret <- lapply(1:nrow(df), function(i) {
    BART.Regression.Config(num_trees = df$num_trees[i], k = df$k[i], q = df$q[i], nu = df$nu[i])
  })
}

setMethod("BaseLearner.Fit", "BART.Regression.Config",
  function(object, formula, data, tmpfile=NULL, print.level=1) {
    mf <- model.frame(formula, data, drop.unused.levels=TRUE, na.action = na.fail)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)
    y <- model.response(mf, "numeric")
    est <- bartMachine(X = as.data.frame(X), y = y, num_trees = object@num_trees, k = object@k, q = object@q, nu = object@nu
                       , verbose=print.level>=1, mem_cache_for_speed = FALSE)
    pred <- as.numeric(predict(est, new_data = as.data.frame(X)))
    gc()
    
    if (!is.null(tmpfile)) {
      save(est, file=tmpfile, compress=FALSE)
      rm(est); gc()
    }
    ret <- BART.Regression.FitObj(config = object
      , est = if (is.null(tmpfile)) est else tmpfile
      , pred = pred
      , mm = list(contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf), terms = mt, colnamesX = colnames(X))
    )
    return (ret)
  }
)

predict.BART.Regression.FitObj <- function(object, newdata=NULL, ...) {
  if (is.null(newdata)) return (object@pred)
  if (is.character(object@est)) object@est <- load.object(object@est)
  
  tt <- object@mm$terms
  Terms <- delete.response(tt)
  
  newdata <- droplevels(newdata)
  mf <- model.frame(Terms, newdata, xlev = object@mm$xlevels)
  X <- model.matrix(Terms, mf, contrasts.arg = object@mm$contrasts)
  
  newpred <- as.numeric(predict(object@est, new_data = as.data.frame(X)))
  #rm(object); gc()
  return (newpred)
}

