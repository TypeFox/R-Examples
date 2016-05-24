elmtrain.formula <-
function(formula,data,nhid,actfun,...) {
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  model <- elmtrain.default(x=x,y=y,nhid=nhid,actfun=actfun, ...)
  model$call <- match.call()
  model$formula <- formula
  model
}
