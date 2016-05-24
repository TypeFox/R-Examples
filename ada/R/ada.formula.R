"ada.formula" <-
function(formula, data,...,subset,na.action=na.rpart){
 ## m = match.call(expand.dots = FALSE)
  ##m[[1]] = as.name("model.frame")
  ##m$...=NULL
  ##m =eval(m,parent.frame())

  m <- match.call(expand.dots = FALSE)
  m$model <- m$method <- m$control <- NULL
  m$x <- m$y <- m$parms <- m$... <- NULL
  m$cost <- NULL
  m$na.action <- na.action
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  
  Terms = attr(m, "terms")
  y = as.vector(model.extract(m,"response"))
  preds<-attr(attributes(m)$terms,"term.labels")
  x<-as.data.frame(m[,!is.na(match(names(m),preds))])
  res = ada.default(x,y,...,na.action=na.action)
  res$terms = Terms
  cl = match.call()
  cl[[1]] = as.name("ada")
  res$call = cl
  res
}



