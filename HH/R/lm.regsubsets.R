"lm.regsubsets" <-
function(object, model.number, ...) {
  sum.reg <- summary(object)
  dim.sum.reg <- dim(sum.reg$outmat)
  vars <- sum.reg$outmat[model.number,]
  rhs <- paste(names(vars)[vars=="*"], collapse="+")
  rhs <- parse(text=paste("~", rhs))[[1]][[2]]
  lm.call <- object$call
  lm.call[[1]] <- as.name("lm")
  lm.call[[2]][[3]] <- rhs
  llmc <- length(lm.call)
  if (llmc > 3) for (i in 4:llmc) lm.call[[i]] <- NULL
  lm.call[[1]] <- as.name("lm")
  eval.parent(lm.call, 2)
}

