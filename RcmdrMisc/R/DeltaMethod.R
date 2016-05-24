# last modified 2016-04-01

DeltaMethod <- function(model, g, level=0.95){
  coefs <- coef(model)
  p <- length(coefs)
  nms <- if (names(coefs)[1] == "(Intercept)") paste0("b", 0:(p - 1)) else paste0("b", 1:p)
  res <- car::deltaMethod(model, g, parameterNames=nms)
  result <- list(test=res, coef=rbind(names(coefs), nms))
  class(result) <- "DeltaMethod"
  result
}

print.DeltaMethod <- function(x, ...){
  coef <- x$coef
  par <- data.frame(t(coef))
  colnames(par) <- c("parameter", "name")
  print(par, row.names=FALSE)
  cat("\n")
  print(x$test)
  invisible(x)
}
