
fitted.dyn <- function(object, ...) {
   formula <- formula(object)
   class(object) <- setdiff(class(object), "dyn")
   series <- attr(model.frame(object), "series")
   .Class <- class(series[[1]])
   NextMethod("fitted")
}

fitted.zooreg <- fitted.ts <- fitted.irts <- fitted.its <-
fitted.zoo <- function(object, ...) {
   series <- attr(model.frame(object), "series")
   idx <- as.numeric(dimnames(as.array(resid(object)))[[1]])
   tt <- time(series[[1]])[idx]
   rval <- zoo(fitted(object), tt, frequency = frequency(series[[1]]))
   as.x <- paste("as", attr(model.frame(object), ".Class"), sep = ".")
   # workaround for bug in its
   # if (length(as.x) == 1 && as.x == "as.its") as.x <- "as.its.zoo" 
   i <- which(sapply(as.x, exists))[1]
   match.fun(as.x[i])(rval)
}

residuals.dyn <- function(object, ...) {
   formula <- formula(object)
   class(object) <- setdiff(class(object), "dyn")
   .Class <- class(attr(model.frame(object), "series")[[1]])
   NextMethod("residuals")
}

residuals.zooreg <- residuals.ts <- residuals.irts <- residuals.its <-
residuals.zoo <- function(object, ...) {
   series <- attr(model.frame(object), "series")
   idx <- as.numeric(dimnames(as.array(resid(object)))[[1]])
   tt <- time(series[[1]])[idx]
   rval <- zoo(residuals(object), tt, frequency = frequency(series[[1]]))
   as.x <- paste("as", attr(model.frame(object), ".Class"), sep = ".")
   # workaround for bug in its
   # if (length(as.x) == 1 && as.x == "as.its") as.x <- "as.its.zoo" 
   i <- which(sapply(as.x, exists))[1]
   match.fun(as.x[i])(rval)
}

predict.dyn <- function(object, newdata, ...) {
   if (missing(newdata)) return(fitted(object))
   class(object) <- setdiff(class(object), "dyn")
   series <- attr(model.frame(object), "series")
   .Class <- class(series[[1]])
   NextMethod("predict")
}

predict.zooreg <- predict.ts <- predict.irts <- predict.its <-
predict.zoo <- function(object, newdata, ...) {
   if (!is.list(newdata)) newdata <- as.list(newdata)
   formula <- formula(object)
   env <- environment(formula)
   # term list except for response
   tl <- as.list(attr(terms(object), "variables")[-1])[-1]
   tt <- eval(as.call(append(merge.zoo, tl)), newdata, env)
   rval <- zoo(predict(object, newdata, ...), time(tt), frequency = frequency(tt))
   as.x <- paste("as", attr(model.frame(object), ".Class"), sep = ".")
   # workaround for bug in its
   # if (length(as.x) == 1 && as.x == "as.its") as.x <- "as.its.zoo" 
   i <- which(sapply(as.x, exists))[1]
   match.fun(as.x[i])(rval)
}
