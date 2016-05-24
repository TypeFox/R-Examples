survReg <- new("StatModel",
     capabilities = new("StatModelCapabilities"),
     name = "survival regression",
     dpp = ModelEnvFormula,
     fit = function(object, weights = NULL, ...){
         
         mydata <- cbind(object@get("response"), object@get("input"))
         names(mydata)[[1]] <- "y"
         if (!is.null(weights)) {
             mydata <- mydata[weights > 0, ]
             weights <- weights[weights > 0]
         }
         RET <- survreg(y ~ ., data = mydata, weights = weights, ...)
         RET$addargs <- list(...)
         RET$ModelEnv <- object
	 RET$weights <- weights
         class(RET) <- c("survReg", "survreg")
         RET
     }
)

fitted.survReg <- function(object, ...) predict(object)

weights.survReg <- function(object, ...) {
    if(is.null(object$weights)) rep(1, NROW(residuals(object))) else object$weights
}

print.survReg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    dist <- x$dist
    substr(dist, 1, 1) <- toupper(substr(dist, 1, 1))
    cat(paste(dist, "survival regression",
        paste("(scale = ", paste(format(x$scale, digits = digits), sep = ", "), ")", sep = ""),
        "with coefficients:\n"))
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    invisible(x)
}

logLik.survReg <- function(object, ...) {
  structure(object$loglik[2], df = NCOL(object$var), class = "logLik")
}

model.matrix.survReg <- function(object, data, ...) {
  if(missing(data)) return(model.matrix(object, model.frame(object), ...))
  NextMethod()
}
