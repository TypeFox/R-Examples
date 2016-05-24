glinearModel <- new("StatModel",
    capabilities = new("StatModelCapabilities"),
    name = "generalized linear regression model",
    dpp = ModelEnvFormula,
    fit = function(object, weights = NULL, ...){

        if (is.null(weights)) {
            z <- glm.fit(x = object@get("designMatrix"),
                         y = object@get("response")[,1],
                         intercept = all(object@get("designMatrix")[,1] == 1),
                         ...)
        } else {
            z <- glm.fit(x = object@get("designMatrix"),
                         y = object@get("response")[,1], 
                         weights = weights, 
                         intercept = all(object@get("designMatrix")[,1] == 1),
                         ...)
        }
        class(z) <- c("glinearModel", "glm", "lm")
        z$offset <- 0
        z$contrasts <- attr(object@get("designMatrix"), "contrasts")

        ## terms should be there, but still need to
	## be worked around in predictions
        z$terms <- attr(object@get("input"), "terms")
        z$predict_response <- function(newdata = NULL) {
            if (!is.null(newdata)) {
                penv <- new.env()
                object@set("input", data = newdata, env = penv)
                dm <- get("designMatrix", envir = penv, inherits = FALSE)
            } else {
                dm <- object@get("designMatrix")
            }
            pr <- z$family$linkinv(drop(dm %*% z$coef))
            return(pr)
        }
        z$addargs <- list(...)
	z$ModelEnv <- object
        z
    },
    predict = function(object, newdata = NULL, ...) 
        object$predict_response(newdata = newdata)
)

predict.glinearModel <- function(object, newdata = NULL, ...)
    object$predict_response(newdata = newdata) 

fitted.glinearModel <- function(object, ...)
    object$predict_response()

print.glinearModel <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    fam <- x$family$family
    substr(fam, 1, 1) <- toupper(substr(fam, 1, 1))
    cat(paste(fam, "GLM with coefficients:\n"))
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    invisible(x)
}

model.matrix.glinearModel <- function(object, ...)
  object$ModelEnv@get("designMatrix")
