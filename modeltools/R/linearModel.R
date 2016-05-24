
### an example for an unfitted statistical model: linear model

lmfit <- function(object, weights = NULL, ...){

    ### extract design and response matrix from the `ModelEnv' object
    ### and call the usual fit methods
    if (is.null(weights)) {
        z <- lm.fit(object@get("designMatrix"),
                    object@get("responseMatrix"),
                    ...)
    } else {
        z <- lm.wfit(object@get("designMatrix"),
                     object@get("responseMatrix"), weights, ...)
    }

    ### returns a model inheriting from `mlm' or / and `lm'
    class(z) <- c("linearModel", if (is.matrix(z$fitted)) "mlm", "lm")
    z$offset <- 0
    z$contrasts <- attr(object@get("designMatrix"), "contrasts")
    z$xlevels <- attr(object@get("designMatrix"), "xlevels")
    z$terms <- attr(object@get("input"), "terms")
    
    ### predict.lm will fails since we cannot provide 
    ### correct $call and $terms elements.
    z$predict_response <- function(newdata = NULL) {
        if (!is.null(newdata)) {
            penv <- new.env()
            object@set("input", data = newdata, env = penv)
            dm <- get("designMatrix", envir = penv, inherits = FALSE)
        } else {
            dm <- object@get("designMatrix")
        }
        pr <- dm %*% coef(z)
        if (ncol(pr) == 1) pr <- drop(pr)
        return(pr)
    }
    z$addargs <- list(...)
    z$ModelEnv <- object
    z$statmodel <- linearModel
    z
}


### an object of class `StatModel' representing unfitted linear models
linearModel <- new("StatModel",
    capabilities = new("StatModelCapabilities"),
    name = "linear regression model",
    dpp = ModelEnvFormula,
    fit = lmfit,
    predict = function(object, newdata = NULL, ...)
        #### simply call the predict_response element
        object$predict_response(newdata = newdata)
)

### we would like to advocate `Predict', but anyway
predict.linearModel <- function(object, newdata = NULL, ...)
    linearModel@predict(object, newdata = newdata)

fitted.linearModel <- function(object, ...)
    object$predict_response()

weights.linearModel <- function(object, ...) {
    if(is.null(object$weights)) rep(1, NROW(object$residuals)) else object$weights
}

print.linearModel <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("Linear model with coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    invisible(x)
}

model.matrix.linearModel <- function(object, ...)
  object$ModelEnv@get("designMatrix")
