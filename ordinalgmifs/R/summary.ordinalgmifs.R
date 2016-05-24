summary.ordinalgmifs <-function (object, model.select = "AIC", ...)
{
	if (is.null(object$x)) {
		coef <- coef(object)
		coef
	} else {
        if (model.select == "AIC") {
            model.select = object$model.select
        }
        else if (model.select == "BIC") {
            model.select = which.min(object$BIC)
        }
        coef <- coef(object, model.select)
        if (object$probability.model == "Cumulative" | object$probability.model ==
        "ForwardCR" | object$probability.model == "BackwardCR") {
            cat(object$probability.model, "model using a ", object$link,
            " link \n")
        }
        else if (object$probability.model == "AdjCategory") {
            cat(object$probability.model, "model using a loge link \n")
        }
        else if (object$probability.model == "Stereotype") {
            cat(object$probability.model, "model using a logit link \n")
        }
        if (is.null(object$x)) {
            cat("logLik     = ", object$logLik, "\n")
        }
        if (!is.null(object$x)) {
            cat("at step    = ", model.select, "\n")
            cat("logLik     = ", object$logLik[model.select], "\n")
            cat("AIC        = ", object$AIC[model.select], "\n")
            cat("BIC        = ", object$BIC[model.select], "\n")
        }
        cat("\n")
        coef
    }
}
