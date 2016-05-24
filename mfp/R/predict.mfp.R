predict.mfp <- function(object, newdata, type = c("link", "response", "lp", "risk", "expected", "terms"), se.fit = FALSE, terms = NULL, 
	dispersion = NULL, na.action = na.pass, collapse, safe = FALSE, ...) 
{
# for glm: type in c("lp", "risk", "expected", "terms")
# for coxph: type in c("lp", "risk", "expected", "terms")
    type <- match.arg(type)
#
if(object$family$family=="Cox") {
    if (is.null(object$terms)) terms = names(object$assign)
	if (!missing(newdata)) 
		if (!missing(collapse)) 
			getFromNamespace("predict.coxph", "survival")(object$fit, newdata = newdata, type = type, se.fit = se.fit, terms = terms, collapse = collapse, safe = safe, ...) 
		else
			getFromNamespace("predict.coxph", "survival")(object$fit, newdata = newdata, type = type, se.fit = se.fit, terms = terms, safe = safe, ...) 
	else
		if (!missing(collapse)) 
			getFromNamespace("predict.coxph", "survival")(object$fit, type = type, se.fit = se.fit, terms = terms, collapse = collapse, safe = safe, ...) 
		else
			getFromNamespace("predict.coxph", "survival")(object$fit, type = type, se.fit = se.fit, terms = terms, safe = safe, ...) 
} else {
	if (!missing(newdata)) 
		predict.glm(object$fit, newdata = newdata, type = type, se.fit = se.fit, dispersion = dispersion, terms = terms, na.action = na.action, ...)
	else 
		predict.glm(object$fit, type = type, se.fit = se.fit, dispersion = dispersion, terms = terms, na.action = na.action, ...)
}
}


