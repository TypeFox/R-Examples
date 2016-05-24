predict.logbin.smooth <- function(object, newdata = NULL, type = c("link", "response",
	"terms"), terms = NULL, na.action = na.pass, ...)
{
	type <- match.arg(type)
	na.act <- object$na.action
	object$na.action <- NULL
	object$terms <- object$terms.full
	
	if (missing(newdata))
		pred <- switch(type, link = object$linear.predictors,
					response = object$fitted.values,
					terms = predict.logbin(object, type = "terms", terms = terms))
	else {
		gp <- interpret.logbin.smooth(object$full.formula)
		smoothterms <- names(gp$smooth.spec)
		knots <- object$knots
		for (smth in smoothterms) {
			smthlabel <- gp$smooth.spec[[smth]]$termlabel
			smthtype <- class(gp$smooth.spec[[smth]])
			x <- newdata[[smth]]
			if (class(x) != "numeric")
				stop(gettextf("error in %s: variable %s in newdata must be numeric",
						smthlabel, smth), domain = NA)
			if (smthtype == "Iso.smooth") {
				if (any(x < knots[[smth]][1] | x > knots[[smth]][length(knots[[smth]])]))
					stop(gettextf("error in %s: the newdata must be in the range %g to %g",
						smthlabel, knots[[smth]][1], knots[[smth]][length(knots[[smth]])]),
						domain = NA)
				x.new <- matrix(0, nrow = nrow(newdata), ncol = length(knots[[smth]]) - 1)
				for(i in 2:length(knots[[smth]]))
					x.new[, (i-1)] <- as.numeric(x >= knots[[smth]][i])
				colnames(x.new) <- paste(smthlabel, 2L:length(knots[[smth]]), sep = "")
			} else if (smthtype == "B.smooth") {
				B <- matrix(NA, nrow = nrow(newdata), ncol = length(knots[[smth]]-3))
				B <- splines::splineDesign(knots[[smth]], x, ord = 3, outer.ok = FALSE)
				colnames(B) <- paste(smthlabel, seq_len(ncol(B)), sep = "")
				x.new <- B[, -1, drop = FALSE]
			}
			newdata <- cbind(newdata, x.new)
		}
		pred <- predict.logbin(object, newdata, type, terms, na.action)
	}
	if (missing(newdata) && !is.null(na.act))
		pred <- napredict(na.act, pred)
	
	pred
}