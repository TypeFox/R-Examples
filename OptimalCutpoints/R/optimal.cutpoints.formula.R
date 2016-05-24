optimal.cutpoints.formula <-
function(X, tag.healthy, methods, data, direction = c("<", ">"), categorical.cov = NULL, pop.prev = NULL, control = control.cutpoints(), ci.fit = FALSE, conf.level = 0.95, trace = FALSE, ...) {
	if(missing(X)) {
		stop("'X' argument required.", call.=FALSE)
	}
	marker <- all.vars(X)[attr(terms(X), "response")]
	status <- attr(terms(X), "term.labels")
	if(length(marker) != 1 | length(status) != 1) {
		stop("Invalid formula. Please correct", call.=FALSE)
	}
	res <- optimal.cutpoints.default(X = marker, status = status, tag.healthy = tag.healthy, methods = methods, data = data, direction = direction, categorical.cov = categorical.cov, pop.prev = pop.prev, control = control, ci.fit = ci.fit, conf.level = conf.level, trace = trace)
	res$call <- match.call()
	res
}
