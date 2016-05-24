# taken from car::utility-functions.R; should be imported into heplots
# but not sure how to do this
# [perhaps need importFrom(car, car:::df.terms, car:::df.terms.default, car:::is.aliased) in NAMESPACE?]

df.terms <- function(model, term, ...){
	UseMethod("df.terms")
}

df.terms.default <- function(model, term, ...){
	if (is.aliased(model)) stop("Model has aliased term(s); df ambiguous.")
	if (!missing(term) && 1 == length(term)){
		assign <- attr(model.matrix(model), "assign")
		which.term <- which(term == labels(terms(model)))
		if (0 == length(which.term)) stop(paste(term, "is not in the model."))
		sum(assign == which.term)
	}
	else {
		terms <- if (missing(term)) labels(terms(model)) else term
		result <- numeric(0)
		for (term in terms) result <- c(result, Recall(model, term))
		names(result) <- terms
		result
	}
}

is.aliased <- function(model){
	!is.null(alias(model)$Complete)
}
