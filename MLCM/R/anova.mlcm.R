anova.mlcm <- function(object, ..., dispersion = NULL, test = NULL){
	dotargs <- list(...)
	m1 <- object$obj
	m2 <- if (length(dotargs) > 0) lapply(dotargs, "[[", "obj")
	if (!(is.null(m1) || 
			((length(m2) > 0) && any(sapply(m2, is.null))))){
	if (length(m2) > 0) return(anova(structure(c(list(m1), m2), class = "glmlist"), dispersion = dispersion, test = test)) else
		return(anova(m1, dispersion = dispersion, test = test))
	} else {
				
		lik1 <- logLik(object)
		lik2 <- if (length(m2) > 0) lapply(dotargs, logLik)
		lik <- c(list(lik1), lik2) 
		ddf <- -diff(sapply(lik, attr, "df"))
		dlik <- -2 * diff(unlist(lik))
		mods <- sapply(c(list(object), dotargs), 
			function(x) deparse(formula(x)))
		cat("Analysis of Deviance Test\n\n")
		for(m in seq_along(mods)) 
			cat("Model", m, ":  ", unlist(mods[m]), "\n")
		p <- pchisq(abs(dlik), abs(ddf), lower.tail = FALSE)
		dd <- data.frame(Df = ddf, Deviance = dlik, p = p)
		print(dd)
	}
}

formula.mlcm <- function(x, ...){
	if (x$method == "glm") formula(x$obj) else x$formula
	}