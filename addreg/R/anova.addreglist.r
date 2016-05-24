anova.addreglist <- function(object, ..., test = NULL)
{
	doscore <- !is.null(test) && test == "Rao"
	fam <- lapply(object, "[[", "family")
	fams <- sapply(fam, "[[", "family")
	links <- sapply(fam, "[[", "link")
	if(any(fams != fams[1]) || any(links != links[1]))
		stop("model(s) do not all have the same family and link function")
	responses <- as.character(lapply(object, function(x) {
        deparse(formula(x)[[2L]])
    }))
    sameresp <- responses == responses[1L]
	if(!all(sameresp)) {
		object <- object[sameresp]
		warning(gettextf("models with response %s removed because response differs from model 1", 
            sQuote(deparse(responses[!sameresp]))), domain = NA)
	}
    ns <- sapply(object, function(x) length(x$residuals))
    if (any(ns != ns[1L]))
        stop("models were not all fitted to the same size of dataset")
	boundary <- sapply(object, "[[", "boundary")
	if(any(boundary)) {
		object <- object[!boundary]
		warning(paste("model(s) ",paste(which(boundary),collapse=", ")," removed because MLE is on the boundary of the parameter space",sep=""))
	}
	nmodels <- length(object)
	if(nmodels == 1)
		stop('anova.addreg does not support anova for a single model. Fit nested models manually and input to anova.addreg')
	resdf <- as.numeric(lapply(object, function(x) x$df.residual))
	resdev <- as.numeric(lapply(object, function(x) x$deviance))
	if (doscore) {
		score <- numeric(nmodels)
		score[1] <- NA
		df <- -diff(resdf)
		for(i in seq_len(nmodels - 1)) {
			m1 <- if(df[i] > 0)
				object[[i]]
			else object[[i + 1]]
			m2 <- if(df[i] > 0)
				object[[i + 1]]
			else object[[i]]
			r <- m1$residuals
			w <- m1$weights * m1$family$mu.eta(m1$linear.predictors)^2/m1$family$variance(m1$fitted.values)
			zz <- glm.fit(x = m2$x, y = r, weights = w)
			score[i + 1] <- zz$null.deviance - zz$deviance
			if (df[i] < 0) score[i + 1] <- -score[i + 1]
		}
	}
	table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(NA, -diff(resdev)))
	variables <- lapply(object, function(x) paste(deparse(x$formula),
        collapse = "\n"))
	dimnames(table) <- list(1L:nmodels, c("Resid. Df", "Resid. Dev",
		"Df", "Deviance"))
	if (doscore) table <- cbind(table, Rao = score)
	title <- "Analysis of Deviance Table\n"
	topnote <- paste("Model ", format(1L:nmodels), ": ", variables, sep = "", collapse = "\n")
	if(!is.null(test)) {
		bigmodel <- object[[order(resdf)[1L]]]
		dispersion <- summary(bigmodel)$dispersion
		df.dispersion <- if (dispersion == 1)
			Inf
		else min(resdf)
		if (test == "F" && df.dispersion == Inf) {
			bigfam <- bigmodel$family$family
			if (bigfam == "binomial" || bigfam == "poisson")
				warning(gettextf("using F test with a '%s' family is inappropriate",
					bigfam), domain = NA, call. = FALSE)
			else warning("using F test with a fixed dispersion is inappropriate")
		}
		table <- stat.anova(table = table, test = test, scale = dispersion,
			df.scale = df.dispersion, n = length(bigmodel$residuals))
	}
	structure(table, heading = c(title, topnote), class = c("anova", "data.frame"))
}