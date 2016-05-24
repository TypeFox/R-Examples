# partial eta^2 measures of association for multivariate tests
# (mod of car:::print.Anova.mlm, just for testing)
# added etasq.lm 7/29/2010
# fixed buglet in etasq.lm 12/5/2010
# copied stats:::Pillai, etc. to utility.R to avoid using :::

etasq <- function(x, ...){
	UseMethod("etasq", x)
}

etasq.mlm <- function(x, ...) {
	etasq(Anova(x, ...))
}

etasq.Anova.mlm <- function(x, anova=FALSE, ...){
	test <- x$test
	if (test == "Roy") warning("eta^2 not defined for Roy's test")
	repeated <- x$repeated
	ntests <- length(x$terms)
	tests <- matrix(NA, ntests, 4)
	assoc <- rep(NA, ntests)
	if (!repeated) SSPE.qr <- qr(x$SSPE) 
	for (term in 1:ntests){
		# some of the code here adapted from stats:::summary.manova
		eigs <- Re(eigen(qr.coef(if (repeated) qr(x$SSPE[[term]]) else SSPE.qr,
								x$SSP[[term]]), symmetric = FALSE)$values)
		tests[term, 1:4] <- switch(test,
				Pillai = Pillai(eigs, x$df[term], x$error.df),
				Wilks = Wilks(eigs, x$df[term], x$error.df),
				"Hotelling-Lawley" = HL(eigs, x$df[term], x$error.df),
				Roy = Roy(eigs, x$df[term], x$error.df))
		s <- min(length(eigs), x$df[term])
		assoc[term] <- switch(test,
			Pillai = tests[term,1] / s,
			Wilks = 1 - tests[term,1] ^(1/s),
			"Hotelling-Lawley" = tests[term,1] / (tests[term,1] + s),
			Roy = tests[term,1] / (tests[term,1] + 1))
			
	}
	ok <- tests[, 2] >= 0 & tests[, 3] > 0 & tests[, 4] > 0
	ok <- !is.na(ok) & ok
	if(anova) {
  	result <- cbind(assoc, x$df, tests, pf(tests[ok, 2], tests[ok, 3], tests[ok, 4], 
  					lower.tail = FALSE))
  	rownames(result) <- x$terms
  	colnames(result) <- c("eta^2", "Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)")
  	result <- structure(as.data.frame(result), 
  			heading = paste("\nType ", x$type, if (repeated) " Repeated Measures",
  					" MANOVA Tests: ", test, " test statistic", sep=""), 
  			class = c("anova", "data.frame"))
	}
	else {
		result <- data.frame(assoc)
  	rownames(result) <- x$terms
  	colnames(result) <- "eta^2"
		
	}
	
	result      
}

etasq.lm <- function(x, anova=FALSE, partial=TRUE, ...) {
	aov <-Anova(x, ...)
	neff <- nrow(aov)
	SSH <- aov[-neff,1]
	SSE <- aov[neff,1]
	SST <- sum(SSH) + SSE
	eta2 <- if (partial) c(SSH / (SSH + SSE), NA) else c(SSH / SST, NA)
	etalab <- if (partial) "Partial eta^2" else "eta^2"
	if (anova) {
		result <- cbind(eta2, aov)
		rownames(result) <- rownames(aov)
		colnames(result) <- c(etalab, colnames(aov))
		result <- structure(as.data.frame(result), 
				heading = attr(aov, "heading"), 
				class = c("anova", "data.frame"))
	}
	else {
		result <- data.frame(eta2)
		rownames(result) <- rownames(aov)
		colnames(result) <- etalab
	}
	result      
}

