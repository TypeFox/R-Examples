utils::globalVariables(c('Outcome','ID','Treatment'))

#' Provides a summary of the matched triplets including analysis of outcome
#' measure if provided.
#' 
#' If an outcome measure is provided this function will perform a Freidman
#' Rank Sum Test and repeated measures ANOVA. If either test has a statistically
#' significant difference (as determined by the value of the \code{p} parameter),
#' a Pairwise Wilcoxon Rank Sum Test will also be provided.
#' 
#' @param object result of \code{\link{trimatch}}.
#' @param outcome vector representing the outcome measure.
#' @param p threshold of the p value to perform a 
#' @param ... parameters passed to other statistical tests.
#' @param ordering specify the order for doing the paired analysis, that is
#'        analysis will be conducted as:
#'        \code{ordering[1] - ordering[2]}, \code{ordering[1] - ordering[3]},
#'        and \code{ordering[2] - ordering[3]}.
#' @seealso \code{\link{friedman.test}}, \code{\link{ezANOVA}}, 
#'        \code{\link{pairwise.wilcox.test}}
#' @return a trimatch.summary object.
#' @method summary triangle.matches
#' @export
summary.triangle.matches <- function(object, outcome, p=.05, 
									 ordering = attr(object, 'match.order'), ...) {
	results <- list()
	
	tpsa <- attr(object, 'triangle.psa')
	
	tab <- table(tpsa$treat)
	tab2 <- c(length(unique(object[,names(tab)[1]])),
		length(unique(object[,names(tab)[2]])),
		length(unique(object[,names(tab)[3]])) )
	
	results$PercentMatched <- as.numeric(tab2 / tab)
	names(results$PercentMatched) <- names(tab)
		
	if(missing(outcome)) {
		stop('Outcome measure not specified.')
	} else {
		tmatch.out <- merge(object, outcome)
		#outcomes <- grep(".out$", names(tmatch.out), perl=TRUE)
		outcomes <- sapply(ordering, function(x) { 
				which(names(tmatch.out) == paste0(x, '.out')) })
		tmatch.out$id <- 1:nrow(tmatch.out)
		out <- melt(tmatch.out[,c(outcomes, which(names(tmatch.out) == 'id'))],id.vars='id')
		names(out) <- c('ID','Treatment','Outcome')
		out$ID <- as.factor(out$ID)
		results$friedman.test <- friedman.test(Outcome ~ Treatment | ID, out, ...)
		
		#Repeated measures ANOVA
		results$rmanova <- ezANOVA(data=out, dv=Outcome, wid=ID, 
								   within=Treatment, ...)
		
		if(results$rmanova$ANOVA$p <= p || results$friedman.test$p.value <= p) {
			#Possible approach for post-hoc test
			results$pairwise.wilcox.test <- pairwise.wilcox.test(
						x=out$Outcome, g=out$Treatment, paired=TRUE, 
						p.adjust.method='bonferroni', ...)
			
			# Paired t-tests
			diffcols <- c(
				paste(names(tmatch.out)[outcomes[c(1,2)]], collapse='-'),
				paste(names(tmatch.out)[outcomes[c(1,3)]], collapse='-'),
				paste(names(tmatch.out)[outcomes[c(2,3)]], collapse='-')
			)
			
			tmatch.out[, diffcols[1] ] <- tmatch.out[,outcomes[1]] - tmatch.out[,outcomes[2]]
			tmatch.out[, diffcols[2] ] <- tmatch.out[,outcomes[1]] - tmatch.out[,outcomes[3]]
			tmatch.out[, diffcols[3] ] <- tmatch.out[,outcomes[2]] - tmatch.out[,outcomes[3]]
			
			#Individual t-tests
			t1 <- t.test(x=tmatch.out[,outcomes[1]], y=tmatch.out[,outcomes[2]], paired=TRUE)
			t2 <- t.test(x=tmatch.out[,outcomes[1]], y=tmatch.out[,outcomes[3]], paired=TRUE)
			t3 <- t.test(x=tmatch.out[,outcomes[2]], y=tmatch.out[,outcomes[3]], paired=TRUE)
			
			ci <- data.frame(row.names=1:3)
			ci$Treatments <- diffcols
			ci$t <- c(t1$statistic, t2$statistic, t3$statistic)
			ci$df <- c(t1$parameter, t2$parameter, t3$parameter)
			ci$p.value <- c(t1$p.value, t2$p.value, t3$p.value)
			ci$sig <- star(ci$p.value)
			ci$mean.diff <- c(t1$estimate, t2$estimate, t3$estimate)
			ci$ci.min <- c(t1$conf.int[1], t2$conf.int[1], t3$conf.int[1])
			ci$ci.max <- c(t1$conf.int[2], t2$conf.int[2], t3$conf.int[2])
			
			results$t.tests <- ci
		}		
	}
	
	class(results) <- c("trimatch.summary", "list")
	return(results)
}

#' Prints the results of \code{\link{summary.triangle.matches}}.
#' 
#' This is an S3 generic function to print the results of 
#' \code{\link{summary.triangle.matches}}.
#' 
#' @param x results of \code{\link{summary.triangle.matches}}.
#' @param ... multiple results of \code{\link{summary.triangle.matches}}. These
#'        must be named. For example, \code{"Method 1" = summary(tmath, outcome)}.
#' @method print trimatch.summary
#' @export
print.trimatch.summary <- function(x, ...) {
	parms <- list(...)
	if(!missing(x)) {
		NextMethod(x) #TODO: Provide a better print method.
	} else {
		df <- data.frame(
			Method=character(),
			Friedman.chi2=numeric(), 
			Friedman.p=numeric(), 
			Friedman.sig=character(),
			rmANOVA.F=numeric(), 
			rmANOVA.p=numeric(), 
			rmANOVA.sig=character(),
			stringsAsFactors=FALSE)
		for(i in seq_along(parms)) {
			df <- rbind(df, data.frame(
				Method = names(parms)[i],
				Friedman.chi2 = parms[[i]]$friedman.test$statistic,
				Friedman.p = parms[[i]]$friedman.test$p.value,
				Friedman.sig = star(parms[[i]]$friedman.test$p.value),
				rmANOVA.F = parms[[i]]$rmanova[[1]]$F,
				rmANOVA.p = parms[[i]]$rmanova[[1]]$p,
				rmANOVA.sig = star(parms[[i]]$rmanova[[1]]$p),
				stringsAsFactors=FALSE))
		}
		names(df)[7] <- ' '
		names(df)[4] <- ' '
		row.names(df) <- 1:nrow(df)
		class(df) <- c('tritable','data.frame')
		return(df)
	}
}
