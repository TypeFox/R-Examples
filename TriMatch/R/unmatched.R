#' Returns rows from \code{\link{trips}} that were not matched by \code{\link{trimatch}}.
#' 
#' This function returns a subset of \code{\link{trips}} that were not matched by
#' \code{\link{trimatch}}. All data frame methods work with the returned object
#' but special \code{summary} function will provided relevant information.
#' 
#' @param tmatch the results of \code{\link{trimatch}}.
#' @return a data frame of unmatched rows.
#' @export
unmatched <- function(tmatch) {
	tpsa <- attr(tmatch, 'triangle.psa')
	unmatched <- tpsa[!(tpsa$id %in% c(tmatch[,1], tmatch[,2], tmatch[,3])),]
	class(unmatched) <- c('unmatched', 'triangle.psa', 'data.frame')
	attr(unmatched, 'triangle.psa') <- tpsa
	return(unmatched)
}

#' Provides a summary of unmatched subjects.
#' 
#' Will return as a list and print the percentage of total unmatched rows and
#' percent by treatment.
#' 
#' @param object results of \code{\link{unmatched}}
#' @param digits number of digits to print.
#' @param ... currently unused.
#' @return a list of summary results.
#' @export
#' @method summary unmatched
summary.unmatched <- function(object, digits=3, ...) {
	tpsa <- attr(object, 'triangle.psa')
	results <- list()
	results$PercentUnmatched <- nrow(object) / nrow(tpsa)
	results$UnmatchedByTreat <- table(object$treat)
	results$PercentByTreat <- table(object$treat) / table(tpsa$treat)
	results$Unmatched <- table(object$treat)
	cat(paste0(nrow(object), ' (',
			prettyNum(nrow(object) / nrow(tpsa) * 100, digits=digits), 
			'%) of ', nrow(tpsa), ' total data points were not matched.'))
	cat('\nUnmatched by treatment:\n')
	tmp <- paste0(results$UnmatchedByTreat, ' (', prettyNum(results$UnmatchedByTreat / 
							table(tpsa$treat) * 100, digits=digits), '%)')
	names(tmp) <- names(results$UnmatchedByTreat)
	print(tmp, quote=FALSE)
	invisible(results)
}
