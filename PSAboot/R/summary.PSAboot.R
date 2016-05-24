#' Summary of pooled results from PSAboot
#' 
#' @param object result of \code{\link{PSAboot}}.
#' @param ... currently unused.
#' @return a list with pooled summary statistics.
#' @method summary PSAboot
#' @export
summary.PSAboot <- function(object, ...) {
	sum <- list()
	bal <- balance(object, ...)
	for(i in unique(object$pooled.summary$method)) {
		sum2 <- list()
		rows <- object$pooled.summary[object$pooled.summary$method == i,]
		
		sig.pos <- rows$ci.min > 0
		sig.neg <- rows$ci.max < 0
		sum2[['sig.pos.per']] <- prop.table(table(factor(sig.pos, levels=c('TRUE','FALSE')))) * 100
		sum2[['sig.neg.per']] <- prop.table(table(factor(sig.neg, levels=c('TRUE','FALSE')))) * 100
		sum2[['sig.tot.per']] <- prop.table(table(factor(sig.pos | sig.neg,
											   levels=c('TRUE','FALSE')))) * 100
		
		m <- mean(rows$estimate, na.rm=TRUE)
		wm <- weighted.mean(rows[!is.na(rows$estimate),]$estimate, 
							# HACK: Occassionally there will be NA estimates
							1 / apply(bal$balances[[i]], 1, mean), 
							na.rm=TRUE)
		ci <- c(m - 2 * sd(rows$estimate, na.rm=TRUE),
				m + 2 * sd(rows$estimate, na.rm=TRUE) )
		
		complete <- object$complete.summary[object$complete.summary$method == i,]
		
		sum2[['bootstrap.mean']] <- m
		sum2[['bootstrap.ci']] <- ci
		sum2[['bootstrap.weighted.mean']] <- wm
		sum2[['percent.sig']] <- table(sig.pos | sig.neg)
		sum2[['complete']] <- complete
		sum[[i]] <- sum2
	}
	
	class(sum) <- c('PSAbootSummary', 'list')
	return(sum)
}

#' Print method for PSAboot Summary.
#' 
#' @param x result of \code{\link{summary.PSAboot}}
#' @param digits desired number of digits after the decimal point.
#' @param ... unused.
#' @method print PSAbootSummary
#' @export
print.PSAbootSummary <- function(x, digits=3, ...) {
	for(i in names(x)) {
		sum2 <- x[[i]]
		complete <- x[[i]][['complete']]
		m <- x[[i]][['bootstrap.mean']]
		wm <- x[[i]][['bootstrap.weighted.mean']]
		ci <- x[[i]][['bootstrap.ci']]
		sig.tot.per <- x[[i]][['sig.tot.per']]
		sig.pos.per <- x[[i]][['sig.pos.per']]
		sig.neg.per <- x[[i]][['sig.neg.per']]
		cat(paste0(i, ' Results:'))
		cat(paste0('\n   Complete estimate = ', prettyNum(complete$estimate, digits=digits)))
		cat(paste0('\n   Complete CI = [', prettyNum(complete$ci.min, digits=digits), ', ',
				   prettyNum(complete$ci.max, digits=digits), ']'))
		cat(paste0('\n   Bootstrap pooled estimate = ', prettyNum(m, digits=digits), 
				   '\n   Bootstrap weighted pooled estimate = ', prettyNum(wm, digits=digits),
				   '\n   Bootstrap pooled CI = [', prettyNum(ci[1], digits=digits), ', ', 
				   prettyNum(ci[2], digits=digits), ']\n'))
		
		cat(paste0('   ',
				   prettyNum(unname(sig.tot.per['TRUE']), 
				   		  digits=digits),
				   '% of bootstrap samples have confidence intervals that do not span zero.\n',
				   '      ', prettyNum(unname(sig.pos.per['TRUE']), digits=digits), '% positive.\n',
				   '      ', prettyNum(unname(sig.neg.per['TRUE']), digits=digits), '% negative.\n'))
	}
}

#' Convert the results of PSAboot summary to a data frame.
#' 
#' @param x results of \code{\link{summary.PSAboot}}
#' @param row.names row names.
#' @param optional unused.
#' @param ... unused.
#' @method as.data.frame PSAbootSummary
#' @export
as.data.frame.PSAbootSummary <- function(x, row.names=NULL, optional=FALSE, ...) {
	df <- data.frame()
	for(i in names(x)) {
		complete <- x[[i]]$complete
		df <- rbind(df, data.frame(
			method=i, 
			bootstrap.estimate=x[[i]]$bootstrap.mean,
			bootstrap.ci.min=x[[i]]$bootstrap.ci[1],
			bootstrap.ci.max=x[[i]]$bootstrap.ci[2],
			complete.estimate=complete[1,]$estimate,
			complete.ci.min=complete[1,]$ci.min,
			complete.ci.max=complete[1,]$ci.max,
			stringsAsFactors=FALSE
			))
	}
	return(df)
}
