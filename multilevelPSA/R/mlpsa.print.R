#' Returns true if the object is of type \code{mlpsa}
#'
#' @param x the object to test
#' @export
is.mlpsa <- function(x) {
	inherits(x, 'mlpsa')
}

#' Prints basic information about a \code{mlpsa} class.
#'
#' @param x the \code{mlpsa} class.
#' @param ... unused.
#' @method print mlpsa
#' @export
print.mlpsa <- function(x, ...) {
	#TODO: Create a better print function
	cat('The following fields are available:\n')
	print(ls(x))
}

#' Provides a summary of a \code{mlpsa} class.
#'
#' @param object the mlpsa object.
#' @param overall.label the label to place in the strata column for the overall results.
#' @param ... unused.
#' @method summary mlpsa
#' @export
summary.mlpsa <- function(object, overall.label='Overall', ...) {
	message(paste("Multilevel PSA Model of ", nrow(object$level1.summary), ' strata for ',
			  nrow(object$level2.summary), ' levels.\n',
		'Approx t: ', round(object$approx.t, digits=2), '\n',
		'Confidence Interval: ', round(object$overall.ci[1], digits=2), ', ', round(object$overall.ci[2], digits=2),
		'\n\n',
	sep=''))
	
	l1 <- object$level1.summary
	l2 <- object$level2.summary
	
	xdf <- data.frame(level2=character(), strata=character(),
					  Treat=numeric(), Treat.n=integer(),
					  Control=numeric(), Control.n=integer(), 
					  ci.min=numeric(), ci.max=numeric(),
					  stringsAsFactors=FALSE)
	for(i in 1:nrow(l2)) {
		l2.row <- l2[i,]
		l1.rows <- l1[l1$level2 == l2.row$level2[1],]
		xdf <- rbind(xdf, data.frame(
			level2 = l2.row$level2,
			strata = overall.label,
			Treat = l2.row[,object$x.label],
			Treat.n = l2.row[,paste(object$x.label, '.n', sep='')],
			Control = l2.row[,object$y.label],
			Control.n = l2.row[,paste(object$y.label, '.n', sep='')],
			ci.min = l2.row$ci.min,
			ci.max = l2.row$ci.max,
			stringsAsFactors=FALSE
		))
		for(j in 1:nrow(l1.rows)) {
			l1.row <- l1.rows[j,]
			xdf <- rbind(xdf, data.frame(
				level2 = NA,
				strata = j,
				Treat = l1.row[,object$x.label],
				Treat.n = l1.row[,paste(object$x.label, '.n', sep='')],
				Control = l1.row[,object$y.label],
				Control.n = l1.row[,paste(object$y.label, '.n', sep='')],
				ci.min = NA,
				ci.max = NA,
				stringsAsFactors=FALSE
			))			
		}
	}
	return(xdf)
}

