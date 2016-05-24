utils::globalVariables(c('estimate','method','bootstrap.estimate','bootstrap.ci.min',
						 'bootstrap.ci.max','variable'))

#' Boxplot of PSA boostrap results.
#' 
#' @param x result of \code{\link{PSAboot}}.
#' @param tufte use Tufte's boxplot style. Requires the \code{ggthemes} package.
#' @param coord.flip Whether to flip the coordinates.
#' @param bootstrap.mean.color the color of the point for the boostrap mean, or NA
#'        to omit.
#' @param bootstrap.ci.color the color of the confidence intervals of the bootstrap
#'        samples, or NA to omit.
#' @param bootstrap.ci.width the width of the confidence interval lines at the end.
#' @param bootstrap.ci.size the size of the confidence interval lines.
#' @param overall.mean.color the color of the point for the overall (before bootstrapping)
#'        mean, or NA to omit.
#' @param ... unused
#' @method boxplot PSAboot
#' @export
boxplot.PSAboot <- function(x,
							bootstrap.mean.color='blue',
							bootstrap.ci.color='green',
							bootstrap.ci.width=0.5,
							bootstrap.ci.size=3,
							overall.mean.color='red',
							tufte=FALSE, 
							coord.flip=TRUE,
							...) {
	sum <- as.data.frame(summary(x))
	pooled.mean <- mean(x$pooled.summary$estimate, na.rm=TRUE)
	pooled.sd <- sd(x$pooled.summary$estimate, na.rm=TRUE)
	pooled.ci <- c(ci.min=pooled.mean - (qnorm(0.975) * pooled.sd/sqrt(x$M)),
				   ci.max=pooled.mean + (qnorm(0.975) * pooled.sd/sqrt(x$M)))
	p <- ggplot(x$pooled.summary, aes(y=estimate, x=method)) +
		geom_hline(yintercept=0, alpha=.5, size=2)
	if(!is.na(bootstrap.ci.color)) {
		p <- p + geom_hline(yintercept=pooled.ci, color=bootstrap.ci.color)
	}
	if(!is.na(bootstrap.mean.color)) {
		p <- p + geom_hline(yintercept=pooled.mean, color=bootstrap.mean.color)
	}
	if(!is.na(bootstrap.ci.color)) {
		p <- p + geom_errorbar(data=sum, aes(x=method, y=bootstrap.estimate, ymin=bootstrap.ci.min, 
					ymax=bootstrap.ci.max), color=bootstrap.ci.color, 
							   width=bootstrap.ci.width, size=bootstrap.ci.size)
	}
	if(tufte) {
		p <- p + geom_tufteboxplot()		
	} else {
		p <- p + geom_boxplot(alpha=.5)		
	}
	if(!is.na(bootstrap.mean.color)) {
		p <- p + geom_point(data=sum, aes(y=bootstrap.estimate, x=method), 
				   color=bootstrap.mean.color, size=5, alpha=.5)
	}
	if(!is.na(overall.mean.color)) {
		p <- p + geom_point(data=x$complete.summary, aes(y=estimate, x=method), 
				   color=overall.mean.color, size=3, alpha=.5)
	}
	if(coord.flip) {
		p <- p + coord_flip()
	}
	p <- p + xlab('')
	return(p)
}
