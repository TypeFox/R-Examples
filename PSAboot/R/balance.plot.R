#' Plot method for balance.
#' 
#' @param x results from \code{\link{balance}}
#' @param unadjusted.color color of the vertical line representing the mean 
#'        unadjusted effect size for all covariates.
#' @param complete.color color of the vertical line representing the mean adjusted
#'        effect size for all covariates using the complete dataset.
#' @param pooled.color color of the vertical line represeting the mean adjusted
#'        effect size for all covariates across all bootstrapped samples.
#' @param ... currently unused.
#' @method plot PSAboot.balance
#' @export
plot.PSAboot.balance <- function(x, 
								 unadjusted.color='red', 
								 complete.color='blue', 
								 pooled.color='black', 
								 ...) {
	df.complete <- x$complete
	df.complete <- reshape2::melt(apply(df.complete, 1, x$pool.fun, na.rm=TRUE))
	df.complete$color <- 'Complete'
	df.complete$X2 <- row.names(df.complete)
	df.pool <- x$pooled
	df.pool <- reshape2::melt(apply(df.pool, 2, x$pool.fun, na.rm=TRUE))
	df.pool$X2 <- row.names(df.pool)
	df.pool$color <- 'Pooled'
	df.unadj <- data.frame(method='Unadjusted', value=x$pool.fun(x$unadjusted))
	df.unadj$color <- 'Unadjusted'
	tmp <- reshape2::melt(x$pooled)
	names(tmp) <- c('Var1', 'X2', 'value')
	tmp$color <- 'Pooled'
	p <- ggplot(tmp) + 
		geom_vline(data=df.pool, aes(xintercept=value, color=color), alpha=1) +
		geom_vline(data=df.unadj, aes(xintercept=value, color=color), alpha=.75) + 
		geom_vline(data=df.complete, aes(xintercept=value, color=color), alpha=.75) +
		geom_density(aes(x=value), color='black') +
		scale_color_manual('Mean Balance', values=c('Pooled'=pooled.color,
													'Unadjusted'=unadjusted.color,
													'Complete'=complete.color)) +
		xlab('Balance (Effect Size)') + ylab('') + 
		xlim(c(0, 1.05 * max(c(x$pooled, .1, df.unadj$value)))) +
		facet_wrap(~ X2, ncol=1)
	return(p)
}
