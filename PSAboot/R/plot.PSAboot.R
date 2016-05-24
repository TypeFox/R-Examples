utils::globalVariables(c('iter','estimate','sig','bootstrap.estimate','bootstrap.ci.min',
						 'bootstrap.ci.max','value','color'))

#' Plot the results of PSAboot
#' 
#' @param x result of \code{\link{PSAboot}}.
#' @param sort how the sort the rows by mean difference. Options are to sort
#'        using the mean difference from matching, stratificaiton, both 
#'        individually, or no sorting.
#' @param ci.sig.color the color used for confidence intervals that do not span zero.
#' @param plot.overall whether to plot vertical lines for the overall (non-bootsrapped)
#'        estimate and confidence interval.
#' @param plot.bootstrap whether to plot vertical lines for the bootstrap pooled
#'        estimate and confidence interval.
#' @param ... currently unused.
#' @method plot PSAboot
#' @export
plot.PSAboot <- function(x, sort='all', ci.sig.color='red', 
						 plot.overall=FALSE,
						 plot.bootstrap=TRUE,
						 ...) {
	results <- x$pooled.summary
	results$sig <- results$ci.min > 0 | results$ci.max < 0
	ci.min <- mean(results$estimate, na.rm=TRUE) - 2 * sd(results$estimate, na.rm=TRUE)
	ci.max <- mean(results$estimate, na.rm=TRUE) + 2 * sd(results$estimate, na.rm=TRUE)
	
	if(sort == 'all') {
		results <- results[order(results$estimate),]
		for(i in unique(results$method)) {
			rows <- which(results$method == i)
			results[rows,]$iter <- 1:length(rows)
		}
	} else if(sort %in% unique(results$method)) {
		results.estimate <- reshape2::dcast(results[,c('iter','method','estimate')], 
								 iter ~ method, value='estimate')
		results.ci.min <- reshape2::dcast(results[,c('iter','method','ci.min')], 
							   iter ~ method, value='ci.min')
		results.ci.max <- reshape2::dcast(results[,c('iter','method','ci.max')], 
							   iter ~ method, value='ci.max')		
		o <- order(results.estimate[,sort])
		results.estimate <- results.estimate[o,]
		results.ci.min <- results.ci.min[o,]
		results.ci.max <- results.ci.max[o,]
		results.estimate$iter <- results.ci.min$iter <- 
			results.ci.max$iter <- 1:nrow(results.estimate)
		results.estimate <- reshape2::melt(results.estimate, id='iter')
		results.ci.min <- reshape2::melt(results.ci.min, id='iter')
		results.ci.max <- reshape2::melt(results.ci.max, id='iter')
		results <- cbind(results.estimate, results.ci.min$value, results.ci.max$value)
		names(results) <- c('iter', 'estimate','method','ci.min','ci.max')
		results$sig <- results$ci.min > 0 | results$ci.max < 0
		row.names(results) <- 1:nrow(results)
	}
	
	bootsum <- as.data.frame(summary(x))
	
	p <- ggplot(results, aes(y=iter, xmin=ci.min, xmax=ci.max, x=estimate, color=sig)) +
		geom_vline(xintercept=0, size=1.5, alpha=.25) + 
		geom_errorbarh(height=0, alpha=.5) + 
		geom_point(color='blue')
	if(plot.overall) {
		p <- p + geom_vline(data=x$complete.summary, aes(xintercept=estimate), 
				   color='blue', linetype=2) + 
 		geom_vline(data=x$complete.summary, aes(xintercept=ci.min), 
 				   color='green', linetype=2) +
 		geom_vline(data=x$complete.summary, aes(xintercept=ci.max), 
 				   color='green', linetype=2)
	}
	if(plot.bootstrap) {
		p <- p + geom_vline(data=bootsum, aes(xintercept=bootstrap.estimate), 
				   color='blue', linetype=1) +
		geom_vline(data=bootsum, aes(xintercept=bootstrap.ci.min), 
				   color='green', linetype=1) +
		geom_vline(data=bootsum, aes(xintercept=bootstrap.ci.max), 
				   color='green', linetype=1)
	}
	p <- p + scale_y_continuous() + 
		scale_color_manual(values=c("TRUE"=ci.sig.color, "FALSE"="grey")) +
		theme(legend.position='none', axis.ticks.y=element_blank(), 
			  axis.text.y=element_blank()) +
		xlab('Mean Difference') + ylab('') +
		facet_wrap(~ method, nrow=1)
	return(p)
}
