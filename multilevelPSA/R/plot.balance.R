utils::globalVariables(c('value','covariate','EffectSize','level2','group'))

#' Multiple covariate blance assessment plot.
#' 
#' A graphic based upon \code{\link{cv.bal.psa}} function in the \code{PSAgraphics}
#' package. This graphic plots the effect sizes for multiple covariated before and
#' after propensity score andjustement.
#'
#' @param x results of \code{\link{covariate.balance}}.
#' @param plot.strata whether individual strata should be plotted.
#' @param order how to order the y-axis. Possible values are adjusted,
#'        unadjusted, or NULL (don't reorder).
#' @param strata.size text size for strata if plotted.
#' @param strata.legend.guide guide for legend placement for strata.
#' @param point.size size of the overall effect size points.
#' @param point.alpha transparency level of the overall effect size points.
#' @param line.color the color of the line connecting the overall effect size ponts.
#' @param line.alpha transparency level of the line connecting the overall effect size points.
#' @param legend.position where to position the legend.
#' @param ... currently unused.
#' @return a ggplot2 with an attribute, \code{effects}, that is the data frame
#'        used to create the plot.
#' @method plot covariate.balance
#' @export
plot.covariate.balance <- function(x, 
								   plot.strata=FALSE, 
								   order=c('unadjusted','adjusted'),
								   strata.size=3,
								   strata.legend.guide='none', 
								   point.size=3,
								   point.alpha=1,
								   line.color='black',
								   line.alpha=.2,
								   legend.position=c(.8,.2),
								   ...) {
	plot.weighted <- FALSE
	if('plot.weighted' %in% names(list(...))) {
		plot.weighted <- list(...)[['plot.weighted']]
	}
	
	bal <- x$effects[,c('covariate',ifelse(plot.weighted, 'es.adj.wtd','es.adj'),'es.unadj')]
	strata <- x$strata.effects
	diff <- bal$es.unadj - bal$es.adj
	cvorder <- bal[order(diff, decreasing=FALSE),]$covariate
	tmp <- melt(bal, id='covariate', variable_name='EffectSize')
	tmp$covariate <- factor(tmp$covariate, levels=cvorder, ordered=TRUE)
	tmp$EffectSize <- as.character(tmp$EffectSize)
	
	if(!is.null(order)) {
		if(order[1] == 'adjusted') {
			ordering <- tmp[tmp$EffectSize == 'es.adj',]
		} else if(order[1] == 'unadjusted') {
			ordering <- tmp[tmp$EffectSize == 'es.unadj',]
		}
		ordering <- ordering[order(ordering$value),]$covariate
		tmp$covariate <- factor(tmp$covariate, levels=ordering, ordered=TRUE)
	}
	
	p <- ggplot(tmp, aes(x=value, y=covariate, group=covariate))
	if(plot.strata) {
		strata.effects <- data.frame()
		for(i in seq_along(strata)) {
			cov <- names(strata)[i]
			eff <- strata[[i]]
			eff$covariate <- cov
			strata.effects <- rbind(strata.effects, eff)
		}
		strata.effects$covariate <- factor(strata.effects$covariate, levels=cvorder, ordered=TRUE)
		p <- p + geom_text(data=strata.effects, aes(x=abs(es), y=covariate, 
							color=level2, label=strata), size=2, alpha=.3)
	}
	
	p <- p + geom_line(alpha=line.alpha, color=line.color)
	if(plot.strata) {
		p <- p + geom_point(size=point.size, alpha=point.alpha, aes(shape=EffectSize))
		p <- p + scale_color_hue(guide=strata.legend.guide)
	} else {
		p <- p + geom_point(size=strata.size, aes(shape=EffectSize, color=EffectSize))
		p <- p + scale_color_hue(' ', labels=c('Adjusted', 'Unadjusted'))
	}
	p <- p + scale_shape(' ', labels=c('Adjusted', 'Unadjusted')) +
		theme(legend.position=legend.position) +
		ylab('') + xlab('Effect Size')
	
	return(p)
}
