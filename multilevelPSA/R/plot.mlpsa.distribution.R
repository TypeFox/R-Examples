#' Plots distribution for either the treatment or comparison group.
#' 
#' @param x the results of \code{\link{mlpsa}}.
#' @param treat the group to plot. This must be one of the two levels of the treatment variable.
#' @param fill.colours if specified, the colors to use for level 2 points.
#' @param flip if TRUE, the level 2 clusters will be on the y-axis and the outcome
#'        variable on the x-axis. Otherwise reversed.
#' @param label the label to use for the axis.
#' @param level2.label the axis label for the level 2 indicators.
#' @param legendlab the label for the legend, or NULL to exclude a legend.
#' @param axis.text.size the size of the axis text
#' @param fill.colors if specified, the colors to use for level 2 points.
#' @param ... currently unused.
#' @seealso plot.mlpsa
#' @export 
mlpsa.distribution.plot <- function(x, 
									treat, 
									fill.colours=NULL, 
									flip=TRUE, 
									label=treat, 
									level2.label=NULL, 
									legendlab=NULL, 
									axis.text.size=8,
									fill.colors=NULL, 
									...) {
	stopifnot(is.mlpsa(x))
	multilevelPSA = x
	
	if(is.na(treat) | !treat %in% names(multilevelPSA$level2.summary)[4:5]) {
		stop(paste('treat parameter must be specified. Possible values are ', 
				   names(multilevelPSA$level2.summary)[4], ' or ',
				   names(multilevelPSA$level2.summary)[4], sep=''))
	}
	
	level1.summary = multilevelPSA$level1.summary
	level2.summary = multilevelPSA$level2.summary
	unweighted.summary = multilevelPSA$unweighted.summary
	plot.range = multilevelPSA$plot.range
	overall = NA
	if(treat == multilevelPSA$x.label) {
		overall = multilevelPSA$overall.mnx
	}  else {
		overall = multilevelPSA$overall.mny
	}
	
	# This is possibly a bug with ggplot2. Smaller test reveal that ggplot2 has 
	# problems with column names of TRUE and FALSE.
	if(treat == 'TRUE') {
		names(level1.summary)[names(level1.summary) == 'TRUE'] <- 'Treatment'
		names(level2.summary)[names(level2.summary) == 'TRUE'] <- 'Treatment'
		treat <- 'Treatment'
	} else if(treat == 'FALSE') {
		names(level1.summary)[names(level1.summary) == 'FALSE'] <- 'Control'
		names(level2.summary)[names(level2.summary) == 'FALSE'] <- 'Control'
		treat <- 'Control'
	}
	
	xname = 'level2'
	yname = treat
	fillname = 'level2'
	
	if(flip) {
		xname = treat
		yname = 'level2'
		p = ggplot(level1.summary, aes_string(x=xname, y=yname))
		#This is a bit of a hack. I renamed the mnx and mny columns in the mlpsa
		#function to use the treatment levels. This will duplicate those columns.
		level2.summary$mnx = multilevelPSA$level2.summary[,multilevelPSA$x.lab]
		level2.summary$mny = multilevelPSA$level2.summary[,multilevelPSA$y.lab]
		p = p + scale_x_continuous(limits=plot.range)
		p = p + theme(legend.position='none', 
					 axis.text.y=element_text(size=axis.text.size, angle=0, hjust=.5))
		p = p + ylab(level2.label)+ xlab(label)
		p = p + geom_rug(data=level1.summary, aes_string(x=treat, y=NULL, colour=fillname), 
						 alpha=.6, size=.5)
		p = p + geom_vline(xintercept=overall, colour='blue', size=.6)
	} else {
		p = ggplot(level1.summary, aes_string(x=xname, y=yname))
		level2.summary$mnx = multilevelPSA$level2.summary[,multilevelPSA$x.lab]
		level2.summary$mny = multilevelPSA$level2.summary[,multilevelPSA$y.lab]
		p = p + scale_y_continuous(limits=plot.range)
		p = p + theme(legend.position=c(-1,-1), 
					 axis.text.x=element_text(size=axis.text.size, angle=-90, hjust=0, vjust=.5))
		p = p + ylab(label) + xlab(level2.label)
		p = p + geom_rug(data=level1.summary, aes_string(x=NULL, y=treat, colour=fillname), 
						 alpha=.6, size=.5)
		p = p + geom_hline(yintercept=overall, colour='blue', size=.6)
	}
	
	p = p + geom_point(stat='identity', alpha=.3, size=1.3)
	if(!is.null(fill.colors)) {
		p = p + scale_color_manual(guide='none', values=fill.colors) + 
			scale_fill_manual(guide='none', values=fill.colors)
	} else if(length(unique(level2.summary$level2)) > 20) {
		#No legend since the legend would be bigger than the plot
		p = p + scale_color_hue(guide='none') + scale_fill_hue(guide='none')
	} else if(length(unique(level1.summary$level2)) > 8) {
		p = p + scale_color_hue(legendlab) + scale_fill_hue(legendlab)
	} else {
		p = p + scale_color_brewer(legendlab, type='qual') + 
			scale_fill_brewer(legendlab, type='qual')
	}
	p = p + geom_point(data=level2.summary, aes_string(x=xname, y=yname, size='n', fill=fillname), 
					   stat='identity', shape=21, colour='black')
	if(flip) {
		labeling = level2.summary[,c(yname, treat)]
		names(labeling) = c('yname', 'label')
		p = p + geom_text(data=labeling, x=plot.range[1], 
						  aes(y=yname, label=prettyNum(label, digits=3, drop0trailing=FALSE)), 
						  size=2.5, hjust=0)
	} else {
		labeling = level2.summary[,c(xname, treat)]
		names(labeling) = c('xname', 'label')
		p = p + geom_text(data=labeling, y=plot.range[1], 
						  aes(x=xname, label=prettyNum(label, digits=3, drop0trailing=FALSE)), 
						  size=2.5, hjust=1, angle=-90)
	}
	return(p)
}
