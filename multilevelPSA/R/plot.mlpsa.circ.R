utils::globalVariables(c('mnx','mny','Diff','strata2','xmark','ymark','n','y','es'))

#' Plots the results of a multilevel propensity score model.
#'
#' The plot created uses the \code{ggplot2} framework. As such, additional modificaitons
#' can be made. This plot is an extension of the \code{circ.psa} function in the
#' \code{PSAgraphics} package for multilevel models.
#'
#' @param x the results of \code{\link{mlpsa}}.
#' @param xlab label for the x-axis.
#' @param ylab label for the y-axis.
#' @param legendlab the label for the legend, or NULL to exclude.
#' @param title title for the figure.
#' @param overall.col the color used for the overall results.
#' @param overall.ci.col the color used for the confidence intervals.
#' @param level1.plot logical value indicating whether level 1 points should be plotted.
#' @param level1.point.size the size of level 1 points
#' @param level1.rug.plot the placement for plotting a level 2 rug. Possible values
#'        are \code{bl} (for left and bottom), \code{tr} (for top and right), or
#'        NULL (to exclude).
#' @param level1.projection.lines logical value indicating whether level 1 project lines
#'        (parallel to the unit line) are drawn.
#' @param level2.plot logical value indicating whether level 2 points should be plotted.
#' @param level2.point.size the size of level 2 points
#' @param level2.rug.plot the placement for plotting a level 2 rug. Possible values
#'        are \code{bl} (for left and bottom), \code{tr} (for top and right), or
#'        NULL (to exclude).
#' @param level2.projection.lines logical value indicating whether level 2 project lines
#'        (parallel to the unit line) are drawn.
#' @param level2.label logical value indicating whether level 2 points should be labeled.
#' @param unweighted.means logical value indicating whether horizontal and vertical
#'        lines are drawn representing the unweighted (i.e. unadjusted from phase I
#'        of PSA) means for each level 2, or cluster.
#' @param weighted.means logical value indicating whether horizontal and vertical
#'        lines are drawn representing the weighted means for each level 2, or cluster.
#' @param fill.colors if specified, the colors to use for level 2 points.
#' @param ... currently unused.
#' @seealso plot.mlpsa
#' @export
#' @examples
#' \dontrun{
#' data(pisana)
#' data(pisa.colnames)
#' data(pisa.psa.cols)
#' mlctree = mlpsa.ctree(pisana[,c('CNT','PUBPRIV',pisa.psa.cols)], 
#'                       formula=PUBPRIV ~ ., level2='CNT')
#' student.party = getStrata(mlctree, pisana, level2='CNT')
#' student.party$mathscore = apply(student.party[,paste0('PV', 1:5, 'MATH')], 1, sum) / 5
#' results.psa.math = mlpsa(response=student.party$mathscore, 
#'        treatment=student.party$PUBPRIV, 
#'        strata=student.party$strata, 
#'        level2=student.party$CNT, minN=5)
#' mlpsa.circ.plot(results.psa.math, legendlab=FALSE)
#' }
mlpsa.circ.plot <- function(x,
		xlab=names(multilevelPSA$level2.summary)[4], 
		ylab=names(multilevelPSA$level2.summary)[5], 
		legendlab='Level 2', 
		title=NULL,
		overall.col="blue", 
		overall.ci.col='green',
		level1.plot=FALSE, 
		level1.point.size=NULL, 
		level1.rug.plot=NULL, 
		level1.projection.lines=FALSE,
		level2.plot=TRUE, 
		level2.point.size=NULL,
		level2.rug.plot='tr', 
		level2.projection.lines=TRUE,
		level2.label=FALSE, 
		unweighted.means=FALSE, 
		weighted.means=FALSE,
		fill.colors=NULL, 
		...
) {
	stopifnot(is.mlpsa(x))
	multilevelPSA = x
	ggplot.alpha <- function(...) get("alpha", grep("package:ggplot2$", search()))(...)

	if(missing(multilevelPSA)) {
		stop('Must provide multilevelPSA from multilevel.psa')
	}
	
	level1.summary = multilevelPSA$level1.summary
	level2.summary = multilevelPSA$level2.summary
	unweighted.summary = multilevelPSA$unweighted.summary
	plot.range = multilevelPSA$plot.range
	overall.ci = multilevelPSA$overall.ci
	overall.wtd = multilevelPSA$overall.wtd
	overall.mnx = multilevelPSA$overall.mnx
	overall.mny = multilevelPSA$overall.mny
	projection.intercept = multilevelPSA$projection.intercept
	
	#This is a bit of a hack. I renamed the mnx and mny columns in the mlpsa
	#function to use the treatment levels. This will duplicate those columns.
	level2.summary$mnx = multilevelPSA$level2.summary[,4]
	level2.summary$mny = multilevelPSA$level2.summary[,5]
	level2.summary$panel = 'Circular'
	
	#Setup ggplot2
	p = ggplot(level1.summary, aes_string(x=names(multilevelPSA$level2.summary)[4], 
										  y=names(multilevelPSA$level2.summary)[5]))
	p = p + coord_fixed(ratio=1) + 
			scale_x_continuous(limits=plot.range) +
			scale_y_continuous(limits=plot.range) +
			theme(axis.text=element_text(margin=ggplot2::unit(.1, "cm")))
	#Draw dashed lines for unweighted means
	if(unweighted.means) {
		p = p + geom_segment(data=unweighted.summary, 
							 aes_string(x=names(unweighted.summary)[3], 
							 	xend=names(unweighted.summary)[3], 
							 	yend=names(unweighted.summary)[2], 
							 	color='level2'), 
							 y=plot.range[1], 
							 alpha=.4, linetype='dashed', size=.5) +
				geom_segment(data=unweighted.summary, 
							 aes_string(y=names(unweighted.summary)[2], 
							 	xend=names(unweighted.summary)[3], 
							 	yend=names(unweighted.summary)[2], 
							 	color='level2'),
							 x=plot.range[1], alpha=.4, linetype='dashed', size=.5)
	}
	#Draw solid lines for weighted means
	if(weighted.means) {
		p = p + geom_segment(data=level2.summary, 
							 aes(x=mnx, xend=mnx, yend=mny, color=level2),
							 y=plot.range[1], alpha=.7, size=.5) +
				geom_segment(data=level2.summary, 
							 aes(y=mny, xend=mnx, yend=mny, color=level2),
							 x=plot.range[1], alpha=.7, size=.5)
	}
	#Rug plots
	if(!is.null(level1.rug.plot)) {
		p = p + level1.rug.plot(data=level1.summary, 
								sides = level1.rug.plot,
								aes_string(x=names(level1.summary)[5], 
									y=names(level1.summary)[4], color='level2'), 
								alpha=.5, size=1)
	}
	if(!is.null(level2.rug.plot)) {
		p = p + geom_rug(data=level2.summary, sides = level2.rug.plot,
						aes(x=mnx, y=mny, color=level2), alpha=.6, size=1)
	}
	#Projection lines
	if(level1.projection.lines) {
		p = p + geom_abline(data=level1.summary, 
							aes(intercept=Diff, slope=1, color=strata2), alpha=1, size=1)
	}
	if(level2.projection.lines) {
		tmp = level2.summary[order(level2.summary$diffwtd),]
		p = p + geom_segment(data=tmp, aes(x=mnx, y=mny, xend=xmark, yend=ymark, color=level2), 
							 size=1, alpha=1, linetype=1)
	}
	#Unit line
	p = p + geom_abline(slope=1, intercept=0, alpha=.7, size=1.4)
	#Overall multilevelPSA
	p = p + geom_abline(slope=1, intercept=overall.ci[1], 
						color=overall.ci.col, linetype=3, size=.6, alpha=.9) +
			geom_abline(slope=1, intercept=overall.ci[2], 
						color=overall.ci.col, linetype=3, size=.6, alpha=.9)
	#Overall difference line (parallel to the unit line)
	p = p + geom_abline(slope=1, intercept=overall.wtd, 
						color=overall.col, linetype='dashed', size=.6, alpha=.9)
	#Overall results (vertical line)
	p = p + geom_vline(xintercept=overall.mnx, 
					   color=overall.col, size=.6, alpha=.7) +
			geom_hline(yintercept=overall.mny, 
					   color=overall.col, size=.6, alpha=.7)
	#Point for each level 1 stratum
	if(level1.plot) {
		#TODO: WARNING can't seem to specify both size and fill for secondary data set 
		#(seems to be Windows only)
		if(is.null(level1.point.size)) {
			p = p + geom_point(data=level1.summary, 
							   aes_string(x=names(level1.summary)[5], 
								   	y=names(level1.summary)[4], 
								   	fill='level2', 
								   	size='n'), 
								   alpha=.6)
		} else {
			p = p + geom_point(data=level1.summary, 
							   aes_string(x=names(level1.summary)[5], 
								   	y=names(level1.summary)[4], 
								   	fill='level2'), 
							   size=1, alpha=.6)
		}
	}
	#Level 2 points
	if(level2.plot) {
		if(is.null(level2.point.size)) {
			p = p + geom_point(data=level2.summary, 
							   aes(x=mnx, y=mny, fill=level2, size=n), 
							   shape=21, color='black') 
		} else {
			p = p + geom_point(data=level2.summary, 
							   aes(x=mnx, y=mny, fill=level2, size=n), 
							   size=multilevelPSA$level2.point.size, shape=21, color='black')
		}
	}
	#Label level 2 points
	if(level2.label) { 
		p = p + geom_text(data=level2.summary, 
						  aes(x=mnx, y=mny, label=level2, hjust=.5, vjust=.5), 
						  stat='identity', size=4, color='black')
	}
	#Projected difference distribution
	p = p + geom_abline(slope=-1, 
						intercept=(projection.intercept - .03 * diff(plot.range)), 
						color='black', size=.5, alpha=.7)
	#Label the mean lines
	labeling = rbind(
		data.frame(x = plot.range[1] , y = overall.mny, 
				   label = prettyNum(overall.mny, digits=2)),
		data.frame(x = overall.mnx, y = plot.range[1], 
				   label = prettyNum(overall.mnx, digits=2))
	)
	p = p + geom_text(data=labeling[1,], aes(x=x, y=y, label=label), 
					  color=overall.col, vjust=-.5, hjust=-1, size=3)
	p = p + geom_text(data=labeling[2,], aes(x=x, y=y, label=label), 
					  color=overall.col, vjust=-.5, hjust=1.5, angle=-90, size=3)
	#Labels
 	p = p + xlab(xlab) + ylab(ylab) + scale_size_continuous('Size')
	#Difference disttribution (as x's)
	p = p + geom_point(data=level2.summary, 
					   aes(x=xmark, y=ymark, color=level2), # TODO: make the shape a parameter
					   stat='identity', size=4, shape=3, alpha=1)
	#Set color scheme and legend
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
	if(!is.null(title)) {
		p = p + theme(title=title)
	}
	
	return(p)
}
