utils::globalVariables(c('Diff','ci.min','ci.max','n','ci.min.adjust','ci.max.adjust'))

#' Creates a graphic summarizing the differences between treatment and comparison
#' groups within and across level two clusters.
#'
#' @param x the results of \code{\link{mlpsa}}.
#' @param xlab label for the x-axis, or NULL to exclude.
#' @param ylab label for the y-aixs, or NULL to exclude.
#' @param title title of the figure, or NULL to exclude.
#' @param overall.col the color of the overall results line.
#' @param overall.ci.col the color of the overall confidence interval.
#' @param level2.point.size the point size of level 2 points.
#' @param level1.points logical value indicating whether level 1 strata should be plotted.
#' @param errorbars logical value indicating whether error bars should be plotted for
#'        for each level 1.
#' @param errorbars.adjusted.ci whether the Bonferonni adjusted error bars should
#'        be plotted (these will be dashed lines).
#' @param level2.rug.plot logical value indicating whether a rug plot should be
#'        plotted for level 2.
#' @param jitter logical value indicating whether level 1 points should be jittered.
#' @param reorder logical value indicating whether the level two clusters should be
#'        reordered from largest difference to smallest.
#' @param labelLevel2 logical value indicating whether the difference for each level 2
#'        should be labeled.
#' @param sd If specified, effect sizes will be plotted instead of difference in the
#'        native unit.
#' @param xlim the limits of the x-axis.
#' @param ... currently unused.
#' @seealso plot.mlpsa
#' @export
#' @examples
#' \dontrun{
#' data(pisana)
#' data(pisa.colnames)
#' data(pisa.psa.cols)
#' mlctree = mlpsa.ctree(pisana[,c('CNT','PUBPRIV',pisa.psa.cols)], formula=PUBPRIV ~ ., level2='CNT')
#' student.party = getStrata(mlctree, pisana, level2='CNT')
#' student.party$mathscore = apply(student.party[,paste0('PV', 1:5, 'MATH')], 1, sum) / 5
#' results.psa.math = mlpsa(response=student.party$mathscore, 
#'        treatment=student.party$PUBPRIV, 
#'        strata=student.party$strata, 
#'        level2=student.party$CNT, minN=5)
#' mlpsa.difference.plot(results.psa.math, sd=mean(student.party$mathscore, na.rm=TRUE))
#' }
mlpsa.difference.plot <- function(x,
		xlab,
		ylab=NULL,
		title=NULL,
		overall.col="blue",
		overall.ci.col='green',
		level2.point.size=NULL,
		level1.points=TRUE,
		errorbars=TRUE,
		errorbars.adjusted.ci=TRUE,
		level2.rug.plot=TRUE,
		jitter=TRUE,
		reorder=TRUE,
		labelLevel2=TRUE,
		sd=NULL,
		xlim,
		...
) {
	stopifnot(is.mlpsa(x))
	multilevelPSA = x
	#ggplot.alpha <- function(...) get("alpha", grep("package:ggplot2$", search()))(...)

	if(reorder) {
		multilevelPSA$level2.summary = multilevelPSA$level2.summary[
			order(multilevelPSA$level2.summary$diffwtd),]
		ord.level2 = multilevelPSA$level2.summary$level2[
			order(multilevelPSA$level2.summary$diffwtd)]
		multilevelPSA$level1.summary$level2 = factor(multilevelPSA$level1.summary$level2, 
													 levels=ord.level2)
		multilevelPSA$level2.summary$level2 = factor(multilevelPSA$level2.summary$level2, 
													 levels=ord.level2)
	}

	if(missing(xlab)) {
		if(is.null(sd)) {
			xlab <- 'Difference Score '
		} else {
			xlab <- 'Effect Size'
		}
		if(TRUE == all.equal( (multilevelPSA$overall.mnx - multilevelPSA$overall.mny),
					 multilevelPSA$overall.wtd)) {
			xlab <- paste0(xlab, ' (', multilevelPSA$x.label, ' - ',
						   multilevelPSA$y.label, ')')
		} else if(TRUE == all.equal( (multilevelPSA$overall.mny - multilevelPSA$overall.mnx),
							 multilevelPSA$overall.wtd)) {
			xlab <- paste0(xlab, ' (', multilevelPSA$y.label, ' - ',
						   multilevelPSA$x.label, ')')
		} else {
			warning('Cannot determine subtraction order.')
		}
	}
	
	if(!is.null(sd)) {
		multilevelPSA$level1.summary$Diff = multilevelPSA$level1.summary$Diff / sd
		multilevelPSA$level2.summary$diffwtd = multilevelPSA$level2.summary$diffwtd / sd
		multilevelPSA$level2.summary$ci.min = multilevelPSA$level2.summary$ci.min / sd
		multilevelPSA$level2.summary$ci.max = multilevelPSA$level2.summary$ci.max / sd
		multilevelPSA$level2.summary$ci.min.adjust = multilevelPSA$level2.summary$ci.min.adjust / sd
		multilevelPSA$level2.summary$ci.max.adjust = multilevelPSA$level2.summary$ci.max.adjust / sd
		multilevelPSA$level1.summary$Diff = multilevelPSA$level1.summary$Diff / sd
		multilevelPSA$overall.ci = multilevelPSA$overall.ci / sd
		multilevelPSA$overall.wtd = multilevelPSA$overall.wtd /sd
		multilevelPSA$plot.range = multilevelPSA$plot.range / sd
	}
	
	p = ggplot(multilevelPSA$level1.summary, aes(x=level2, y=Diff)) + coord_flip() +
			geom_hline(aes(yintercept=0), colour='black', size=1, alpha=.7) +
			geom_hline(yintercept=multilevelPSA$overall.wtd, colour=overall.col, size=1) + 
			geom_hline(yintercept=multilevelPSA$overall.ci, colour=overall.ci.col, size=1) + 
			theme(axis.text=element_text(margin=ggplot2::unit(0, "cm")), 
				  axis.text.y=element_text(size=8, angle=0, hjust=.5))
	if(errorbars.adjusted.ci) {
		p = p + geom_errorbar(data=multilevelPSA$level2.summary, 
							  aes(x=level2, y=NULL, ymin=ci.min.adjust, ymax=ci.max.adjust), 
							  colour='green', alpha=.6, linetype=2)		
	}
	if(errorbars) {
		p = p + geom_errorbar(data=multilevelPSA$level2.summary, 
							  aes(x=level2, y=NULL, ymin=ci.min, ymax=ci.max), 
							  colour='green', alpha=.6, linetype=1)
	}
	if(level1.points) {
		if(jitter) {
			p = p + geom_point(stat='identity', alpha=.3, size=.8, position='jitter')
		} else {
			p = p + geom_point(stat='identity', alpha=.3, size=.8)
		}
	}
	p = p + geom_point(data=multilevelPSA$level2.summary, aes(x=level2, y=diffwtd, size=n), 
					   fill='blue', alpha=.6, stat='identity', shape=21, colour='black')
	if(level2.rug.plot) {
		p = p + geom_rug(data=multilevelPSA$level2.summary, aes(x=NULL, y=diffwtd), 
						 alpha=.6, size=.5, colour='blue')
	}
	p = p + xlab(ylab) + ylab(xlab) + scale_size_continuous('Size')
	if(!is.null(title)) {
		p = p + ggtitle(title)
	}
	
	if(!missing(xlim)) {
		p <- p + ylim(xlim)
	}
	
	if(labelLevel2) {
		if(!missing(xlim)) {
			labelPos <- min(c(multilevelPSA$level2.summary$ci.min,
							  multilevelPSA$level1.summary$Diff,
							  xlim))
		} else {
			labelPos <- min(c(multilevelPSA$level2.summary$ci.min,
							  multilevelPSA$level1.summary$Diff))
		}
			#.1 * (max(multilevelPSA$level2.summary$ci.max) - min(multilevelPSA$level2.summary$ci.min))
		p = p + geom_text(data=multilevelPSA$level2.summary, aes(x=level2, 
						label=prettyNum(diffwtd, digits=2, drop0trailing=FALSE)), 
						y=labelPos, size=3, hjust=0)
	}
			
	return(p)
}
