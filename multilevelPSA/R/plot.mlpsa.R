#' Plots the results of a multilevel propensity score model.
#'
#' The plot created uses the \code{ggplot2} framework. As such, additional modificaitons
#' can be made. This plot is an extension of the \code{circ.psa} function in the
#' \code{PSAgraphics} package for multilevel models.
#'
#' @param x the results of \code{\link{mlpsa}}.
#' @param ratio the ratio of the size of the distribution plots (left and bottom)
#'        to the circular plot.
#' @param plotExtra a plot to place in the lower left corner.
#' @param ... parameters passed to \code{\link{mlpsa.circ.plot}} and 
#'        \code{\link{mlpsa.distribution.plot}}
#' @method plot mlpsa
#' @export
#' @examples
#' \dontrun{
#' require(multilevelPSA)
#' require(party)
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
#' plot(results.psa.math)
#' }
plot.mlpsa <- function(x, ratio=c(1,2), plotExtra=NULL, ...) {
	stopifnot(is.mlpsa(x))
	mlpsa = x
	
	pcirc = mlpsa.circ.plot(mlpsa, legendlab=FALSE, ...) + 
				theme(legend.position='none') +
				xlab(NULL) + ylab(NULL)
	px = mlpsa.distribution.plot(mlpsa, treat=mlpsa$x.lab,
						flip=TRUE, label=mlpsa$x.lab, ...) +
				theme(legend.position='none')#, axis.text.x=element_blank())
	py = mlpsa.distribution.plot(mlpsa, treat=mlpsa$y.lab,
						flip=FALSE, label=mlpsa$y.lab, ...) +
				theme(legend.position='none')#, axis.text.y=element_blank())
	
	grid_layout = grid.layout(nrow=2, ncol=2, widths=c(ratio[1:2]), heights=ratio[2:1], respect=TRUE)
	grid.newpage()
	pushViewport( viewport( layout=grid_layout ) )
	align.plots(grid_layout, list(pcirc, 1, 2), list(px, 2, 2), list(py, 1, 1))
	if(!is.null(plotExtra)) {
		pushViewport(viewport(layout.pos.row=2, layout.pos.col=1, just=c("center", "center")))
		grid.draw(ggplotGrob(plotExtra))
	}
}
