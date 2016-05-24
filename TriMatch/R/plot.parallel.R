utils::globalVariables(c('Outcome','ID','group1'))

#' Parallel coordinate plot for the three groups and dependent variable.
#' 
#' Creates a \code{ggplot2} figure of a parallel coordinate plot.
#' 
#' @param tmatch results from \code{\link{trimatch}}.
#' @param outcome vector of the outcome
#' @export
parallel.plot <- function(tmatch, outcome) {
	tmatch <- merge(tmatch, outcome)
	outcomes <- grep("*.out", names(tmatch))
	tmatch$id <- 1:nrow(tmatch)
	out <- melt(tmatch[,c(outcomes, ncol(tmatch))],id.vars='id')
	names(out) <- c('ID','Treatment','Outcome')
	out2 <- describeBy(out$Outcome, group=out$Treatment, mat=TRUE)
	p <- ggplot(out) + 
		geom_line(aes(x=Treatment, y=Outcome, group=ID, colour=ID), alpha=.1) +
		theme(legend.position='none') + scale_colour_gradient() +
		geom_point(data=out2, aes(x=group1, y=mean), color='red', shape=20, size=4)
	return(p)
}
