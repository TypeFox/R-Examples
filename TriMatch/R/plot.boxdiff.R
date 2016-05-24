utils::globalVariables(c('Treatments','Difference','V1','V2','estimate'))

#' Returns a \code{ggplot2} box plot of the differences.
#' 
#' A boxplot of differences between each pair of treatments.
#' 
#' @param tmatch the results from \code{\link{trimatch}}.
#' @param out a vector of the outcome measure of interest.
#' @param plot.mean logical indicating whether the means should be plotted.
#' @param ordering specify the order for doing the paired analysis, that is
#'        analysis will be conducted as:
#'        \code{ordering[1] - ordering[2]}, \code{ordering[1] - ordering[3]},
#'        and \code{ordering[2] - ordering[3]}.
#' @param ci.width the width for the confidence intervals.
#' @return a \code{ggplot2} boxplot of the differences.
#' @export
boxdiff.plot <- function(tmatch, out, plot.mean=TRUE, 
						 ordering = attr(tmatch, 'match.order'),
						 ci.width=.5) {
	tmatch.out <- merge(tmatch, out)
	outcomes <- sapply(ordering, function(x) { which(names(tmatch.out) == paste0(x, '.out')) })
	tmatch.out$id <- 1:nrow(tmatch.out)
	
	diffcols <- c(
		paste(names(tmatch.out)[outcomes[c(1,2)]], collapse='-'),
		paste(names(tmatch.out)[outcomes[c(1,3)]], collapse='-'),
		paste(names(tmatch.out)[outcomes[c(2,3)]], collapse='-')
	)
	
	tmatch.out[, diffcols[1] ] <- tmatch.out[,outcomes[1]] - tmatch.out[,outcomes[2]]
	tmatch.out[, diffcols[2] ] <- tmatch.out[,outcomes[1]] - tmatch.out[,outcomes[3]]
	tmatch.out[, diffcols[3] ] <- tmatch.out[,outcomes[2]] - tmatch.out[,outcomes[3]]
	
	#Individual t-tests
	t1 <- t.test(x=tmatch.out[,outcomes[1]], y=tmatch.out[,outcomes[2]], paired=TRUE)
	t2 <- t.test(x=tmatch.out[,outcomes[1]], y=tmatch.out[,outcomes[3]], paired=TRUE)
	t3 <- t.test(x=tmatch.out[,outcomes[2]], y=tmatch.out[,outcomes[3]], paired=TRUE)
	
	ci <- as.data.frame(rbind(t1$conf.int, t2$conf.int, t3$conf.int))
	ci$Treatments <- diffcols
	ci$estimate <- c(t1$estimate, t2$estimate, t3$estimate)

	out.box <- melt(tmatch.out[,c('id',diffcols)], id.vars='id')
	names(out.box) <- c('ID','Treatments','Difference')
	
	#Create new labels for each box
	labels <- gsub('-', ' - ', gsub('.out', '', diffcols))
	names(labels) <- diffcols
	
	p <- ggplot(out.box, aes(x=Treatments, y=Difference)) + 
		geom_boxplot() + 
		geom_crossbar(data=ci, aes(x=Treatments, ymin=V1, ymax=V2, y=estimate), 
					  color='green', fill='green', width=ci.width, alpha=.6) +
		geom_hline(color='blue', yintercept=0)
	if(plot.mean) {
		p <- p + 
			geom_point(data=ci, aes(x=Treatments, y=estimate), color='red', size=3) +
			geom_text(data=ci, aes(x=Treatments, y=estimate, 
					  label=prettyNum(estimate, digits=2)), vjust=-1)
	}
	
	p <- p + scale_x_discrete(NULL, labels=labels) + xlab(NULL)
	
	return(p)
}
