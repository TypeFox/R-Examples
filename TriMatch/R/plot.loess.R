utils::globalVariables(c('ps','treat','id'))

#' Loess plot for matched triplets.
#' 
#' This function will create a \code{ggplot2} figure with propensity scores on the
#' x-axis and the outcome on the y-axis. Three Loess regression lines will be plotted
#' based upon the propensity scores from \code{model}. Since each model produces
#' propensity scores for two of the three groups, the propensity score for the third
#' group in each matched triplet will be the mean of the other two. If \code{model}
#' is not specified, the default will be to use the model that estimates the propensity
#' scores for the first two groups in the matching order.
#' 
#' @param tmatch the results of \code{\link{trimatch}}.
#' @param outcome a vector representing the outcomes.
#' @param model an integer between 1 and 3 indicating from which model the propensity
#'        scores will be used.
#' @param ylab the label for the y-axis.
#' @param plot.connections boolean indicating whether lines will be drawn connecting
#'        each matched triplet.
#' @param connections.color the line color of connections.
#' @param connections.alpha number between 0 and 1 representing the alpha levels for 
#'        connection lines.
#' @param plot.points a \code{ggplot2} function for plotting points. Usually 
#'        \code{\link{geom_point}} or \code{\link{geom_jitter}}. If \code{NULL} no points
#'        will be drawn.
#' @param points.alpha number between 0 and 1 representing the alpha level for the points.
#' @param points.palette the color palette to use. See \code{\link{scale_colour_brewer}}
#'        and \url{http://colorbrewer2.org/} for more information.
#' @param ... other parameters passed to \code{\link{geom_smooth}} and 
#'        \code{\link{stat_smooth}}.
#' @return a \code{ggplot2} figure.
#' @export
loess3.plot <- function(tmatch, outcome, model,
						ylab='Outcome', 
						plot.connections=FALSE,
						connections.color='black',
						connections.alpha=.2,
						plot.points=geom_point,
						points.alpha=.1,
						points.palette='Dark2', ...) {
	tpsa <- attr(tmatch, 'triangle.psa')
	tmatch2 <- merge(x=tmatch, y=outcome)
	groups <- names(tmatch2)[1:3]
	
	if(missing(model)) {
		for(i in 1:3) {
			if(length(which(is.na(tpsa[tpsa$treat %in% groups[1:2],
									   paste('model', i, sep='')]))) == 0) {
				model <- i
				break;
			}
		}
		if(model == 0) {
			stop('Could not find model. There are missing propensity scores in all models.')
		}
	} else {
		groups <- c(
			as.character(tpsa[which(!tpsa[,paste0('model',model)]),'treat'][1]),
			as.character(tpsa[which(tpsa[,paste0('model',model)]),'treat'][1]),
			as.character(tpsa[which(is.na(tpsa[,paste0('model',model)])),'treat'][1])
		)
	}

	xlab <- paste0('Propensity Score (0=', 
				   tpsa[which(!tpsa[,paste0('model',model)]),'treat'][1],
				   ', 1=',
				   tpsa[which(tpsa[,paste0('model',model)]),'treat'][1],
				   ')')
	
	tmatch2 <- merge(tmatch2, tpsa[which(tpsa$treat == groups[1]), 
								   c('id',paste('ps', model, sep=''))], 
					 by.x=groups[1], by.y='id', all.x=TRUE)
	names(tmatch2)[ncol(tmatch2)] <- paste(groups[1], '.ps', sep='')
	tmatch2 <- merge(tmatch2, tpsa[which(tpsa$treat == groups[2]), 
								   c('id',paste('ps', model, sep=''))], 
					 by.x=groups[2], by.y='id', all.x=TRUE)
	names(tmatch2)[ncol(tmatch2)] <- paste(groups[2], '.ps', sep='')
	tmatch2[,paste(groups[3], '.ps', sep='')] <- apply(
			tmatch2[,(ncol(tmatch2)-1):ncol(tmatch2)], 1, mean)
	
	tmatch2$id <- 1:nrow(tmatch2)
	
	out <- cbind(melt(tmatch2[,c(paste(groups, '.out', sep=''), 'id')], id.vars='id'),
	             melt(tmatch2[,c(paste(groups, '.ps', sep=''), 'id')], id.vars='id') )[,c(1,2,3,6)]
	names(out) <- c('id','treat','out','ps')
	out$treat <- as.character(out$treat)
	out$treat <- substr(out$treat, 1, (sapply(out$treat, nchar)-4)) #Strip .out from value
	
	p <- ggplot(out, aes(x=ps, y=out, group=treat, fill=treat, colour=treat))
	if(plot.connections) {
		p <- p + geom_path(aes(group=id), colour=connections.color, alpha=connections.alpha)
	}
	if(!is.null(plot.points)) {
		p <- p + plot.points(alpha=points.alpha)
	}
	p <- p + geom_smooth(...)
	p <- p + ylab(ylab) + xlab(xlab)
	p <- p + scale_color_brewer('Treatment', palette=points.palette)
	p <- p + scale_fill_brewer('Treatment', palette=points.palette)
	p
}
