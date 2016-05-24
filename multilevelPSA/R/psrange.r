#' Estimates models with increasing number of comparision subjects starting from
#' 1:1 to using all available comparison group subjects.
#' 
#' @param df data frame with variables to pass to glm
#' @param treatvar vector representing treatment placement. Should be coded as
#'        0s (for control) and 1s (for treatment).
#' @param formula formula for logistic regression model
#' @param nsteps number of steps to estimate from 1:1 to using all control records.
#' @param nboot number of models to execute for each step.
#' @param samples the sample sizes to draw from control group for each step.
#' @param type either logistic for Logistic regression (using \code{glm} 
#'        function) or ctree for Conditional Inference Trees (using the 
#'        \code{ctree} function).
#' @param ... other parameters passed to glm.
#' @return a class of psrange that contains a summary data frame, a details data
#'         frame, and a list of each individual result from glm.
#' @export
psrange <- function(df, treatvar, formula, nsteps=10, nboot=10, samples, 
					type=c('logistic','ctree'), ... ) {
	if(missing(samples)) {
		samples <- seq(length(which(treatvar==1)), length(which(treatvar==0)),
					   (length(which(treatvar==0)) - length(which(treatvar==1))) / nsteps )
	}
	
	results <- list()
	
	ncontrol <- length(which(treatvar == 0))
	ntreat <- length(which(treatvar == 1))
	dfrange <- data.frame(p=integer(), i=integer(),
						 ntreat=integer(), ncontrol=integer(), 
						 psmin=numeric(), psmax=numeric())
	pb <- txtProgressBar(min=1, max=(length(samples)*nboot), style=3)
	fittedValues <- list()
	models <- list()
	for(i in 1:length(samples)) {
		tosample <- samples[i]
		models[[i]] <- list()
		density.df <- data.frame(iter=integer(), treat=integer(), ps=numeric())
		for(j in 1:nboot) {
			rows <- c(which(treatvar == 1),
					 sample(which(treatvar == 0), tosample))
			model.fit <- NA
			ps <- NA
			if(type[1] == 'logistic') {
				model.fit <- glm(formula, data=df[rows,], family='binomial', ...)
				ps <- fitted(model.fit)
			} else if(type[1] == 'ctree') {
				model.fit <- ctree(formula, data=df[rows,])
				ps <- sapply(treeresponse(model.fit), function(x) { x[[1]] })
			} else {
				stop('Unsupported type!')
			}
			dfrange <- rbind(dfrange, data.frame(ind=i, p=tosample/ncontrol*100, i=j,
												ntreat=ntreat, ncontrol=tosample, 
												psmin=range(ps)[1],
												psmax=range(ps)[2]))
			density.df <- rbind(density.df, 
					data.frame(iter=j, treat=treatvar[rows], ps=ps))
			models[[i]][[j]] <- model.fit
			setTxtProgressBar(pb, (((i-1)*nboot) + j))
		}
		fittedValues[[i]] <- density.df
	}
	
	results$fittedValues <- fittedValues
	results$models <- models
	dfrange$ratio <- dfrange$ncontrol / dfrange$ntreat
	results$details <- dfrange
	smin <- describeBy(dfrange$psmin, group=dfrange$p, mat=TRUE)[,
							c('mean','sd','median','se','min','max')]
	names(smin) <- paste('min', names(smin), sep='.')
	smax <- describeBy(dfrange$psmax, group=dfrange$p, mat=TRUE)[,
							c('mean','sd','median','se','min','max')]
	names(smax) <- paste('max', names(smax), sep='.')
	results$summary <- cbind(dfrange[!duplicated(dfrange$p),
					c('p','ntreat','ncontrol','ratio')], smin, smax)
	class(results) <- c('psrange')
	return(results)
}

#' Prints the summary results of psrange.
#' 
#' @param object psrange to print summary of.
#' @param ... currently unused.
#' @export
#' @method summary psrange
summary.psrange <- function(object, ...) {
	return(object$summary)
}

#' Prints information about a psrange result.
#' 
#' @param x psrange to print info about.
#' @param ... currently unused
#' @export
#' @method print psrange
print.psrange <- function(x, ...) {
	cat(paste('Executed ', length(x$models[[1]]), ' for ', length(x$models), 
		  ' ranging from 1:', min(x$summary$ratio), ' to 1:', max(x$summary$ratio),
		  '\n', sep=''))
	print(x$summary[,c('ntreat','ncontrol','ratio','p','min.mean','max.mean')])
}

utils::globalVariables(c('ps', '..scaled..','treat','max.min','max.max','min.min','min.max',
						 'psmin','psmax','min.mean','max.mean','ncontrol','ratio'))

#' Plots densities and ranges for the propensity scores.
#' 
#' @param x the result of psrange.
#' @param xlab label for x-axis.
#' @param ylab label for y-axis.
#' @param labels labels for the comparison and treatment legend.
#' @param text.ratio.size size of the text for the ratio.
#' @param text.ncontrol.size size of the text for the number of control units.
#' @param point.size size of the points for the minimum and maximum ranges for
#'        each model.
#' @param point.alpha the alpha (transparency) level for the points.
#' @param line.width the width of the line between the median of the minimum
#'        and maximum ranges.
#' @param density.alpha the alpha (transparency) level of the density curves.
#' @param rect.color the color of the rectangle surrounding the range of minimum
#'        and maximum ranges.
#' @param rect.alpha the alpha (transparency) level of the rectangle.
#' @param ... currently unused.
#' @return a ggplot2 object
#' @method plot psrange
#' @export
plot.psrange <- function(x,
		 xlab=NULL,
		 ylab=NULL,
		 labels=c('Comparison','Treatment'),
		 text.ratio.size = 4,
		 text.ncontrol.size = 3,
		 point.size = 1, 
		 point.alpha = .6,
		 line.width = 6,
		 density.alpha = .2,
		 rect.color = 'green',
		 rect.alpha = .2,
		 ...
) {
	if(missing(xlab)) {
		xlab <- paste('Propensity Score Range (ntreat = ', 
				   prettyNum(x$summary[1,'ntreat'], big.mark=','), ')', sep='')
	}
	
	densities.df <- data.frame(p=numeric(), treat=integer(), ps=numeric())
	for(i in seq_len(length(x$fittedValues))) {
		densities.df <- rbind(densities.df, cbind(p=x$summary[i,'p'], 
												  x$fittedValues[[i]]))
	}
	#densities.df$treat = as.integer(densities.df$treat)
	if(is.logical(densities.df$treat)) {
		densities.df$treat = factor(as.integer(densities.df$treat))
	} else {
		densities.df$treat = factor(densities.df$treat)
	}
	
	text.vjust = -.4
	bar.factor = 1

	psrange_labeller <- function(variable, value) {
		if(variable == 'p') {
			return(paste(round(as.numeric(as.character(value))), '%', sep=''))
		} else {
			return(value)
		}
	}
	
	p <- ggplot() + 
			xlim(c(-.05,1.05)) + ylim(c(-1,1)) +
			stat_density(data=densities.df[densities.df$treat==1,], 
				aes(x=ps, ymax=-..scaled.., fill=treat, ymin = 0),
				geom="ribbon", position="identity", alpha=density.alpha) +
			stat_density(data=densities.df[densities.df$treat==0,], 
				aes(x=ps, ymax=..scaled.., fill=treat, ymin = 0),
				geom="ribbon", position="identity", alpha=density.alpha) +
			geom_rect(data=x$summary, aes(group=p, xmin=(max.min-.005), xmax=(max.max+.005)),
				ymin=0.25, ymax=.75, fill=rect.color, alpha=rect.alpha) +
			geom_rect(data=x$summary, aes(group=p, xmin=(min.min-.005), xmax=(min.max+.005)), 
				ymin=-.75, ymax=-0.25, fill=rect.color, alpha=rect.alpha) +
			geom_point(data=x$details, aes(y=-.5, x=psmin), 
				size=point.size, alpha=point.alpha, shape=23) +
  		  	geom_point(data=x$details, aes(y=.5, x=psmax), 
  		  		size=point.size, alpha=point.alpha, shape=22) +
			geom_errorbarh(data=x$summary, 
				aes(y=0, x=min.mean + (max.mean-min.mean)/ 2, xmin=min.mean, xmax=max.mean), 
				colour='black', width=line.width) + 
		   	geom_text(data=x$summary,
		   		aes(label=paste(prettyNum(floor(ncontrol), big.mark=','), sep='')), 
		   		x=0, y=0, size=text.ncontrol.size, hjust=1.1, vjust=-0.2) +
   		  	geom_text(data=x$summary,
   		  		aes(label=paste('1:', round(ratio, digits=1), sep=''), 
   		  		x=(min.mean + (max.mean-min.mean)/2)), y=0, 
   		  		size=text.ratio.size, vjust=text.vjust) +
		  	facet_grid(p ~ ., as.table=FALSE, labeller=psrange_labeller) +
			theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
				  strip.background=element_rect(fill='#EFEFEF'),
				  strip.text.y=element_text(angle=360)) +
			ylab(ylab) + xlab(xlab) +
			scale_fill_hue('', limits=c(0,1), labels=labels)
	
	return(p)
}

