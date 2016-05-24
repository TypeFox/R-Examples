utils::globalVariables(c('Treatment','Covariate','Mean','Strata','ymin','ymax',
						 '..count..','ID','describeBy','cstrata.psa','cv.trans.psa'))

#' Balance plot for the given covariate.
#'
#' If the covariate is numeric, boxplots will be drawn with red points for the mean
#' and green error bars for the standard error. For non-numeric covariates a barplot
#' will be drawn.
#' 
#' A Friedman rank sum test will be performed for all covariate types, printed,
#' and stored as an attribute to the returned object named \code{friedman}. If
#' a continuous covariate a repeated measures ANOVA will also be performed, printed,
#' and returned as an attribute named \code{rmanova}.
#' 
#' @param x results from \code{\link{trimatch}}.
#' @param covar vector of the covaraite to check balance of.
#' @param model an integer between 1 and 3 indicating from which model the 
#'        propensity scores will be used.
#' @param nstrata number of strata to use.
#' @param ylab label of the y-axis.
#' @param xlab label of the x-axis.
#' @param se.ratio a multiplier for how large standard error bars will be.
#' @param label label for the legend.
#' @param print print the output if the Freidman Rank Sum Test and repeated
#'        measures ANOVA (for continuous variables).
#' @param legend.position the position of the legend. See \code{\link{theme}}.
#' @param x.axis.labels labels for the x-axis.
#' @param x.axis.angle angle for x-axis labels.
#' @param ... parameters passed to \code{\link{plot.balance.plots}}.
#' @return a \code{ggplot2} figure or a list of \code{ggplot2} figures if \code{covar}
#'        is a data frame.
#' @export
balance.plot <- function(x, covar, model,
					     nstrata=attr(attr(tmatch, 'triangle.psa'), 'nstrata'),
					     label='Covariate',
					     ylab='',
						 xlab=NULL,
						 se.ratio = 2,
						 print=TRUE,
						 legend.position='top',
						 x.axis.labels,
						 x.axis.angle = -45,
						 ...) {
	tmatch <- x
	if(is.data.frame(covar)) {
		#If covar is a data frame, create a grid of figures.
		plots <- list()
		for(i in seq_along(covar)) {
			psub <- balance.plot(tmatch,
								 covar[,i], 
								 model = model,
								 nstrata = nstrata,
								 label = names(covar)[i],
								 ylab = names(covar)[i],
								 se.ratio = se.ratio,
								 legend.position = legend.position,
								 x.axis.labels = x.axis.labels,
								 x.axis.angle = x.axis.angle,
								 print = FALSE, ...)
			plots[[names(covar)[i]]] <- psub
		}
		class(plots) <- 'balance.plots'
		return(plots)
	} else {	
		if(!is.numeric(covar)) {
			covar <- as.character(covar)
		}
		
		#TODO: Much of this code is shared with plot.loess3
		tpsa <- attr(tmatch, 'triangle.psa')
		tmatch2 <- merge(x=tmatch, y=covar)
		groups <- names(tmatch2)[1:3]
		
		if(missing(model)) {
			for(i in 1:3) {
				if(length(which(is.na(tpsa[tpsa$treat %in% 
								groups[1:2],paste('model', i, sep='')]))) == 0) {
					model <- i
					break;
				}
			}
			if(model == 0) {
				stop('Could not find model. There are missing propensity scores in all models.')
			}
			if(print) {
				message(paste('Using propensity scores from model ', model, 
							  ' for evaluating balance.', sep=''))
			}
		} else {
			#Need to determine which groups these models represent
			if(length(which(is.na(tpsa[tpsa$treat %in% groups[1:2],
									paste('model', model, sep='')]))) == 0) {
				groups <- groups[c(1,2,3)]
			} else if(length(which(is.na(tpsa[tpsa$treat %in% groups[2:3],
									paste('model', model, sep='')]))) == 0) {
				groups <- groups[c(2,3,1)]
				
			} else if(length(which(is.na(tpsa[tpsa$treat %in% groups[c(1,3)],
									paste('model', model, sep='')]))) == 0) {
				groups <- groups[c(1,3,2)]			
			} else {
				stop(paste0('Could not use model ', model, 
						   '. There is are unexpected missing propensity scores'))
			}
		}
	
		tmatch2 <- merge(tmatch2, tpsa[which(tpsa$treat == groups[1]), 
									   c('id',paste('ps', model, sep=''))], 
						 by.x=groups[1], by.y='id', all.x=TRUE)
		names(tmatch2)[ncol(tmatch2)] <- paste(groups[1], '.ps', sep='')
		tmatch2 <- merge(tmatch2, tpsa[which(tpsa$treat == groups[2]), 
									   c('id',paste('ps', model, sep=''))], 
						 by.x=groups[2], by.y='id', all.x=TRUE)
		names(tmatch2)[ncol(tmatch2)] <- paste(groups[2], '.ps', sep='')
		tmatch2[,'mean.ps'] <- apply(tmatch2[,(ncol(tmatch2)-1):ncol(tmatch2)], 1, mean)
	
		breaks <- quantile(tmatch2$mean.ps, probs=seq(0,1,1/nstrata), na.rm=TRUE)
		tmatch2$strata <- cut(tmatch2$mean.ps, 
							  breaks=breaks, 
							  labels=1:nstrata, 
							  include.lowest=TRUE)
		
		tmatch2$id <- 1:nrow(tmatch2)
		
		badrows <- which(is.na(tmatch2$strata))
		if(length(badrows) > 0) {
			warning(paste('Could not determine strata for the following rows: ', 
						  paste(badrows, collapse=', '), sep=''))
			tmatch2 <- tmatch2[-badrows,] 
		}
		
		out <- melt(tmatch2[,c(paste(groups, '.out', sep=''), 'id', 'strata')], 
					id.vars=c('id','strata'))
		names(out) <- c('ID','Strata','Treatment','Covariate')	
		out$Treatment <- as.character(out$Treatment)
		#Strip .out from value
		out$Treatment <- substr(out$Treatment, 1, (sapply(out$Treatment, nchar)-4))
		
		p <- ggplot(out, aes(x=Treatment)) + facet_grid(~ Strata) 
		if(is.numeric(out$Covariate)) {
			df <- describeBy(out$Covariate, group=list(out$Treatment, out$Strata), 
							 mat=TRUE, skew=FALSE)[,c('group1','group2','mean','sd','se')]
			names(df) <- c('Treatment','Strata','Mean','SD','SE')
			df$ymin <- (df$Mean - se.ratio * df$SE)
			df$ymax <- (df$Mean + se.ratio * df$SE)
			p <- p + geom_boxplot(aes(y=Covariate)) +
				geom_point(data=df, aes(y=Mean), color='red', size=3) +
				geom_line(data=df, aes(y=Mean, group=Strata)) +
				geom_errorbar(data=df, aes(ymin=ymin, ymax=ymax), color='green', width=.5)
		} else { #Categorical varaible
			p <- p + geom_bar(aes(fill=factor(Covariate), y=(..count..)/sum(..count..)), 
							  position='fill')  +
				scale_fill_hue(label) + ylab('Percent') +
				scale_y_continuous(labels = percent) +
				theme(legend.position=legend.position)
		}
		
		p <- p + theme(axis.text.x=element_text(angle=x.axis.angle, vjust=.5),
					   panel.background=element_rect(color='black', fill='#F9F3FD'))
		p <- p + xlab(xlab) + ylab(ylab)
		
		if(!missing(x.axis.labels)) {
			p <- p + scale_x_discrete(labels=c('C','T1','T2'))
		}
		
		ft <- friedman.test(Covariate ~ Treatment | ID, out)
		if(print) { print(ft) }
		attr(p, 'friedman') <- ft
		
		if(is.numeric(out$Covariate)) {
			#We can also use repeated measure ANOVA for continuous covariates
			out$Treatment <- as.factor(out$Treatment)
			out$ID <- as.factor(out$ID)
			rmanova <- ezANOVA(data=out, dv=Covariate, wid=ID, within=Treatment)
			if(print) {
				cat(' Repeated measures ANOVA\n\n')
				print(rmanova$ANOVA)
			}
			attr(p, 'rmanova') <- rmanova
		}
		
		return(p)
	}
}

#' Print the results of \code{\link{balance.plot}} for a data frame of covariates.
#' 
#' @param x the results of \code{\link{balance.plot}} when a data frame is specified.
#' @param ... parameters passed to \code{\link{plot.balance.plots}} and
#'        \code{\link{summary.balance.plots}}.
#' @method print balance.plots
#' @export
print.balance.plots <- function(x, ...) {
	print(summary(x, ...))
	plot(x, ...)
}

#' Prints a summary table of the test statistics of each balance plot.
#' 
#' The \code{\link{balance.plot}} function will create a grid of balance plots
#' if a data frame is provided. The returned object is a list of \code{ggplot2}
#' figures with the statistical tests (i.e. Friedmen Rank Sum tests and if a
#' continuous variable, repeated measures ANOVA as well) saved as attributes.
#' This function will return a data frame combining all of those results.
#' 
#' @param object the results of \code{\link{balance.plot}} when a data frame is specified.
#' @param ... currenlty unused.
#' @return a data frame
#' @export
#' @method summary balance.plots
summary.balance.plots <- function(object, ...) {
	friedman.results <- data.frame(Covariate=character(), 
								   Friedman=numeric(), Friedman.p=numeric(),
								   rmANOVA=numeric(), rmANOVA.p=numeric(),
								   stringsAsFactors=FALSE)
	for(i in seq_along(object)) {
		f <- attr(object[[i]], 'friedman')
		a <- attr(object[[i]], 'rmanova')
		friedman.results <- rbind(friedman.results, data.frame(
			Covariate = names(object)[i],
			Friedman = f$statistic,
			Friedman.p = f$p.value,
			Friedman.sig = star(f$p.value),
			rmANOVA = ifelse(is.null(a), NA, a[[1]]$F),
			rmANOVA.p = ifelse(is.null(a), NA, a[[1]]$p),
			rmANOVA.sig = star(ifelse(is.null(a), NA, a[[1]]$p)),
			stringsAsFactors=FALSE
		))
	}
	row.names(friedman.results) <- 1:nrow(friedman.results)
	return(friedman.results)
}

#' Prints a grid of balance plots.
#' 
#' @param x the results of \code{\link{balance.plot}} when a data frame is specified.
#' @param rows if \code{covar} is a data frame of covariates, the number of
#'        rows in the grid of figures.
#' @param cols if \code{covar} is a data farme of covariates, the number of
#'        columns in the grid of figures.
#' @param byrow if TRUE (default), plots will be drawn by rows, otherwise by columns.
#' @param plot.sequence the sequence (or subset) of plots to draw.
#' @param ... currenlty unused.
#' @method plot balance.plots
#' @export
plot.balance.plots <- function(x, rows, cols, byrow = TRUE, 
							   plot.sequence=seq_along(bplots), ...) {
	bplots <- x
	grid.newpage()
	if(missing(rows) & missing(cols)) {
		rows <- ceiling(sqrt(length(bplots)))
		cols <- ceiling(length(bplots) / rows)
	} else if(missing(rows)) {
		rows <- ceiling(length(bplots) / cols)
	} else if(missing(cols)) {
		cols <- ceiling(length(bplots) / rows)
	}
	pushViewport(viewport(layout = grid.layout(rows, cols)))
	vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
	for(i in seq_along(plot.sequence)) {
		psub <- bplots[[plot.sequence[i]]]
		if(byrow) {
			print(psub, vp=vplayout( ceiling(i / cols) , (i-1) %% cols + 1 ))
		} else {
			print(psub, vp=vplayout( (i-1) %% rows + 1 , ceiling(i / rows) ))
		}
	}
}

