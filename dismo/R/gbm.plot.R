# function to plot gbm response variables, with the option
# of adding a smooth representation of the response if requested
# additional options in this version allow for plotting on a common scale
# note too that fitted functions are now centered by subtracting their mean
# 
# version 2.9
#
# j. leathwick/j. elith - March 2007
#


gbm.plot <-
function(gbm.object,                # a gbm object - could be one from gbm.step
     variable.no = 0,               # the var to plot - if zero then plots all
     smooth = FALSE,                # should we add a smoothed version of the fitted function 
     rug = TRUE,                    # plot a rug of deciles
     n.plots = length(pred.names),  # plot the first n most important preds
     common.scale = TRUE,           # use a common scale on the y axis
     write.title = TRUE,            # plot a title above the plot
     y.label = "fitted function",   # the default y-axis label
     x.label = NULL,                # the default x-axis label
     show.contrib = TRUE,           # show the contribution on the x axis
     plot.layout = c(3,4),          # define the default layout for graphs on the page
     ...                            # other arguments to pass to the plotting 
                                    # useful options include cex.axis, cex.lab, etc.
     )
{

	if (! requireNamespace('gbm') ) { stop ('you need to install the gbm package to run this function') }
	requireNamespace('splines')
	
	gbm.call <- gbm.object$gbm.call
	gbm.x <- gbm.call$gbm.x
	pred.names <- gbm.call$predictor.names
	response.name <- gbm.call$response.name
	
	data <- gbm.call$dataframe

	max.plots <- plot.layout[1] * plot.layout[2]
	plot.count <- 0
	n.pages <- 1

	if (length(variable.no) > 1) { stop("only one response variable can be plotted at a time") }

	if (variable.no > 0) {   #we are plotting all vars in rank order of contribution
		n.plots <- 1
	}

	max.vars <- length(gbm.object$contributions$var)
	if (n.plots > max.vars) {
		n.plots <- max.vars
		warning("reducing no of plotted predictors to maximum available (",max.vars,")")
	}

	predictors <- list(rep(NA,n.plots)) # matrix(0,ncol=n.plots,nrow=100)
	responses <- list(rep(NA,n.plots)) # matrix(0,ncol=n.plots,nrow=100)

	for (j in c(1:n.plots)) {  #cycle through the first time and get the range of the functions
		if (n.plots == 1) {
			k <- variable.no
		} else {
			k <- match(gbm.object$contributions$var[j],pred.names)
		}
		if (is.null(x.label)) {
			var.name <- gbm.call$predictor.names[k]
		} else {
			var.name <- x.label
		}
		pred.data <- data[ , gbm.call$gbm.x[k]]

		response.matrix <- gbm::plot.gbm(gbm.object, k, return.grid = TRUE)

		predictors[[j]] <- response.matrix[,1]
		if (is.factor(data[,gbm.call$gbm.x[k]])) {
			predictors[[j]] <- factor(predictors[[j]],levels = levels(data[,gbm.call$gbm.x[k]]))
		}
		responses[[j]] <- response.matrix[,2] - mean(response.matrix[,2])

		if(j == 1) {
			ymin = min(responses[[j]])
			ymax = max(responses[[j]])
		} else {
			ymin = min(ymin,min(responses[[j]]))
			ymax = max(ymax,max(responses[[j]]))
		}
	}

# now do the actual plots

	op <- graphics::par(no.readonly = TRUE) 
	graphics::par(mfrow = plot.layout)

	for (j in c(1:n.plots)) {

		if (plot.count == max.plots) {
			plot.count = 0
			n.pages <- n.pages + 1
		}

		plot.count <- plot.count + 1

		if (n.plots == 1) {
			k <- match(pred.names[variable.no],gbm.object$contributions$var)
			if (show.contrib) {
				x.label <- paste(var.name,"  (",round(gbm.object$contributions[k,2],1),"%)",sep="")
			}
		} else {
			k <- match(gbm.object$contributions$var[j],pred.names)
			var.name <- gbm.call$predictor.names[k]
			if (show.contrib) {
				x.label <- paste(var.name,"  (",round(gbm.object$contributions[j,2],1),"%)",sep="")
			} else x.label <- var.name
		}

		if (common.scale) {
			plot(predictors[[j]],responses[[j]],ylim=c(ymin,ymax), type='l',
			xlab = x.label, ylab = y.label, ...)
		} else {
			plot(predictors[[j]],responses[[j]], type='l', 
			xlab = x.label, ylab = y.label, ...)
		}
		if (smooth & is.vector(predictors[[j]])) {
			temp.lo <- loess(responses[[j]] ~ predictors[[j]], span = 0.3)
			lines(predictors[[j]],fitted(temp.lo), lty = 2, col = 2)
		}
		if (plot.count == 1) {
			if (write.title) {
				title(paste(response.name," - page ",n.pages,sep=""))
			}
			if (rug & is.vector(data[,gbm.call$gbm.x[variable.no]])) {
				rug(quantile(data[,gbm.call$gbm.x[variable.no]], probs = seq(0, 1, 0.1), na.rm = TRUE))
			}
		} else {
			if (write.title & j == 1) {
				title(response.name)
			}
			if (rug & is.vector(data[,gbm.call$gbm.x[k]])) {
				rug(quantile(data[,gbm.call$gbm.x[k]], probs = seq(0, 1, 0.1), na.rm = TRUE))
			}
		}
	}
	graphics::par(op)

}

