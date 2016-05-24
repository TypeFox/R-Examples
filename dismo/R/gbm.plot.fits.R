# j leathwick, j elith - 7th January 2005
#
# version 2.0 - developed in R 2.0
#
# to plot distribution of fitted values in relation to ydat from mars or other p/a models
# allows masking out of absences to enable focus on sites with high predicted values
# fitted values = those from model; raw.values = original y values
# label = text species name; ydat = predictor dataset
# mask.presence forces function to only plot fitted values for presences
# use.factor forces to use quicker printing box and whisker plot
# file.name routes to a pdf file of this name
#

gbm.plot.fits <- function(gbm.object,  
	v = 0,
	mask.presence = FALSE, 
	use.factor = FALSE )
{


    gbm.call <- gbm.object$gbm.call	#and the call details
    gbm.x <- gbm.call$gbm.x    
    gbm.y <- gbm.call$gbm.y
    family <- gbm.call$family
    dat <- gbm.call$dataframe
    n.cases <- nrow(dat)

    xdat <- as.data.frame(dat[,gbm.x])
    ydat <- as.data.frame(dat[,gbm.y])

    n.preds <- ncol(xdat)

	if (v==0) {
		v <- 1:n.preds
		if (length(v) > 16) {
			warning('only first 16 layers are plotted')
			v <- 1:16
		}
	}
	nl <- length(v)
	nc <- ceiling(sqrt(nl))
	nr <- ceiling(nl / nc)
	old.par <- graphics::par(no.readonly = TRUE) 
	on.exit(graphics::par(old.par))
	graphics::par(mfrow=c(nr, nc))

    fitted.values <- gbm.object$fitted

    pred.names <- names(dat)[gbm.x]
    sp.name <- names(dat)[gbm.y]

    if (mask.presence) {
		mask <- ydat == 1 
	} else {
		mask <- rep(TRUE, length = n.cases) 
	}

    robust.max.fit <- approx(ppoints(fitted.values[mask]), sort(fitted.values[mask]), 0.99) #find 99%ile value
	
    for (j in v) {

		if (is.numeric(xdat[mask,j])) {
			wt.mean <- zapsmall(mean((xdat[mask, j] * fitted.values[mask]^5)/mean(fitted.values[mask]^5),na.rm=TRUE),2)
		} else {
			wt.mean <- "na"
		}
		
		if (use.factor) {
			temp <- factor(cut(xdat[mask, j], breaks = 12))
			if (family == "binomial") {
				plot(temp, fitted.values[mask], xlab = pred.names[j], ylab = "fitted values", ylim = c(0, 1))
			} else {
				plot(temp, fitted.values[mask], xlab = pred.names[j], ylab = "fitted values")
			}
		} else {
			if (family == "binomial") {
				plot(xdat[mask, j], fitted.values[mask], xlab = pred.names[j], ylab = "fitted values", ylim = c(0, 1))
			} else {
				plot(xdat[mask, j], fitted.values[mask], xlab = pred.names[j], ylab = "fitted values")}
		    }
			abline(h = (0.333 * robust.max.fit$y), lty = 2.)
		if (j == 1) { 
			title(paste(sp.name, ", wtm = ", wt.mean))
		} else {
			title(paste("wtm = ", wt.mean))
		}
	}
}


