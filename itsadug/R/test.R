#' Function for comparing two GAMM models.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @param model1 First model.
#' @param model2 Second model.
#' @param signif.stars Logical (default = TRUE). Whether or not to display 
#' stars indicating the level of significance on 95\% confidence level.
#' @param print.output Logical: whether or not to print the output. 
#' By default controlled globally in the package options:  
#' If the function \code{\link{infoMessages}} is set to TRUE, the output 
#' will be automatically printed.
#' Could be also set by explicitly providing TRUE or FALSE. See notes.
#' @details
#' 
#' As an Chi-Square test is performed on two times the difference in 
#' minimized smoothing parameter selection score (GCV, fREML, REML, ML), 
#' and the difference in degrees of freedom specified in the model. 
#' The degrees of freedom of the model terms are the sum of
#' 1) the number of estimated smoothing parameters for the model, 
#' 2) number of parametric (non-smooth) model terms including the intercept, 
#' and 3) the sum of the penalty null space dimensions of each smooth object.
#'
#' This method is preferred over other functions such as \code{\link{AIC}} for 
#' models that include an AR1 model or random effects (especially nonlinear 
#' random smooths using \code{bs="fs"}). CompareML also reports the AIC 
#' difference, but that value should be treated with care.
#' 
#' Note that the Chi-Square test will result in a very low p-value
#' when the difference in degrees of freedom approaches zero. Use common sense 
#' to determine if the difference between the two models is meaningful. 
#' A warning is presented when the difference in score is smaller 
#' than 5.
#'
#' The order of the two models is not important.
#' Model comparison is only implemented for the methods GCV, fREML, REML, and ML.
#'
#' @section Notes:
#' If no output is provided in the command window, set info messages to TRUE: 
#' \code{infoMessages("on")} and try again.
#' For suppressing the output and all warnings, set infoMessages to FALSE 
#' (\code{infoMessages("on")} ) and use the function 
#' \code{\link{suppressWarnings}} to suppress warning messages.
#' @return Optionally returns the Chi-Square test table.
#' @author Jacolien van Rij. With many thanks to Simon N. Wood for his feedback.
#' @seealso For models without AR1 model or random effects \code{\link{AIC}} can be used.
# help function
#' @examples
#' data(simdat)
#'
#' \dontrun{
#' infoMessages("on")
#' # some arbitrary models:
#' m1 <- bam(Y~Group + s(Time, by=Group), method="fREML", data=simdat)
#' m2 <- bam(Y~Group + s(Time), method="fREML", data=simdat)
#' 
#' compareML(m1, m2)
#'
#' # exclude significance stars:
#' compareML(m1, m2, signif.stars=FALSE)
#'
#' m3 <- bam(Y~Group + s(Time, by=Group, k=25), method="fREML", 
#'     data=simdat)
#' compareML(m1, m3)
#' 
#' # do not print output, but save table for later use:
#' cml <- compareML(m1, m2, print.output=FALSE)$table
#' cml
#' # alternative way:
#' infoMessages('off')
#' cml <- compareML(m1, m2, print.output=FALSE)$table
#' infoMessages('on')
#'
#' # Use suppressWarnings to also suppress warnings:
#' suppressWarnings(cml <- compareML(m1, m2, print.output=FALSE)$table)
#' 
#' }
#' @family Testing for significance
compareML <- function(model1, model2,
    signif.stars=TRUE,
    print.output=getOption('itsadug_print')) {
    # check gam or bam model:
    if((!"gam" %in% class(model1)) | (!"gam" %in% class(model2))){
        stop("Models are not gam objects (i.e., build with bam()/gam()).")
    }
    
    # check whether models are comparable:
    if (model1$method != model2$method) {
        stop(sprintf("Models are incomparable: method model1 = %s, method model2 = %s", model1$method, model2$method))
    }
    
    
    type <- model1$method
    
    ml1 <- model1$gcv.ubre[1]
    ml2 <- model2$gcv.ubre[1]
    
    ### OLD METHOD, SIMON SAYS NOT OK! ### edf1 <- sum(model1$edf) edf2 <- sum(model2$edf) NEW METHOD: ###
    ndf1 <- length(model1$sp) + model1$nsdf + ifelse(length(model1$smooth)>0,
        sum(sapply(model1$smooth, FUN = function(x) {
            x$null.space.dim
            }, USE.NAMES = FALSE)), 0)
    ndf2 <- length(model2$sp) + model2$nsdf + ifelse(length(model2$smooth)>0,
        sum(sapply(model2$smooth, FUN = function(x) {
            x$null.space.dim
            }, USE.NAMES = FALSE)), 0)
    
    if (! model1$method %in%  c("fREML", "REML", "ML")) {
        type <- "AIC"
        ml1 <- AIC(model1)
        ml2 <- AIC(model2)
        
	    ndf1 <- length(model1$sp) + model1$nsdf + ifelse(length(model1$smooth)>0,
	        sum(sapply(model1$smooth, FUN = function(x) {
	            x$null.space.dim
	            }, USE.NAMES = FALSE)), 0)
	    ndf2 <- length(model2$sp) + model2$nsdf + ifelse(length(model2$smooth)>0,
	        sum(sapply(model2$smooth, FUN = function(x) {
	            x$null.space.dim
	            }, USE.NAMES = FALSE)), 0)
        warning(sprintf("\nCompareML is not implemented for smoothing parameter estimation method %s. AIC scores are used for model comparison. Consider running the model with REML, fREML, or ML as method.\n-----\n", model1$method))
    }
    
    
    # pchisq(4, .5, lower.tail=F) # p < .1 pchisq(-4, .5, lower.tail=F) # p = 1 pchisq(4, -.5, lower.tail=F) # NaN
    
    # Book keeping;
    info <- sprintf("%s: %s\n\n%s: %s\n", deparse(substitute(model1)), paste(deparse(model1$formula), collapse="\n"),
        deparse(substitute(model2)), paste(deparse(model2$formula), collapse="\n"))
    if(print.output){
        cat(info)
    }
    
    out      <- NULL
    advice   <- NULL
    warning  <- NULL
    
    # if (type != 'AIC') {
    # Situation 1: model 1 has lower score, but model 2 has lower df. Is it significantly better model than model 2?
	# Situation 0: equal df
	if(abs(round(ndf2 - ndf1)) < .5){
		if( ml1 < ml2){
            advice <- sprintf("\nModel %s preferred: lower %s score (%.3f), and equal df (%.3f).\n-----\n", 
            	deparse(substitute(model1)), 
                type, ml2 - ml1, ndf2 - ndf1)
            out <- data.frame(Model = c(deparse(substitute(model2)), deparse(substitute(model1))), 
            	Score = c(ml2, ml1), 
            	Edf = c(ndf2, ndf1), 
            	Difference = c("", sprintf("%.3f", ml2 - ml1)), 
            	Df = c("", sprintf("%.3f", ndf1 - ndf2)),
                p.value = c("",NA),
                Sign. = c("", ""))
		}else{
            advice <- sprintf("\nModel %s preferred: lower %s score (%.3f), and equal df (%.3f).\n-----\n", 
            	deparse(substitute(model2)), 
                type, ml1 - ml2, ndf1 - ndf2)
            out <- data.frame(Model = c(deparse(substitute(model1)), deparse(substitute(model2))), 
            	Score = c(ml1, ml2), 
            	Edf = c(ndf1, ndf2),
            	Difference = c("", sprintf("%.3f", ml2 - ml1)), 
            	Df = c("", sprintf("%.3f", ndf1 - ndf2)),
                p.value = c("",NA),
                Sign. = c("", ""))
		}
	# Situation 1: model 1 has lower score, but model 2 has lower df. Is it significantly better model than model 2?
    }else if ((ml1 < ml2) & (ndf2 < ndf1)) {
        
        # twice the amount of difference in likelihood
        h1 <- pchisq(2 * (ml2 - ml1), abs(ndf1 - ndf2), lower.tail = F)
        
        out <- data.frame(Model = c(deparse(substitute(model2)), deparse(substitute(model1))), 
        	Score = c(ml2, ml1), 
        	Edf = c(ndf2, ndf1), 
            Chisq = c("", sprintf("%.3f", ml2 - ml1)), 
            Df = c("", sprintf("%.3f", abs(ndf1 - ndf2))), 
            p.value = c("", ifelse(h1 < 2e-16, sprintf(" < 2e-16"), 
            	ifelse(h1 < 0.001, sprintf("%.3e", h1), 
            	ifelse(h1 < 0.01, sprintf("%.3f", h1), 
            	ifelse(h1 < 0.05, sprintf("%.3f", h1), sprintf("%.3f", h1)))))), 
            Sig. = c("", ifelse(h1 < 0.001, sprintf("***", h1), 
            	ifelse(h1 < 0.01, sprintf("** ", h1), 
            	ifelse(h1 < 0.05, sprintf("*  ", h1), sprintf("   ", h1))))))
        
    # Situation 2: model 2 has lower score, but model 1 has lower df. Is it significantly better model than model 1?
    } else if ((ml2 < ml1) & (ndf1 < ndf2)) {
        
        h1 <- pchisq(2 * (ml1 - ml2), abs(ndf1 - ndf2), lower.tail = F)
        
        out <- data.frame(Model = c(deparse(substitute(model1)), deparse(substitute(model2))), 
        	Score = c(ml1, ml2), 
        	Edf = c(ndf1, ndf2), 
        	Chisq = c("", sprintf("%.3f", ml1 - ml2)), 
        	Df = c("", sprintf("%.3f", abs(ndf1 - ndf2))), 
        	p.value = c("", ifelse(h1 < 2e-16, sprintf(" < 2e-16"), 
        		ifelse(h1 < 0.001, sprintf("%.3e", h1), 
            	ifelse(h1 < 0.01, sprintf("%.3f", h1), 
            	ifelse(h1 < 0.05, sprintf("%.3f", h1), sprintf("%.3f", h1)))))), 
            Sig. = c("", ifelse(h1 < 0.001, sprintf("***", h1), 
            	ifelse(h1 < 0.01, sprintf("** ", h1), 
            	ifelse(h1 < 0.05, sprintf("*  ", h1), sprintf("   ", h1))))))
        
    # Situation 3: model 1 has lower score, and also lower df.
    } else if ((ml1 < ml2) & (ndf1 < ndf2)) {
        advice <- sprintf("\nModel %s preferred: lower %s score (%.3f), and lower df (%.3f).\n-----\n", 
        	deparse(substitute(model1)), 
        	type, ml2 - ml1, ndf2 - ndf1)
        out <- data.frame(Model = c(deparse(substitute(model2)), deparse(substitute(model1))), 
        	Score = c(ml2, ml1), 
        	Edf = c(ndf2, ndf1), 
        	Difference = c("", sprintf("%.3f", ml2 - ml1)), 
        	Df = c("", sprintf("%.3f", ndf1 - ndf2)),
            p.value = c("",NA),
            Sign. = c("", ""))
        
    # Situation 4: model 2 has lower score, and also lower df.
    } else if ((ml2 < ml1) & (ndf2 < ndf1)) {
        advice <- sprintf("\nModel %s preferred: lower %s score (%.3f), and lower df (%.3f).\n-----\n", 
        	deparse(substitute(model2)), 
            type, ml1 - ml2, ndf1 - ndf2)
        out <- data.frame(Model = c(deparse(substitute(model1)), deparse(substitute(model2))), 
        	Score = c(ml1, ml2), 
        	Edf = c(ndf1, ndf2), 
        	Difference = c("", sprintf("%.3f", ml2 - ml1)), 
        	Df = c("", sprintf("%.3f", ndf1 - ndf2)),
            p.value = c("",NA),
            Sign. = c("", ""))
    # Other cases:
    } else {
        advice <- "No preference:\n-----\n"
        out <- data.frame(Model = c(deparse(substitute(model1)), deparse(substitute(model2))), 
        	Score = c(ml1, ml2), Edf = c(ndf1, ndf2), 
        	Difference = c("", sprintf("%.3f", ml2 - ml1)), 
        	Df = c("", sprintf("%.3f", ndf1 - ndf2)),
            p.value = c("",NA),
            Sign. = c("", ""))
    }
    rownames(out) <- NULL
    if(signif.stars==FALSE){
        out <- out[,1:6]
    }
    if(print.output){
        if(is.null(advice)){
            cat(sprintf("\nChi-square test of %s scores\n-----\n", type))
        }else{
            cat(advice)
        }
        if(is.na(out$p.value[2])){
            print(out[,1:5])
        }else{
            print(out)
        }
        
        cat("\n")
    }
    
    if (type != "AIC") {
        if (is.null(model1$AR1.rho)) {
            rho1 = 0
        } else {
            rho1 = model1$AR1.rho
        }
        
        if (is.null(model2$AR1.rho)) {
            rho2 = 0
        } else {
            rho2 = model2$AR1.rho
        }
        
        if (rho1 == 0 & rho2 == 0) {
            # AIC is useless for models with rho
            if (AIC(model1) == AIC(model2)) {
              warning <- sprintf("AIC difference: 0.\n\n")
            } else {
              warning <- sprintf("AIC difference: %.2f, model %s has lower AIC.\n\n", 
                AIC(model1) - AIC(model2), 
                ifelse(AIC(model1) >= AIC(model2), 
                    deparse(substitute(model2)), 
                    deparse(substitute(model1))))
            }
            if(print.output){
                cat(warning)
            }
        } else {
            warning(warning <- sprintf(" AIC is not reliable, because an AR1 model is included (rho1 = %f, rho2 = %f). ", 
                rho1, rho2))
        }
    }
    
    if (abs(ml1 - ml2) <= 5) {
        warning(sprintf("Only small difference in %s...\n", type))
    }
    
    invisible( list(method=type,
    	m1=list(Model=model1$formula, Score=ml1, Df=ndf1),
    	m2=list(Model=model2$formula, Score=ml2, Df=ndf2),
    	table = out,
        advice = ifelse(is.null(advice), NA, advice),
        AIC = ifelse(is.null(warning), NA, warning) ) )
}
 





#' Plot difference curve based on model predictions.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @aliases plotDiff
#' @param model A GAMM model, resulting from the functions
#' \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}.
#' @param  view Name of continuous predictor that should be plotted on the x-
#' axis.
#' @param comp Named list with the grouping predictor (categorical variable)
#' and the 2 levels to calculate the difference for.
#' @param cond A named list of the values to use for the predictor terms. 
#' Variables omitted from this list will have the closest observed value to 
#' the median for continuous variables, or the reference level for factors. 
#' @param plotCI Logical: whether or not to plot confidence intervals. 
#' Default is TRUE.
#' @param f A number to scale the standard error. Defaults to 1.96, resulting 
#' in 95\% confidence intervals. For 99\% confidence intervals use a value of 
#' 2.58.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is FALSE. Alternatively a string (or vector of strings) with the 
#' name of the random effect(s) to remove.
#' @param mark.diff Logical: whether or not marking where the difference 
#' is significantly different from 0.
#' @param eegAxis Logical: whether or not to reverse the y-axis, plotting 
#' negative values upwards. Default is FALSE.
#' @param col Line color. Shading color is derived from line color.
#' @param shade Logical: plot shaded confidence interval (TRUE) 
#' or dashed lines that indicate confidence region (FALSE).
#' @param xlim Range of x-axis. If not specified, the function automatically 
#' generates an appropriate x-axis.
#' @param ylim Range of y-axis. If not specified, the function automatically 
#' generates an appropriate y-axis.
#' @param main Text string, alternative title for plot.
#' @param xlab Text string, alternative label for x-axis. 
#' @param ylab Text string, alternative label for y-axis. 
#' @param n.grid Number of data points sampled as predictions. Defaults to 100.
#' @param add Logical: whether or not to add the line to an existing plot. 
#' Default is FALSE.
#' When no plot window is available and \code{add=TRUE}, 
#' the function will generate an error.
#' @param plot Logical: whether or not to plot the difference. If FALSE, then 
#' the output is returned as a list, with the estimated difference 
#' (\code{est}) and the standard error over the estimate (\code{se.est}) and 
#' the x-values (\code{x}). Default is TRUE.
#' @param transform.view Function for transforming 
#' the values on the x-axis. Defaults to NULL (no transformation). 
#' (See \code{\link{plot_smooth}} for more info.)
#' @param hide.label Logical: whether or not to hide the label 
#' (i.e., "difference"). Default is FALSE.
#' @param print.summary Logical: whether or not to print the summary. 
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param ... Optional arguments for plot.
#' @return If the result is not being plotted, a list is 
#' returned with the estimated difference (\code{est}) and the standard error 
#' over the estimate (\code{se}) and the x-values (\code{x}) is returned.
#' @author Martijn Wieling, Jacolien van Rij
#' @examples
#' data(simdat)
#' \dontrun{
#' m1 <- bam(Y ~ Group + te(Time, Trial, by=Group),
#'     data=simdat)
#' plot_diff(m1, view='Time', comp=list(Group=c("Children", "Adults")))
#' # Reversed y-axis (for EEG data):
#' plot_diff(m1, view='Time', comp=list(Group=c("Children", "Adults")), 
#'     eegAxis=TRUE)
#' # retrieving plot values...
#' out <- plot_diff(m1, view='Time', comp=list(Group=c("Children", "Adults")), 
#'    plot=FALSE)
#' #... which might be used for indicating differences:
#' x <- find_difference(out$est, out$se, f=1.96, xVals=out$xVals)
#' # add lines:
#' arrows(x0=x$start, x1=x$end, y0=0, y1=0,code=3, length=.1, col='red')
#' }
#'
#' @family Testing for significance
plot_diff <- function(model, view, comp, cond=NULL, plotCI=TRUE, f=1.96, 
	eegAxis=FALSE, col="black", shade=TRUE, n.grid=100, add=FALSE,
	print.summary=getOption('itsadug_print'), plot=TRUE, rm.ranef=NULL,
	main=NULL, ylab=NULL, xlab=NULL, xlim=NULL, ylim=NULL, 
	transform.view=NULL, mark.diff=TRUE, 
	hide.label=FALSE, ...) { 
	dat = model$model
	xvar <- NULL
	by_predictor <- NULL
	# check view
	if(length(view) > 1){
		warning("Only first element of 'view' is being used. Use plot_diff2 for plotting difference surfaces.")
	}else{
		xvar <- view[1]
		if(xvar %in% names(cond)){
			warning(sprintf('Predictor %s specified in view and cond. Values in cond being used, rather than the whole range of %s.', xvar))
		}else{
			cond[[xvar]] <- seq(min(na.exclude(dat[,xvar])), max(na.exclude(dat[,xvar])), length=n.grid)
		}
	}
	if(!is.null(xlim)){
        if(length(xlim) != 2){
            warning("Invalid xlim values specified. Argument xlim is being ignored.")
        }else{ 
        	cond[[xvar]] <- seq(xlim[1], xlim[2], length=n.grid)
        }
    }
	newd <- c()
	newd <- get_difference(model, comp=comp, cond=cond, 
		print.summary=print.summary, rm.ranef=rm.ranef, f=f)
	# transform values x-axis:
    errormessage <- function(){
        return("Error: the function specified in transformation.view cannot be applied to x-values, because infinite or missing values are not allowed.")
    
    }
    if(!is.null(transform.view)){
        tryCatch(newd[,xvar] <- sapply(newd[,xvar], transform.view), 
                error=function(x){},
                warning=function(x){})
        
        if(any(is.infinite(newd[,xvar])) | any(is.nan(newd[,xvar])) | any(is.na(newd[,xvar]))){
            stop(errormessage())
        }
            
        if(print.summary){
            cat("\t* Note: x-values are transformed.\n")
        }
    }
	if (is.null(main)) {
		levels1 <- paste(sapply(comp, function(x) x[1]), collapse='.')
		levels2 <- paste(sapply(comp, function(x) x[2]), collapse='.')
		main = sprintf('Difference between %s and %s', levels1, levels2)
	} 
	if(is.null(ylab)) {
		ylab = sprintf("Est. difference in %s", as.character(model$formula[[2]]))
	}
	if(is.null(xlab)) {
		xlab = xvar
	}
	if(is.null(ylim)){
		ylim <- range(newd$difference)
		if (plotCI) { 
			ylim <- with(newd, range(c(difference+CI, difference-CI)))
		}
	}
	out <- data.frame(est=newd$difference,
		x=newd[,xvar])
	names(out)[2] <- xvar
	if(plotCI){
		out$CI <- newd$CI
		out$f <- f
	}
	out$comp=list2str(names(comp), comp)
	if(plot==TRUE){
		if(add==FALSE){
			emptyPlot(range(newd[,xvar]), ylim, 
				main=main, xlab=xlab, ylab=ylab, h0=0,
				eegAxis=eegAxis, ...)
			if(hide.label==FALSE){
	            addlabel = "difference"
	            if(!is.null(rm.ranef)){
	                if(rm.ranef !=FALSE){
	                    addlabel = paste(addlabel, "excl. random", sep=", ")
	                }
	            }
	            mtext(addlabel, side=4, line=0, adj=0, 
	                cex=.75, col='gray35', xpd=TRUE)
	        }        			
		}
		if(plotCI==TRUE){
			plot_error(newd[,xvar], newd$difference, newd$CI, shade=shade, col=col, ...)
		}else{
			lines(newd[,xvar], newd$difference, col=col, ...)
		}
		if(mark.diff==TRUE){
			diff <- find_difference(newd$difference, newd$CI, newd[,xvar])
			addInterval(pos=getFigCoords('p')[3], diff$start, diff$end, col="red", lwd=2*par()$lwd, length=0, xpd=TRUE)
			abline(v=c(diff$start, diff$end), lty=3, col='red')
			if(length(diff$start) > 0){
				tmp <- data.frame(Significant = sprintf("%f - %f", diff$start, diff$end))
				if(print.summary==TRUE){
					print(tmp)
				}
			}
		}
		invisible(out)
	}else{
		return(out)
	}
}





#' Plot difference surface based on model predictions.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @aliases plotDiff2D
#' @param model A GAMM model, resulting from the functions
#' \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}.
#' @param  view Name of continuous predictors that should be plotted on the x-
#'  and y-axes. Vector of two values.
#' @param comp Named list with the grouping predictor (categorical variable)
#' and the 2 levels to calculate the difference for.
#' @param cond Named list of the values to use for the other predictor terms 
#' (not in view). 
#' @param color Colorpalette
#' @param nCol Range of colors of background of contour plot.
#' @param col Line color.
#' @param add.color.legend Logical: whether or not to add a color legend. 
#' Default is TRUE. If FALSE (omitted), one could use the function
#' \code{\link{gradientLegend}} to add a legend manually at any position.
#' @param plotCI Logical: whether or not to plot confidence intervals.
#' @param f A number to scale the standard error. Defaults to 1.96, resulting 
#' in 95\% confidence intervals. For 99\% confidence intervals use a value of 
#' 2.58.
#' @param n.grid Resolution.
#' @param nlevels Levels of contour lines.
#' @param zlim A two item array giving the lower and upper limits for the z-
#' axis scale. NULL to choose automatically.
#' @param xlim A two item array giving the lower and upper limits for the x-
#' axis scale. NULL to choose automatically.
#' @param ylim A two item array giving the lower and upper limits for the y-
#' axis scale. NULL to choose automatically.
#' @param main Title of plot.
#' @param xlab Label x-axis.
#' @param ylab Label y-axis.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is FALSE. Alternatively a string (or vector of strings) with the 
#' name of the random effect(s) to remove.
#' @param transform.view List with two functions for transforming 
#' the values on the x- and y-axis respectively. If one of the axes 
#' need to be transformed, set the other to NULL (no transformation). 
#' (See \code{\link{fvisgam}} for more info.)
#' @param print.summary Logical: whether or not to print a summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param hide.label Logical: whether or not to hide the label 
#' (i.e., "difference"). Default is FALSE.
#' @param dec Numeric: number of decimals for rounding the color legend. 
#' When NULL (default), no rounding. If -1 (default), automatically determined. 
#' Note: if value = -1 (default), rounding will be applied also when 
#' \code{zlim} is provided.
#' @param ... Optional arguments for \code{\link{plotsurface}}.
#' @return If the result is not being plotted, a list is 
#' returned with the estimated difference (\code{est}) and the standard error 
#' over the estimate (\code{se.est}) and the x-values (\code{x}) is returned.
#' @author Martijn Wieling, reimplemented by Jacolien van Rij
#'
#' @examples
#' data(simdat)
#' \dontrun{
#' m1 <- bam(Y ~ Group + te(Time, Trial, by=Group),
#'     data=simdat)
#' plot_diff2(m1, view=c('Time', 'Trial'), 
#'     comp=list(Group=c("Children", "Adults")))
#' }
#' @family Testing for significance
# plots differences in 2D plot
plot_diff2 <- function(model, view, comp, cond=NULL, 
	color='topo', nCol=100, col=NULL, add.color.legend=TRUE,
	plotCI=FALSE, f=1.96, n.grid=30, nlevels=10, 
	zlim=NULL, xlim=NULL, ylim=NULL, 
	main=NULL, xlab=NULL, ylab=NULL,
	rm.ranef=NULL, transform.view=NULL,
	hide.label=FALSE,
	dec=NULL,
	print.summary=getOption('itsadug_print'), ...) { 
	dat = model$model
	xvar <- NULL
	yvar <- NULL
	by_predictor <- NULL
	# check view
	if(length(view) < 2){
		stop('Provide predictors for x- and y-axes in view.')
	}else{
		xvar <- view[1]
		yvar <- view[2]
		if(xvar %in% names(cond)){
			warning(sprintf('Predictor %s specified in view and cond. Values in cond being used, rather than the whole range of %s.', xvar, xvar))
		}else{
			cond[[xvar]] <- seq(min(na.exclude(dat[,xvar])), max(na.exclude(dat[,xvar])), length=n.grid)
			if(!is.null(xlim)){
        		if(length(xlim) != 2){
            		warning("Invalid xlim values specified. Argument xlim is being ignored.")
        		}else{ 
            		cond[[xvar]] <- seq(xlim[1], xlim[2], length=n.grid)
        		}
    		}
		}
		if(yvar %in% names(cond)){
			warning(sprintf('Predictor %s specified in view and cond. Values in cond being used, rather than the whole range of %s.', yvar, yvar))
			cond[[yvar]] <- NULL
		}else{
			cond[[yvar]] <- seq(min(na.exclude(dat[,yvar])), max(na.exclude(dat[,yvar])), length=n.grid)
			if(!is.null(ylim)){
        		if(length(ylim) != 2){
            		warning("Invalid ylim values specified. Argument ylim is being ignored.")
        		}else{ 
            		cond[[yvar]] <- seq(ylim[1], ylim[2], length=n.grid)
        		}
    		}
		}
	}
	newd <- c()
	newd <- get_difference(model, comp=comp, cond=cond, 
		print.summary=print.summary, rm.ranef=rm.ranef, f=f)
	# transform values x- and y-axes:
    errormessage <- function(name){
        return(sprintf("Error: the function specified in transformation.view cannot be applied to %s-values, because infinite or missing values are not allowed.", name))
    
    }
    if(!is.null(transform.view)){
        if(length(transform.view)==1){
            tryCatch(newd[,xvar] <- sapply(newd[,xvar], transform.view), 
                error=function(x){},
                warning=function(x){})
            tryCatch(newd[,yvar] <- sapply(newd[,yvar], transform.view), 
                error=function(x){},
                warning=function(x){})
            if(any(is.infinite(newd[,xvar])) | any(is.nan(newd[,xvar])) | any(is.na(newd[,xvar]))){
                stop(errormessage("x"))
            }
            if(any(is.infinite(newd[,yvar])) | any(is.nan(newd[,yvar])) | any(is.na(newd[,yvar]))){
                stop(errormessage("y"))
            }
            if(print.summary){
                cat("\t* Note: The same transformation is applied to values of x-axis and y-axis.\n")
            }
        }else if(length(transform.view) >= 2){
            if(!is.null(transform.view[[1]])){
                tryCatch(newd[,xvar] <- sapply(newd[,xvar], transform.view[[1]]), 
                    error=function(x){},
                    warning=function(x){})
                if(any(is.infinite(newd[,xvar])) | any(is.nan(newd[,xvar])) | any(is.na(newd[,xvar]))){
                    stop(errormessage("x"))
                }
            }
            if(!is.null(transform.view[[2]])){
                tryCatch(newd[,yvar] <- sapply(newd[,yvar], transform.view[[2]]), 
                    error=function(x){},
                    warning=function(x){})
                if(any(is.infinite(newd[,yvar])) | any(is.nan(newd[,yvar])) | any(is.na(newd[,yvar]))){
                    stop(errormessage("y"))
                }
            }
            if(print.summary){
                cat("\t* Note: Transformation function(s) applied to values of x-axis and / or y-axis.\n")
            }
        }          
    }
	if (is.null(main)) {
		levels1 <- paste(sapply(comp, function(x) x[1]), collapse='.')
		levels2 <- paste(sapply(comp, function(x) x[2]), collapse='.')
		main = sprintf('Difference between %s and %s', levels1, levels2)
	} 
	if(is.null(ylab)) {
		ylab = view[2]
	}
	if(is.null(xlab)) {
		xlab = view[1]
	}
	if(plotCI){
		p <- plotsurface(newd, view=view, predictor="difference", valCI='CI',
			main=main, xlab=xlab, ylab=ylab, 
			zlim=zlim, 
			col=col, color=color, nCol=nCol, add.color.legend=add.color.legend,
			nlevels=nlevels, dec=dec, ...)
        if(hide.label==FALSE){
            addlabel = "difference"
            if(!is.null(rm.ranef)){
                if(rm.ranef !=FALSE){
                    addlabel = paste(addlabel, "excl. random", sep=", ")
                }
            }
            mtext(addlabel, side=4, line=0, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
        }     		
	}else{
		p <- plotsurface(newd, view=view, predictor="difference", 
			main=main, xlab=xlab, ylab=ylab, 
			zlim=zlim, 
			col=col, color=color, nCol=nCol, add.color.legend=add.color.legend,
			nlevels=nlevels, dec=dec, ...)	
        if(hide.label==FALSE){
            addlabel = "difference"
            if(!is.null(rm.ranef)){
                if(rm.ranef !=FALSE){
                    addlabel = paste(addlabel, "excl. random", sep=", ")
                }
            }
            mtext(addlabel, side=4, line=0, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
        }  	
	}
	
	p[['zlim']] <- zlim
	invisible(p)
}





#' Returns a description of the statistics of the smooth terms for reporting.
#'
#' @export
#' @import mgcv
#' @import stats
#' @description Returns a description of the statistics of the smooth terms 
#' for reporting.
#'
#' @param model A gam or bam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param summary  Optionally include the summary of the model when available, which may speed up the process for large models.
#' @param print.summary Logical: whether or not to print the output.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @examples
#' data(simdat)
#' 
#' # model without random effects:
#' m1 <- bam(Y ~ te(Time, Trial), data=simdat)
#' report_stats(m1)
#' # save report for later use:
#' report <- report_stats(m1, print.summary=FALSE)
#' report[3,2]
#'
#' @author Jacolien van Rij
#' @family Testing for significance
report_stats <- function(model, summary=NULL, print.summary=getOption('itsadug_print')){
	if(!inherits(model, "gam")){
		stop("Function is only implemented for GAMMs.")
	}
	m.sum <- NULL
	if(!is.null(summary)){
		m.sum <- summary$s.table
	}else{
		m.sum <- summary(model)$s.table
	}
	if(colnames(m.sum)[3]=="F"){
		m.sum <-  as.data.frame(m.sum)
		edf2 <- length(model$y) - sum(model$edf)
		p <- ifelse(m.sum[,4]>=.01, sprintf("p=%.3f",m.sum[,4]),
				ifelse(m.sum[,4]>=.001, "p<.01","p<.001"))
		report <- sprintf("F(%.3f, %.3f)=%.2f; %s", 
			m.sum[,1], edf2, m.sum[,3], p)
		report <- data.frame(smooth.term=rownames(m.sum),
			report=report)
		if(print.summary){
			print(report)
		}
		invisible(report)
	}else{
		m.sum <-  as.data.frame(m.sum)
		p <- ifelse(m.sum[,4]>=.01, sprintf("p=%.3f",m.sum[,4]),
				ifelse(m.sum[,4]>=.001, "p<.01","p<.001"))
		report <- sprintf("%s(%.3f)=%.2f; %s", 
			colnames(m.sum)[3], m.sum[,1],m.sum[,3], p)
		report <- data.frame(smooth.term=rownames(m.sum),
			report=report)
		if(print.summary){
			print(report)
		}
		invisible(report)
	}
}





#' Function for post-hoc comparison of the contrasts in a single GAMM model.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @description Function for post-hoc comparison of the intercept differences 
#' for different factors in a single GAMM model.
#' @param model Model, currently only implemented for models generated with 
#' \code{\link[mgcv]{bam}} or \code{\link[mgcv]{gam}}.
#' @param comp Named list with predictors (specified as names) and their levels
#'  to compare. Defaults to NULL, which returns all comparisons, 
#' unless \code{select} is specified.
#' @param select Contrast matrix for manually specified contrasts. 
#' Alternatively, a vector or list could be provided as input. 
#' See examples below. 
#' @param t.test Logical default = FALSE), whether or not to return 
#' the t-test scores instead of the Wald test. Only implemented for 
#' Gaussian models.  This option is not implemented for use with \code{select}.
#' @param null.hypothesis Numeric, value of null hypothesis. Defaults to 0 and 
#' is generally not changed.
#' @param summ Optional summary object. Defaults to NULL. For very large GAMM 
#' models it takes a long time to retrieve the summary. In these cases the 
#' summary could be provided to reduce processing time. 
#' However, it is generally recommended not to specifify a summary object, 
#' to reduce the chance of mismatch errors.
#' @param signif.stars Logical (default = TRUE). Whether or not to display 
#' stars indicating the level of significance on 95\% confidence level.
#' @param print.output Logical: whether or not to print the output. 
#' By default controlled globally in the package options:  
#' If the function \code{\link{infoMessages}} is set to TRUE, the output 
#' will be automatically printed.
#' Could be also set by explicitly providing TRUE or FALSE. See examples.
#' @section Warning:
#' This function is intended for testing intercept differences only.
#' This function compares purely the parametric components, without 
#' considering any interactions with smooth terms. So this could be 
#' considered as a partial effect comparison. For comparing the averages 
#' of conditions use \code{\link{get_difference}}, which outputs the 
#' difference in summed effects for different factor levels. 
#' @return Optionally returns a data frame with test statistics.
#' @author Petar Milin and Jacolien van Rij. 
#' @seealso \code{\link{plot_parametric}}, \code{\link{plot_diff}}, 
#' \code{\link{plot_diff2}}
#' @examples
#' data(simdat)
#' # Convert Condition to factorial predictor for illustration purposes:
#' simdat$Condition <- as.factor(simdat$Condition)
#' 
#' infoMessages("on")
#' 
#' \dontrun{
#' # some arbitrary model:
#' m1 <- bam(Y ~ Condition*Group  
#' 	+ s(Time, by=Condition) 
#' 	+ s(Time, by=Group)
#' 	+ s(Subject, bs='re'), 
#' 	data=simdat)
#' 
#' # print summary to inspect parametric terms:
#' summary(m1)
#'
#' # return all contrasts:
#' wald_gam(m1)
#'
#' # USE OF COMP
#' # return only contrasts for Adults:
#' wald_gam(m1, comp=list(Condition=levels(simdat$Condition)))
#' # return specific contrasts:
#' wald_gam(m1, comp=list(Condition=c("-1", "0", "1"), 
#'     Group=c("Adults", "Children")))
#' 
#' # USE OF SELECT
#' # Specify contrast matrix. 
#' # Note that intercept should be 0.
#' # Example: Compare Condition 0 with Conditions 2 and 3 for children.
#' # Method 1: matrix or vector:
#' R = matrix( c(0,-2,0,1,1,0,0,0,0,0,0,0), nrow=1)
#' wald_gam(m1, select=R) 
#' wald_gam(m1, select=c(0,-2,0,1,1,0,0,0,0,0,0,0)) 
#' # Method 2: list
#' # first list element are reference coefficients, 
#' # second list element are coefficients to compare
#' wald_gam(m1, select=list(2, c(4,5))) 
#' # Replication of contrasts given in summary:
#' wald_gam(m1, select=c(0,1,0,0,0,0,0,0,0,0,0,0))
#'
#' # USE OF T.TEST
#' # This option is not implemented for use with select
#' # Compare with second line of parametric summary:
#' wald_gam(m1, comp=list(Condition=c("-1", "0"), 
#'     Group="Children"), t.test=TRUE)
#' # Compare with Wald test:
#' wald_gam(m1, comp=list(Condition=c("-1", "0"), 
#'     Group="Children"))
#' 
#' # exclude significance stars:
#' wald_gam(m1, comp=list(Condition=c("-1", "0"), 
#'     Group="Children"), signif.stars=FALSE)
#' 
#' # do not print output, but save table for later use:
#' test <- wald_gam(m1, comp=list(Condition=c("-1", "0"), 
#'     Group="Children"), print.output=FALSE)
#' test
#' # alternative way:
#' infoMessages('off')
#' test2 <- wald_gam(m1, comp=list(Condition=c("-1", "0"), 
#'     Group="Children"))
#' infoMessages('on')
#'
#' }
#' @family Testing for significance
wald_gam <- function(model,
	comp=NULL, select=NULL, 
	t.test = FALSE,
	null.hypothesis=0, summ=NULL,
	signif.stars=TRUE,
    print.output=getOption('itsadug_print')) {
	convertPval <- function(x, signif.stars=TRUE){
		out   <-  	sapply(x, function(y){
			if(y < 2e-16){
				return("< 2e-16")
			}else if(y < .001){
				return(sprintf("= %.2e", y))
			}else{
				return(sprintf("= %.3f", y))
			}
		})
		stars = rep("", length(x))
		if(signif.stars==TRUE){
			stars <- sapply(x, function(y){
				if(y < .001){
					return("***")
				}else if(y < .01){
					return("**")
				}else if(y < .05){
					return("*")
				}else if(y < .1){
					return(".")
				}else{
					return("")
				}
			})
		}
		return(cbind(pvalue=out, signif=stars))
	}
	# variables:
	all = TRUE
	par.terms = NULL
	info <- sprintf("%s: %s\n", deparse(substitute(model)), paste(deparse(model$formula), collapse="\n") )
	if(is.null(summ)){
		summ <- summary(model)
	}
		
	if(!all(c("p.coeff", "cov.scaled") %in% names(summ) )){
		stop(sprintf("Function not (yet) implemented for class %s.", gsub("summary\\.", "", class(summ)[1])))
	}
	b <- summ$p.coeff
	numFix <- length(b)
	var = model$var.summary
	if(!is.null(select)){
		if(!is.null(comp)){
			warning("Both comp and select specified. Select is used instead of comp.")
			comp = NULL
		}
		all = FALSE
		if (is.list(select)){
			if(length(select) != 2){
				stop("If select is psecified as list, it much have two elements. See examples in the help file.")
			}
			tmp <- rep(0, numFix)
			el1 <- length(select[[1]])
			el2 <- length(select[[2]])
			if((min(select[[1]]) < 1) | (max(select[[1]]) > numFix)){
				stop("Select should indicate the numbers of the parametric summary. See examples in the help file.")
			}
			if((min(select[[2]]) < 1) | (max(select[[2]]) > numFix)){
				stop("Select should indicate the numbers of the parametric summary. See examples in the help file.")
			}
			tmp[select[[1]]] <- -1*el2
			tmp[select[[2]]] <- el1
			if(sum(tmp != 0) != (el1+el2)){
				stop("Select is specified incorrectly. Rows should not be selected more than once.")
			}
			select = tmp
		}
		if(is.vector(select)){
			select = matrix(select, nrow=1)
			if(select[1]>0){
				warning("Please check contrasts, as intercept is selected. See examples in the help file for using vectors to specify contrasts. If you are not sure, please use comp.")
			}
		}else if (!is.matrix(select)){
			stop("Please specify select as matrix.")
		}
		if(ncol(select) != numFix){
			stop("Number of selection is not equal to the number of parametric coefficients.")
		}	
		# calculate:
		vc <- summ$cov.scaled[1:numFix,1:numFix]
		w  <- t(select%*%b-null.hypothesis)%*%solve(select%*%vc%*%t(select))%*%(select%*%b-null.hypothesis)
		p1 <- pchisq(w[1], length(null.hypothesis), lower.tail=FALSE)
		chisq = data.frame(Chisq=w, Df=length(null.hypothesis), p.value=p1)
		pw <- convertPval(p1, signif.stars)
		# output:
		if(print.output==TRUE){
			cat(sprintf("\nWald test\n-----\n%s\n", info))
			cat("Parametric effects:\n")
			print(b)
			cat('\nContrasts:\n')
			print(select)
			cat('\nNull hypothesis = ', null.hypothesis, "\n\n")
			print(sprintf('Chi-square(%.3f) = %.3f, p %s %s', 
				length(null.hypothesis), 
				w[1], pw[1], pw[2]))
		}
		invisible(chisq)
	}else{
		# organize input:
		if(!is.null(comp)){
			par.terms <- names(comp)
			par.terms.string = list()
			for(i in par.terms){
				if(!i %in% names(var)){
					stop(sprintf("%s is not a modelterm.", i))
				}
				if(!any(inherits(var[[i]], c("factor", "character")))) {
					warning(sprintf("Predictor %s is not a factor and will not be considered.",i))
				}else{
					par.terms.string <- c(par.terms.string, sprintf("%s=c(%s)", i, 
						paste(sprintf("'%s'", comp[[i]]), collapse=',')))
				}
			}
			par.terms.string <- paste(par.terms.string, collapse=",")
			eval(parse(text=sprintf("newdat <- expand.grid(%s)", par.terms.string)))
			names(newdat) <- par.terms
		}else{
			par.terms <- attr(model$pterms, "factors")
			par.terms <- row.names(par.terms)[rowSums(par.terms) > 0]
			par.terms.string <- paste(sprintf("levels(model$var.summary[['%s']])", par.terms), collapse=",")
			eval(parse(text=sprintf("newdat <- expand.grid(%s)", par.terms.string)))
			names(newdat) <- par.terms
		}
		othersettings <- c()
		all.p.terms <- attr(model$pterms, "factors")
		all.p.terms <- row.names(all.p.terms)[rowSums(all.p.terms) > 0]
		for(i in names(var)){
			if(!i %in% par.terms){
				if(inherits(var[[i]], c("numeric", "integer"))){
					newdat[,i] <- var[[i]][2]
				}else if(inherits(var[[i]], c("factor"))){
					newdat[,i] <- as.character(var[[i]])
					if(i %in% all.p.terms){
						othersettings <- c(othersettings, i)
					}
				}else{
					newdat[,i] <- var[[i]][1]
				}
			}
		}
		fv <- predict(model, newdat, type='lpmatrix')
		el <- 1:ncol(fv)
		el <- el[el > numFix]
		fv[, el] <- 0
		out <- c()
		message <- ""
		for(i in 1:(nrow(newdat)-1)){
			for(j in (i+1):nrow(newdat)){
				cond <- data.frame(C1 = as.character(interaction(newdat[i,par.terms])),
					C2 = as.character(interaction(newdat[j,par.terms])) )
				if(length(othersettings)>0){
					cond[, othersettings] <- newdat[i, othersettings]
				}
				est <- (fv[j,]-fv[i,]) %*% coef(model)
				ci  <- sqrt(rowSums(((fv[j,]-fv[i,]) %*%vcov(model))*(fv[j,]-fv[i,]) ))
				diff <- data.frame(Estimate=est, SE=ci, CI=1.96*ci)
				# chisq:
				R  <- matrix( (fv[j,]-fv[i,])[1:numFix], nrow=1)
				vc <- summ$cov.scaled[1:numFix,1:numFix]
				w  <- t(R%*%b-null.hypothesis)%*%solve(R%*%vc%*%t(R))%*%(R%*%b-null.hypothesis)
				p1 <- pchisq(w[1], length(null.hypothesis), lower.tail=FALSE)
				chisq = data.frame(Chisq=w, Df=length(null.hypothesis), p.value=p1)
				# t-test:
				if((t.test==TRUE) & (model$family$family=="gaussian")){
					t1 <- est/ci
					p1 <- 2*pt(-abs(t1),df=summ$n-1)	
					t1 = data.frame(T=t1, n=summ$n, p.value2=p1)	
					if(length(out)> 1){
						out <- rbind(out, cbind(cond, diff, chisq, t1)	)
					}else{
						out <- cbind(cond, diff, chisq, t1)	
					}
				}else{
					t.test = FALSE
					if(length(out)> 1){
						out <- rbind(out, cbind(cond, diff, chisq)	)
					}else{
						out <- cbind(cond, diff, chisq)	
					}	
				}
			}
		}
		out <- out[order(out$Estimate),]
		# output:
		if(print.output==TRUE){
			cat(sprintf("\nWald test\n-----\n%s\n", info))
			cat("Parametric effects:\n")
			print(b)
			cat('\nNull hypothesis = ', null.hypothesis, "\n\n")
			if(length(othersettings)>0){
				for(i in othersettings){
					out[,i] <- as.character(out[,i])
				}
				if(t.test==TRUE){
					pval <- convertPval(out$p.value2, signif.stars)
					for(i in 1:nrow(out)){
						cat(sprintf('Comparing %s with %s (%s):\n\tEst. = %f, SE = %f, t-score = %.3f, p %s %s\n\n', 
						out[i,]$C1, out[i,]$C2, paste(out[i,othersettings], collapse=', '), 
						out[i,]$Estimate, out[i,]$SE, 
						out[i,]$T, pval[i,1], pval[i,2]))
					}
				}else{
					pval <- convertPval(out$p.value, signif.stars)
					for(i in 1:nrow(out)){
						cat(sprintf('Comparing %s with %s (%s):\n\tX2(%.3f) = %.3f, p %s %s\n\n', 
						out[i,]$C1, out[i,]$C2,paste(out[i,othersettings], collapse=', '),
						out[i,]$Df, out[i,]$Chisq, pval[i,1], pval[i,2]))
					}
				}
			}else{
				if(t.test==TRUE){
					pval <- convertPval(out$p.value2, signif.stars)
					cat(sprintf('Comparing %s with %s:\n\tEst. = %f, SE = %f, t-score = %.3f, p %s %s\n\n', 
						out$C1, out$C2, 
						out[i,]$Estimate, out[i,]$SE, 
						out$T, pval[,1], pval[,2]))
				}else{
					pval <- convertPval(out$p.value, signif.stars)
					cat(sprintf('Comparing %s with %s:\n\tX2(%.3f) = %.3f, p %s %s\n\n', 
						out$C1, out$C2,out$Df, out$Chisq, pval[,1], pval[,2]))
				}
			}
		}
		invisible(out)
	}
}





