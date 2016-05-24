#' Inspect residuals of regression models.
#' 
#' @export
#' @param model A regression model, resulting from the functions
#' \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}, or \code{lm},
#' \code{glm}, \code{lmer}, or \code{glmer}.
#' @param  AR_start Defaults to NULL. 
#' Only use this when the model was run in an old versions of package 
#' \code{mgcv} and the function cannot retrieve the used AR.start values from 
#' the \code{model}. When an error is shown with newer versions of 
#' \code{mgcv}, please check the column provided as values of AR.start.
#' when using old versions of package \code{mgcv}.
#' Function will give error when it cannot find AR.start. 
#' @param split_pred A names list indicating time series in the data.
#' @param ask Logical: whether or not to show the plots one by one. 
#' Defaults to TRUE. When set to FALSE, make sure to have specified 
#' sufficient rows and columns to show the X plots. Alternatively, 
#' use \code{select} to plot only specific plots.
#' @param select Vector or numeric value indicating which plots to return 
#' (see Notes). Defaults to 1:4 (all).
#' @section Note:
#' \itemize{
#' \item Plot 1: distribution of residuals with QQ norm plot.
#' \item Plot 2: distribution of residuals with density plot.
#' \item Plot 3: ACF plot of residuals. In case an AR1 model is included, 
#' the gray lines indicate standard residuals, and the thick black lines 
#' indicate AR1 corrected residuals. 
#' \item Plot 4 (optional): In case the \code{split_pred} predictors are 
#' specified an ACF plot averaged over the time series is produced. dashed 
#' lines indicate the maximum and minimum time series (w.r.t. lag 2), the 
#' solid lines the 25% and 75% quantiles (w.r.t. lag 2) and the bars the 
#' mean of all time series.
#' }
#' See the examples on how to specify a selection of these plots.
#' @examples
#' data(simdat)
#'
#' \dontrun{
#' # Add start event column:
#' simdat <- start_event(simdat, event=c("Subject", "Trial"))
#' head(simdat)
#' 
#' # bam model with AR1 model (toy example, not serious model):
#' m1 <- bam(Y ~ Group + te(Time, Trial, by=Group), 
#'    data=simdat, rho=.5, AR.start=simdat$start.event)
#' 
#' # Warning, no time series specified:
#' check_resid(m1)
#' 
#' # Time series specified, results in a "standard" ACF plot, 
#' # treating all residuals as single time seriesand,
#' # and an ACF plot with the average ACF over time series:
#' check_resid(m1, split_pred=list(Subject=simdat$Subject, Trial=simdat$Trial))
#' # Note: residuals do not look very good.
#' # Alternative (results in the same, see help(acf_resid) ):
#' check_resid(m1, split_pred="AR.start")
#' 
#' # Define larger plot window (choose which line you need):
#' dev.new(width=8, height=8) # on windows or mac
#' quartz(,8,8)               # on mac
#' x11(width=8, height=8)     # on linux or mac
#' 
#' par(mfrow=c(2,2), cex=1.1)
#' check_resid(m1, split_pred="AR.start", ask=FALSE)
#' }
#'
#' @author Jacolien van Rij
#' @family Model evaluation
check_resid <- function(model, AR_start = NULL, split_pred=NULL, ask=TRUE, select=1:4){
	
	el.narm <- NULL
	res <- resid(model)
	res.rho <- NULL
	# retrieve res.rho
	if("lm" %in% class(model)){
		el.narm <- missing_est(model)
		if(!is.null(AR_start)){
			if(length(AR_start) > length(res)){
				AR_start <- AR_start[-el.narm]
			}
		}
		if(!is.null(model$AR1.rho)){
			if(model$AR1.rho > 0){
				suppressWarnings( res.rho <- resid_gam(model, incl_na=TRUE) )
			}
		}
	}else if("lmerMod" %in% class(model)){
		res <- resid(model)
	}
	old.ask <- devAskNewPage(ask=ask)
	plot.types <- c("qq", "dens", "acf1", "acf2")
	qua <- function(x){
		q <- quantile(x, probs=c(0,.25,.75,1), na.rm=TRUE)
		return(cbind(y0=q[1], y25=q[2], ym=mean(x, na.rm=TRUE), y75=q[3], y100=q[4]))
	}
	for(i in plot.types[select]){
		if(i == "qq"){
			qqnorm(res)
			qqline(res, col='red')
		}else if(i == "dens"){
			check_normaldist(res, legend.pos=NULL)
		}else if(i == "acf1"){
			acf(res, 
					main=sprintf('ACF resid(%s)', deparse(substitute(model))),
					col="darkgray",
					ylim=range(c(0, acf(res, plot=FALSE)$acf[,,1], acf(res.rho[!is.na(res.rho)], plot=FALSE)$acf[,,1]), na.rm=TRUE),
					bty='L')
			if(!is.null(res.rho)){
				acf.rho <- acf(res.rho[!is.na(res.rho)], plot=FALSE)$acf[,,1]
				lines(0:(length(acf.rho)-1), acf.rho, type='h', lwd=1.5*par()$lwd, col=1,lend=1, xpd=TRUE)
			}
		}else if(i=="acf2"){
			if(is.null(split_pred)){
				warning("Plot 4 is canceled, as no split_pred is defined.")
			}else{
				acf.res <- acf_resid(model, split_pred=split_pred, plot=FALSE, fun=qua)
				emptyPlot(ncol(acf.res)-1, range(acf.res), h0=0, bty='L',
					main=sprintf('ACF resid(%s) - average', 
						deparse(substitute(model))),
					ylab="ACF", xlab="Lag")
				lines(0:(ncol(acf.res)-1),acf.res[5,], type='h', col='darkgray')
				lines(0:(ncol(acf.res)-1), acf.res[4,], 
					col='darkgray', lty=3)
				lines(0:(ncol(acf.res)-1), acf.res[3,], 
					type="h", lwd=1.5*par()$lwd, lend=2)
				lines(0:(ncol(acf.res)-1), acf.res[3,], 
					type="p", pch=15, cex=.5)
				lines(0:(ncol(acf.res)-1), acf.res[2,], lty=3)
				lines(0:(ncol(acf.res)-1),acf.res[1,], type='h', col='darkgray')
			}
		}
	}
	devAskNewPage(ask=old.ask)
}





#' Visualization of the model fit for time series data.
#'
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @description Diagnostic plots for evaluating the model fit.
#'
#' @param model A lm or gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}, \code{\link[stats]{lm}}, \code{\link[stats]{glm}}.
#' @param plot A text string or numeric vector indicating which diagnostic 
#' plots to show. Default is 'all'. See Section 'Details' for the different 
#' options.
#' @param ask Logical: whether the user is prompted before starting a new 
#' page of output. Defaults to TRUE.
#' @param print.summary Logical: whether or not to print summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @section Details:
#' When \code{plot='all'}, the following plots are generated:
#' \enumerate{
#'   \item Residuals by fitted values. Used for inspection of general trends 
#' in the residuals.
#'   \item Residuals ordered by predictor. Useful for checking how the trends 
#' of individual predictors are captured by the model.
#'   \item Distribution of residuals. QQ plot that compares the distribution 
#' of the residuals with the normal distribution. 
#'   \item ACF of residuals. Inspection of autocorrelation in the residuals. 
#' See also \code{\link{acf_resid}}.
#'   \item Trends in the random smooths. Be careful with the interpretation of 
#' the 'fixed' effects and interactions when the random smooths show trends. 
#' See examples below.
#'   \item Printing distributions of numeric predictors.
#' }
#'
#' @examples
#' 
#' data(simdat)
#' \dontrun{
#' # no random smooths:
#' m1 <- bam(Y ~ Group + s(Time, by=Group) + s(Trial) + s(Subject, bs='re'), data=simdat)
#' diagnostics(m1)
#' 
#' # only plot residuals by predictor:
#' diagnostics(m1, plot=2)
#'
#' # without prompts:
#' par(mfrow=c(2,2))
#' diagnostics(m1, plot=1:4, ask=FALSE)
#'
#' # only plot random smooths:
#' diagnostics(m1, plot=5)
#' # Note: the plot does not change,
#' # because there are no random smooths to plot.
#'
#' # with random smooths
#' m2 <- bam(Y ~ Group + s(Time, by=Group) + s(Time, Subject, bs='fs', m=1), data=simdat)  
#' diagnostics(m2)  
#' 
#' ## INSPECTION OF RANDOM SMOOTHS 
#' ## ----------------------------
#' 
#' # In this underspecified model (too much smoothing for the interaction)
#' # part of the effect of Time is captured by the random smooths:
#' m3 <- bam(Y ~ te(Time, Trial, k=c(3,3)) + s(Time, Subject, bs='fs', m=1), data=simdat) 
#' 
#' # The plot shows a clear trend in the average of the random smooths, 
#' # and the amplitude of the mean (!) curve is almost as large as the 
#' # amplitude of the 'fixed' effect of Time:
#' diagnostics(m3, plot=5, ask=FALSE)
#'
#' # Compare with the following models:
#' m4 <- bam(Y ~ te(Time, Trial, k=c(10,5)) + s(Time, Subject, bs='fs', m=1), data=simdat) 
#' diagnostics(m4, plot=5, ask=FALSE)
#'
#' m5 <- bam(Y ~ s(Time) + s(Trial) + ti(Time, Trial) 
#'     + s(Time, Subject, bs='fs', m=1), data=simdat) 
#' diagnostics(m5, plot=5, ask=FALSE)
#' 
#' }
#' @author Jacolien van Rij 
#' @family Model evaluation
diagnostics <- function(model, plot='all', ask=TRUE, 
    print.summary=getOption('itsadug_print')){
    old <- NULL
    if(ask==TRUE){
        old <- devAskNewPage(ask=TRUE)
    }
    plot.selection <- c()
    if(!is.null(plot)){
        if(plot[1]=='all'){
            plot.selection = 1:6
        }else if(is.numeric(plot)){
            if (length(plot[plot >=1 & plot <=6]) > 0){
                plot.selection = plot
            }else{
                stop("Invalid value for argument plot. Choose 'all' or numeric value(s) ranging from 1 to 6.")
            }        
        }
    }
    # plot 1: fitted vs resid
    if( 1 %in% plot.selection ){
        if(print.summary==TRUE){
            cat("1. Residuals by fitted values...\n")
        }
        
        plot(fitted(model), resid(model),
        frame.plot=FALSE,
        main="General trends in residuals",
        xlab=sprintf("fitted(%s)", deparse(substitute(model))),
        ylab=sprintf("resid(%s)", deparse(substitute(model))) )
        abline(h=0, lty=3, col='gray')
        lines(lowess(resid(model) ~ fitted(model)), col='red')
    }
    dat <- model$model
    dat$res <- resid(model)
    # plot 2: residuals of predictors
    if( 2 %in% plot.selection ){
        if(print.summary==TRUE){
            cat("2. Residuals by predictor...\n")
        }
        for(i in names(model$var.summary)){
            if(inherits(dat[,i], c("numeric", 'integer'))){
                dat <- dat[order(dat[,i]),]
                plot(dat[,i], dat$res, 
                    frame.plot=FALSE,
                    main=sprintf("Residuals ~ %s", i),
                    xlab=i, ylab=sprintf("resid(%s)", deparse(substitute(model))))
                abline(h=0, lty=3, col='gray')
                lines(lowess(dat$res ~ dat[,i]), col='red')
            }
        }
    }
    # plot 3: distribution residuals
    if( 3 %in% plot.selection ){
        if(print.summary==TRUE){
            cat("3. Distribution of residuals...\n")
        }
        qqnorm(resid(model))
        qqline(resid(model))
    }
    # plot 4: acf
    if( 4 %in% plot.selection ){
        if(print.summary==TRUE){
            cat("4. ACF of residuals...\n")
        }
        acf_resid(model)
    }
    # get random effects columns:
    if( 5 %in% plot.selection ){
        if(print.summary==TRUE){
            cat("5. Printing averages of random smooths...\n")
        }
        smoothlabels.table <- lapply(model$smooth, 
                function(x){
                    data.frame(Label=x[['label']], 
                        Dim=x[['null.space.dim']], 
                        Class = attr(x, "class")[1],
                        Term1 = x[['term']],
                        stringsAsFactors=FALSE)
                } ) 
        smoothlabels.table <-  as.data.frame( do.call('rbind',
            mapply(function(x,y){
                x$num <- y
                return(x)
            }, smoothlabels.table , as.list(1:length(model$smooth)), 
                SIMPLIFY=FALSE)
            ) )
        # random effects
        el <- unique(smoothlabels.table[smoothlabels.table$Class == "fs.interaction",]$num)
        for ( i in el){
            ir <- inspect_random(model, select=i, fun=mean, 
                print.summary=FALSE,
                col='red', lwd=2)$values
            rangeir <- max(ir$x)- min(ir$x)
            # check if there is a base smooth:
            check <- unique(smoothlabels.table[smoothlabels.table$Term1==model$smooth[[i]]$term[1],]$num)
            check <- check[check != i]
            for(j in check){
                st <- get_modelterm(model, select=j,  print.summary=FALSE)
                range.check <- max(st$fit) - min(st$fit)
                if(print.summary){
                    cat(sprintf("\t* the amplitude of mean( %s ) is %.2f %% of the amplitude of %s\n",
                        model$smooth[[i]]$label, 
                        round(rangeir * 100/ range.check,2),
                        model$smooth[[j]]$label))
                }
            }
        }
    }
    if( 6 %in% plot.selection ){
        # distributions of predictors
        if(print.summary==TRUE){
            cat("6. Printing distributions of numeric predictors...\n")
        }
        for(i in names(model$var.summary)){
            if( inherits(dat[,i], c("numeric", 'integer'))){
                d <- density(dat[,i])
                emptyPlot(range(d$x), range(d$y),
                    main=sprintf('Density of %s',i), xlab=i)
                fill_area(d$x, d$y, col='gray', border=TRUE)
                rug_model(model, view=i)
            }
        }
    }
    if(print.summary==TRUE){
        cat('Done.\n')
    }
    if(ask==TRUE){
        devAskNewPage(ask=old)
    }   
}





#' Visualization of the model fit for time series data.
#'
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @description Plots the fitted values and the data for \code{n} 
#' trials of time series data. For example, plots \code{n} trials 
#' of the same participant.
#'
#' @param x A lm or gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}, \code{\link[stats]{lm}}, \code{\link[stats]{glm}}.
#' @param view Text string containing the predictor or column in the data 
#' to be displayed on the x-axis. 
#' Note that variables coerced to factors in the model formula 
#' won't work as view variables.  
#' @param event column name from the data 
#' that specifies the time series from which \code{n} are being plotted.
#' @param n Number of time series to plot. Default is 3. Set to -1 for plotting 
#' all time series (which may take a considerable time).
#' @param random Numeric: if set to TRUE (default), \code{n} random events are 
#' selected to plot. If set to FALSE, the first \code{n} events are selected 
#' to plot. The events could be precisely controlled with the argument 
#' \code{cond}.
#' @param cond A named list of the values to use for the other predictor terms 
#' (not in view) or to select specific trials or time series to plot. 
#' @param col Two value vector specifiying the colors for the data and 
#' the modelfit respectively.
#' @param add Logical: whether or not to add the lines to an existing plot, or 
#' start a new plot (default).
#' @param fill Logical: whether or not to fill the area between the data and 
#' the fitted values with shading. Default is FALSE.
#' @param eegAxis Logical: whether or not to reverse the y-axis, plotting the 
#' negative amplitudes upwards as traditionally is done in EEG research.
#' If eeg.axes is TRUE, labels for x- and y-axis are provided, when not 
#' provided by the user. Default value is FALSE.
#' @param main Changing the main title for the plot, see also title.
#' @param xlab Changing the label for the x axis, 
#' defaults to a description of x.
#' @param ylab Changing the label for the y axis, 
#' defaults to a description of y.
#' @param ylim the y limits of the plot.
#' @param h0 A vector indicating where to add solid horizontal lines for 
#' reference. By default no values provided.
#' @param v0 A vector indicating where to add dotted vertical lines for 
#' reference. By default no values provided.
#' @param transform Function for transforming the fitted values. 
#' Default is NULL.
#' @param hide.label Logical: whether or not to hide the label 
#' (i.e., "fitted values"). Default is FALSE.
#' @param hide.legend Logical: whether or not to hide the legend. 
#' Default is FALSE.
#' @param print.summary Logical: whether or not to print a summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param ... other options to pass on to lines and plot, 
#' see \code{\link[graphics]{par}}
#' @section Notes:
#' This function plots the fitted effects, including intercept and other 
#' predictors. 
#'
#' @examples
#' data(simdat)
#' 
#' # Create grouping predictor for time series:
#' simdat$Event <- interaction(simdat$Subject, simdat$Trial)
#' 
#' # model without random effects:
#' m1 <- bam(Y ~ te(Time, Trial),
#'     data=simdat)
#' plot_modelfit(m1, view="Time", event=simdat$Event)
#'
#' # colorizing residuals:
#' plot_modelfit(m1, view="Time", event=simdat$Event, fill=TRUE)
#' 
#' # All trials of one subject:
#' \dontrun{
#' # this produces error:
#' plot_modelfit(m1, view="Time", event=simdat$Event, 
#'     cond=list(Subject="a01"), n=-1)
#' }
#' # instead try this:
#' simdat$Subj <- ifelse(simdat$Subject=="a01", TRUE, FALSE)
#' plot_modelfit(m1, view="Time", event=simdat$Event, 
#'     cond=list(Subject=simdat$Subj), n=-1)
#' 
#' \dontrun{
#' # Model with random intercepts for subjects:
#' m2 <- bam(Y ~ te(Time, Trial)+s(Subject, bs='re'),
#'     data=simdat)
#' # now selecting a subject works, because it is in the model:
#' plot_modelfit(m2, view="Time", event=simdat$Event, 
#'     cond=list(Subject="a01"), n=-1, ylim=c(-13,13))
#'
#' # Model with random effect and interactions:
#' m3 <- bam(Y ~ te(Time, Trial)+s(Time, Subject, bs='fs', m=1),
#'     data=simdat)
#' plot_modelfit(m3, view="Time", event=simdat$Event, 
#'     cond=list(Subject="a01"), n=-1, ylim=c(-13,13))
#' }
#' @author Jacolien van Rij
#' @family Model evaluation
plot_modelfit <- function(x, view, event=NULL, 
	n=3, random=TRUE, cond = NULL, 
   	col = c(alpha(1), 'red'), add=FALSE, eegAxis=FALSE, 
   	fill=FALSE, 
    main=NULL, xlab=NULL, ylab=NULL, ylim=NULL, h0=0, v0=NULL, 
    transform=NULL, 
    hide.label=FALSE, hide.legend=FALSE, 
    print.summary=getOption('itsadug_print'), ...) {
       
    dnm <- names(list(...))
    v.names <- names(x$var.summary)
    
    if(is.null(main)){ main <- "Model fit" }
    if(is.null(xlab)){ xlab <- view }
    if(is.null(ylab)){ ylab <- as.character(x$formula[2])}
    dat <- x$model
    y <- as.character(x$formula[2])
    viewcol <- sprintf("%s%d", "tmp", sample.int(1e+06, size = 1L))
    eventcol <- sprintf("%s%d", "ev", sample.int(1e+06, size = 1L))
    plot.events <- NULL
    missing <- missing_est(x)
    if (is.null(view)) {
        stop("Specify one view predictor for the x-axis, either the name of a model predictor or a vector.")
    } else {
        if(length(view)>1){
        	if(view[1] %in% v.names){
        		view=view[1]
        		warning("Only first element of view is being used.")
        	}else{
        		stop("View is not column name of data.")
        	}
        }else{
        	if (sum(view %in% v.names) != 1) {
            	stop(paste(c("View variable must be one of", v.names), collapse = ", "))
        	}
	        if (!inherits(x$var.summary[[view]], c("numeric"))){
	            stop("Don't know what to do with parametric terms that are not simple numeric variables.")
	        }
        }
        dat[,viewcol] <- dat[,view]
    }
    if (!is.null(event)) {
        if(length(event)>1){
        	if(length(event)==nrow(dat)){
        		dat[, eventcol] <- event 
        	}else if(length(event)==(nrow(dat)+length(missing))){
        		dat[, eventcol] <- event[-missing] 
        	}else if(event[1] %in% v.names){
        		event=event[1]
        		dat[, eventcol] <- dat[,event]
        		warning("Only first element of event is being used.")
        	}else{
        		stop("Event is not column name of data.")
        	}
        }else{
            if(event[1] %in% c("AR.start", "start.event")){
                dat[, eventcol] <- derive_timeseries(x)
            }else if (sum(event %in% v.names) != 1) {
            	stop(paste(c("Event variable must be one of", v.names), collapse = ", "))
        	}else{
        		dat[, eventcol] <- dat[,event]
        	}
        }
    }
    dat$fit <- fitted(x)
    if(!is.null(cond)){
        cn <- names(cond)
        for(icn in cn){
        	if(icn %in% v.names){
        		dat <- dat[dat[,icn] %in% cond[[icn]],]
        	}else if(length(cond[[icn]])==nrow(dat)){
        		if(is.logical(cond[[icn]])){
        			dat <- dat[cond[[icn]]==TRUE,]
        		}else{
        			stop(sprintf("%s not found in the model. Provide a column of TRUE (include) and FALSE (do not include) of the size of the data.", icn))
        		}        		
        	}else if(length(cond[[icn]])==(nrow(dat)+length(missing))){
        		if(is.logical(cond[[icn]])){
        			dat <- dat[cond[[icn]][-missing]==TRUE,]
        		}else{
        			stop(sprintf("%s not found in the model. Provide a column of TRUE (include) and FALSE (do not include) of the size of the data.", icn))
        		}    
        	}else{
        		stop(sprintf("%s not found in the model. Provide a column of TRUE (include) and FALSE (do not include) of the size of the data.", icn))
        	}
        }
    }
    # sample events:
    dat <- droplevels(dat)
    if(!is.null(event)){
	    events <- unique(dat[,eventcol], na.rm=TRUE)
	    
	    if(n < 0){
	    	n <- length(events)
	    }
	    plot.events <- events[1:n]
	    if(random){    
	    	plot.events <- sample(events, min(length(events), n))
	    }
		if(is.factor(plot.events)){
	    	plot.events <- as.character(plot.events)
	    }
	    # select events:
	    dat <- droplevels( dat[dat[,eventcol] %in% plot.events,] )
	}
    if(!is.null(transform)){
        dat$fit <- sapply(dat$fit, transform)
    }
    if(is.null(ylim)){ 
        ylim <- range(c(dat$fit, dat[,y]), na.rm=TRUE)
    }
    if(add==FALSE){
        emptyPlot(range(dat[,view]), ylim,
            main=main, xlab=xlab, ylab=ylab,
            h0=h0, v0=v0, eegAxis=eegAxis, ...)
        if(hide.label==FALSE){
            addlabel = "fitted values"
            
            mtext(addlabel, side=4, line=0, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
            if(!is.null(transform)){
                mtext("transformed", side=4, line=.75, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
            }
        }      
    }
    if(names(dev.cur())[1] %in% c("X11", "postscript", "xfig", "pictex") ){
    	if(length(col[grepl("\\#.+", col)]) > 0 ){
    		hexcol <- which( grepl("\\#.+", col) )
    		col[hexcol] <- substr(col[hexcol], 1,7)
    	}
    	if(fill){
    		warning(sprintf("%s device does not allow for transparent colors. Fill color will be set to color 1.", names(dev.cur())[1]))
    	}
    }  
    if(!is.null(event)){
	    for(i in plot.events){
	    	newd <- NULL
	    	newd <- droplevels(dat[dat[, eventcol]==i,])
	    	if(nrow(newd)>1){
	    		# plot data:
	    		if(fill){
	    			fill_area(newd[,viewcol], newd$fit, from=newd[,y], 
	    				col=col[1], border=col[2], outline=FALSE)
	    			lines(newd[,viewcol], newd[,y], col=col[1], ...)
	    		}else{
	    			lines(newd[,viewcol], newd[,y], col=col[1], ...)
	    			lines(newd[,viewcol], newd$fit, col=col[2], ...)
	    		}
	    	} else if(nrow(newd)==1){
	    		# plot data:
	    		points(newd[,viewcol], newd[,y], col=col[1], ...)
	    		points(newd[,viewcol], newd$fit, col=col[2], ...)
	    	} else{
	    		if(print.summary){
	    			message(sprintf("No data for event %s. Ignored.", i))
	    		}
	    	} 
	    }
	}else{
		newd <- droplevels(dat)
    	if(nrow(newd)>1){
    		# plot data:
    		lines(newd[,viewcol], newd[,y], col=col[1], ...)
    		lines(newd[,viewcol], newd$fit, col=col[2], ...)
    	} else if(nrow(newd)==1){
    		# plot data:
    		points(newd[,viewcol], newd[,y], col=col[1], ...)
    		points(newd[,viewcol], newd$fit, col=col[2], ...)
    	} else{
    		if(print.summary){
	    		message("No data to be plotted.")
	    	}
    	} 	
	}
   	# add legend:
   	if(!hide.legend){
	   	gfc <- getFigCoords('f')
	   	legend(gfc[2], gfc[4],
	   		xjust=1, yjust=1,
	   		col=col, lwd=2, seg.len=1.5,
	   		legend=c("data", "model fit"),
	   		bty='n', cex=.85,  xpd=TRUE)
 	}
    #output
    invisible( list(events = ifelse(!is.null(event),plot.events, NA), subdat = dat[,!colnames(dat) %in% c(viewcol, eventcol)], 
    	view = dat[,viewcol], event = ifelse(!is.null(event),dat[, eventcol], NA)  ) )
}
 





