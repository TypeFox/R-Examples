#' Generate N ACF plots of individual or aggregated time series.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @import graphics
#' @param x A vector with time series data, typically residuals of a 
#' regression model.
#' @param n The number of plots to generate.
#' @param split_by List of vectors (each with equal length as \code{x}) that 
#' group the values of \code{x} into trials or timeseries events. 
#' Generally other columns from the same data frame.
#' @param cond Named list with a selection of the time series events
#' specified in \code{split_by}. Default is NULL, indicating that 
#' all time series are being processed, rather than a selection.
#' @param max_lag Maximum lag at which to calculate the acf. 
#' Default is the maximum for the longest time series.
#' @param fun The function used when aggregating over time series 
#' (depending on the value of \code{split_by}).
#' @param random Logical: determine randomly which \code{n} (aggregated) time 
#' series are plotted, or use the \code{\link{quantile}} function to find a 
#' range of different time series to plot. Default is FALSE (not random).
#' @param mfrow A vector of the form c(nr, nc). The figures will be drawn in 
#' an nr-by-nc array on the device by rows.
#' @param add Logical: whether to add the plots to an exiting plot window or 
#' not. 
#' Default is FALSE.
#' @param plot Logical: whether or not to produce plot. Default is TRUE.
#' @param print.summary Logical: whether or not to print summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param ... Other arguments for plotting, see \code{\link[graphics]{par}}.
#' @return \code{n} ACF plots providing information about the autocorrelation 
#' in \code{x}.
#' @author Jacolien van Rij, R. Harald Baayen
#' @seealso Use \code{\link[stats]{acf}} for the original ACF function, 
#' and \code{\link{acf_plot}} for an ACF that takes into account time series 
#' in the data.
#' @examples
#' data(simdat)
#' 
#' # Separate ACF for each time series:
#' acf_n_plots(simdat$Y, split_by=list(simdat$Subject, simdat$Trial))
#' 
#' # Average ACF per participant:
#' acf_n_plots(simdat$Y, split_by=list(simdat$Subject))
#' 
#' \dontrun{
#' # Data treated as single time series. Plot is added to current window.
#' # Note: 1 time series results in 1 plot.
#' acf_n_plots(simdat$Y, add=TRUE)
#' # Plot 4 ACF plots doesn't work without splitting data:
#' acf_n_plots(simdat$Y, add=TRUE, n=4)
#' 
#' # Plot ACFs of 4 randomly selected time series:
#' acf_n_plots(simdat$Y, random=TRUE, n=4, add=TRUE, 
#'     split_by=list(simdat$Subject, simdat$Trial))
#'
#' }
#'
#' #---------------------------------------------
#' # When using model residuals
#' #---------------------------------------------
#' 
#' \dontrun{
#' # add missing values to simdat:
#' simdat[sample(nrow(simdat), 15),]$Y <- NA
#' # simple linear model:
#' m1 <- lm(Y ~ Time, data=simdat)
#'
#' # This will generate an error:
#' # acf_n_plots(resid(m1), split_by=list(simdat$Subject, simdat$Trial))
#'
#' # This should work:
#' el.na <- missing_est(m1)
#' acf_n_plots(resid(m1), 
#'      split_by=list(simdat[-el.na,]$Subject, simdat[-el.na,]$Trial))
#'
#' # This should also work:
#' simdat$res <- NA
#' simdat[!is.na(simdat$Y),]$res <- resid(m1)
#' acf_n_plots(simdat$res, split_by=list(simdat$Subject, simdat$Trial))
#' }
#' 
#' # see the vignette for examples:
#' vignette("acf", package="itsadug")
#' @family functions for model criticism
acf_n_plots <- function(x, n = 5, split_by = NULL, 
    cond = NULL, max_lag = NULL, fun = mean, plot=TRUE,
    random = F, mfrow = NULL, add=FALSE, 
    print.summary=getOption('itsadug_print'),...) {
    
    # get acf data:
    suppressWarnings( acfdat <- acf_plot(x, split_by = split_by, 
        cond=cond, max_lag = max_lag, fun = fun, 
        plot = FALSE, return_all = TRUE) )
    if (!nrow(acfdat$acftable) >= n) {
        warning(sprintf("Number of time series in the data (%d) is smaller than n (%d). N is reduced to %d.\n", 
            nrow(acfdat$acftable), n, nrow(acfdat$acftable)))
        n <- nrow(acfdat$acftable)
    }
    # get lag 1A vector: 
    lag1.all <- acfdat$acftable$"1"
    rn <- rownames(acfdat$acftable)
    lag1 <- lag1.all[!is.na(lag1.all)]
    nevents <- acfdat$n
    findNum <- rep(1, n)
    out <- list()
    
    if (random) {
        findNum <- sample(nrow(acfdat$acftable), size = n)
    } else {
        q <- quantile(lag1, probs = seq(0, 1, length = n))
        
        findClosestElement <- function(el, vals, num = TRUE) {
            y <- abs(vals - el)
            y <- min(y)
            if (num) {
                return(which(abs(vals - el) == y)[1])
            } else {
                return(vals[which(abs(vals - el) == y)[1]])
            }
        }
        
        findNum <- c()
        for (i in 1:n) {
            findNum <- c(findNum, findClosestElement(q[i], lag1))
            if(i > 1){
                out[[i-1]] <- list(quantile=c(q[i-1], q[i]),
                    elements = 
                        data.frame(event=rn[which(lag1.all >= q[i-1] &lag1.all < q[i])],
                            lag1 = lag1.all[which(lag1.all >= q[i-1] &lag1.all < q[i])],
                            stringsAsFactors=FALSE) )
            }
         }
        out[['quantiles']] <- q
        if(print.summary){
            cat("Quantiles to be plotted:\n")
            print(q)
        }
    }
    
    if(plot){    
        cols <- ceiling(sqrt(n))
        rows <- ceiling(n/cols)
        
        if (!is.null(mfrow)) {
            cols = mfrow[2]
            rows = mfrow[1]
        }
        
        if(add==FALSE){
            # dev.new(width = cols * 3, height = rows * 3)
            par(mfrow = c(rows, cols))
            oldmar <- par(mar = rep(3, 4))
        }
        parlist = list(...)
        if('type' %in% names(parlist)){
            parlist[['type']] <- NULL
        }
        xlab <- "Lag"
        ylab <- "ACF"
        ylim <- range(acfdat$acftable[findNum, ], na.rm = T)
        main.new <- NULL
        if('xlab' %in% names(parlist)){
            xlab <- parlist[['xlab']]
            parlist[['xlab']] <- NULL
        }
        if('ylab' %in% names(parlist)){
            ylab <- parlist[['ylab']]
            parlist[['ylab']] <- NULL
        }
        if('ylim' %in% names(parlist)){
            ylim <- parlist[['ylim']]
            parlist[['ylim']] <- NULL
        }
        if('main' %in% names(parlist)){
            main.new <- parlist[['main']]
            parlist[['main']] <- NULL
        }
        other <- paste(sprintf('%s=%s', names(parlist), parlist), collapse=',')
        
        for (j in 1:min(n, nrow(acfdat$acftable)) ) {
            main <- NULL
            if(is.null(main.new)){
                main <- paste("ACF of", row.names(acfdat$acftable[findNum[j], ]))
            }else{
                if(length(main.new) >= j){
                    main <- main.new[j]
                }else{
                    main <- main.new[1]
                }
            }
            eval(parse(text=sprintf(
                "plot(as.numeric(colnames(acfdat$acftable)), acfdat$acftable[findNum[j], ], type = 'h', 
                main=main, xlab=xlab, ylab=ylab, ylim = ylim, %s)", other)))
            ci <- -(1/nevents[findNum[j],'n'])+2/sqrt(nevents[findNum[j],'n'])
            abline(h=c(-1,1)*ci, lty=2, col='blue')
            abline(h = 0)
        }
    }
    invisible(out)    
} 





#' Generate an ACF plot of an aggregated time series.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @import graphics
#' @param x A vector with time series data, typically residuals of a 
#' regression model. 
#' (See examples for how to avoid errors due to missing values.)
#' @param split_by List of vectors (each with equal length as \code{x}) that 
#' group the values of \code{x} into trials or timeseries events. 
#' Generally other columns from the same data frame.
#' @param max_lag Maximum lag at which to calculate the acf. 
#' Default is the maximum for the longest time series.
#' @param cond Named list with a selection of the time series events
#' specified in \code{split_by}. Default is NULL, indicating that 
#' all time series are being processed, rather than a selection.
#' @param plot Logical: whether or not to plot the ACF. Default is TRUE.
#' @param fun The function used when aggregating over time series 
#' (depending on the value of \code{split_by}).
#' @param return_all Returning acfs for all time series.
#' @param ... Other arguments for plotting, see \code{\link[graphics]{par}}.
#' @return An aggregated ACF plot and / or optionally a list with the aggregated ACF values.
#' @author Jacolien van Rij
#' @seealso Use \code{\link[stats]{acf}} for the original ACF function, 
#' and \code{\link{acf_n_plots}} for inspection of individual time series.
#' @examples
#' data(simdat)
#'
#' # Default acf function:
#' acf(simdat$Y)
#' # Same plot with acf_plot:
#' acf_plot(simdat$Y)
#' # Average of ACFs per time series:
#' acf_plot(simdat$Y, split_by=list(simdat$Subject, simdat$Trial))
#' # Median of ACFs per time series:
#' acf_plot(simdat$Y, split_by=list(simdat$Subject, simdat$Trial), fun=median)
#'
#' # extract value of Lag1:
#' lag1 <- acf_plot(simdat$Y, 
#'    split_by=list(Subject=simdat$Subject, Trial=simdat$Trial), 
#'    plot=FALSE)['1']
#'
#' #---------------------------------------------
#' # When using model residuals
#' #---------------------------------------------
#' 
#' # add missing values to simdat:
#' simdat[sample(nrow(simdat), 15),]$Y <- NA
#' # simple linear model:
#' m1 <- lm(Y ~ Time, data=simdat)
#'
#' \dontrun{
#' # This will generate an error:
#' acf_plot(resid(m1), split_by=list(simdat$Subject, simdat$Trial))
#' }
#' # This should work:
#' el.na <- missing_est(m1)
#' acf_plot(resid(m1), 
#'      split_by=list(simdat[-el.na,]$Subject, simdat[-el.na,]$Trial))
#'
#' # This should also work:
#' simdat$res <- NA
#' simdat[!is.na(simdat$Y),]$res <- resid(m1)
#' acf_plot(simdat$res, split_by=list(simdat$Subject, simdat$Trial))
#'
#' # see the vignette for examples:
#' vignette("acf", package="itsadug")
#'
#' @family functions for model criticism
acf_plot <- function(x, split_by = NULL, max_lag = NULL, plot = TRUE, 
    fun = mean, cond=NULL, 
    return_all = FALSE, ...) {
    
    # check x:
    if (!is.vector(x)) {
        if (dim(x)[1] == 1) {
            x <- as.vector(x[1, ])
        } else if (dim(x)[2] == 1) {
            x <- as.vector(x[, 1])
        } else {
            stop(sprintf("%s is not a vector (dim: %d x %d).\n", deparse(substitute(x)), dim(x)[1], dim(x)[2]))
        }
    }
    xname <- deparse(substitute(x))
    plotci <- FALSE
    
    # check length x
    if (length(x) > 1) {
        # find missing values:
        n <- which(is.na(x))
        
        # check split factors:
        if (!is.null(split_by)) {
            for (i in 1:length(split_by)) {
                if (length(split_by[[i]]) != length(x)) {
                  name_split_i <- ifelse(!is.null(names(split_by)[i]), names(split_by)[i], sprintf("split_by[[%d]]", i))
                  errormessage <- sprintf("Split factor %s is not of same length as %s: %s has %d elements, %s has %d elements.\n
                    See help(acf_plot) for examples how to avoid this error.", 
                    name_split_i, 
                    deparse(substitute(x)), 
                    deparse(substitute(x)), 
                    length(x), 
                    name_split_i, 
                    length(split_by[[i]]))
                  stop(errormessage)
                }
            }
        } else {
            # warning(sprintf("No argument for split provided. %s is treated as a single time series.\n", deparse(substitute(x))))
            split_by <- as.factor(rep("Factor", length(x)))
            plotci <- TRUE
        }
        # check condition:
        if(!is.null(cond)){
            if(!is.null(split_by)){
                el <- 1:length(split_by[[1]])
                for(i in names(cond)){
                    if(i %in% names(split_by)){
                        el <- intersect(el, which(split_by[[i]] %in% cond[[i]]))
                    }else{
                        warning(sprintf('Predictor %s not specified in split_by. cond will be ignored.', i))
                    }
                }
                if(length(el) > 0){
                    x <- x[el]
                    for(i in names(split_by)){
                        if(is.factor(split_by[[i]])){
                            split_by[[i]] <- droplevels( split_by[[i]][el] )
                        }else{
                            split_by[[i]] <- split_by[[i]][el] 
                        }
                    }
                }else{
                    warning("Specified conditions not found in values of split_by. cond will be ignored.")
                }
                
            }else{
                warning('Split_by is empty, therefore cond will be ignored. Specify time series in cond for selecting specific time series.')
            }
        }
        
        # split x, and calculate average acf:
        splitdat <- split(x, f = split_by, drop = T)
        acfn <- lapply(splitdat, FUN = function(x) {
            x.rmna <- x[!is.na(x)]
            return( length(x.rmna) )
        })
        splitacf <- lapply(splitdat, FUN = function(x) {
            x.rmna <- x[!is.na(x)]
            return( acf(x.rmna, plot = F)$acf )
        })
        len <- max_lag
        if (is.null(len)) {
            len <- max(unlist(lapply(splitacf, FUN = length)), na.rm = T)
        }
        splitacf <- lapply(splitacf, FUN = function(x, max = len) {
            if (length(x) < max) {
                return(c(x, rep(NA, max - length(x))))
            } else if (length(x) > max) {
                return(x[1:max])
            } else {
                return(x)
            }
        })       
        # create wide format dfr:
        allacf <- as.data.frame(do.call("rbind", splitacf))
        names(allacf) <- (1:ncol(allacf)) - 1
        avgacf <- apply(allacf, 2, FUN = fun)  #apply(allacf, 2, FUN=fun, na.rm=T)
        
        if (plot) {
            # set plot arguments
            plot_default <- list(main = sprintf("ACF of %s", xname) , 
                xlab = "Lag", ylab = ifelse(plotci==TRUE, "ACF function (per time series)", "ACF"),
                ylim = c(min(min(avgacf, na.rm=TRUE), 0), max(max(avgacf, na.rm=TRUE), 1)), col = "black", type = "h")
            
            plot_args <- list(...)
            
            # merge plot info
            plot_args_def <- c()
            for (pa in names(plot_default)[!names(plot_default) %in% names(plot_args)]) {
                value <- plot_default[[pa]]
                if (length(value) > 1) {
                  value <- sprintf("c(%s)", paste(value, collapse = ","))
                  plot_args_def <- c(plot_args_def, paste(pa, value, sep = "="))
                } else {
                  plot_args_def <- c(plot_args_def, paste(pa, ifelse(is.character(value), sprintf("'%s'", value), value), 
                    sep = "="))
                }
            }
            
            if (length(plot_args_def) > 0) {
                eval(parse(text = paste("plot(0:(len-1), avgacf, ", paste(plot_args_def, collapse = ","), ", ...)")))
            } else {
                plot(0:(len - 1), avgacf, ...)
            }
            if(plotci){
                ci <- -(1/acfn[[1]])+2/sqrt(acfn[[1]])
                abline(h=c(-1,1)*ci, lty=2, col='blue')   
            }
            # }else{
            #     tmpn <- length(allacf[!is.na(allacf)])
            #     ci <- -(1/tmpn)+2/sqrt(tmpn)
            #     abline(h=c(-1,1)*ci, lty=2, col='blue')
            # }
            abline(h = 0)
        }
        
        # set output:
        acf_out <- avgacf
        if (return_all) {
            # create long format dfr:
            dfracf <- do.call('rbind',
                mapply(function(x, y, z){
                    data.frame(acf=x, 
                        lag=0:(length(x)-1),
                        n = rep(y, length(x)),
                        event = rep(z, length(x))) 
                    }, splitacf, acfn, names(splitacf), SIMPLIFY=FALSE, USE.NAMES=FALSE) )
            dfracf$ci <- -(1/dfracf$n) + 2/sqrt(dfracf$n)
            # add event info:
            events <- as.data.frame(split_by)
            events$event <- apply(events, 1, function(x){gsub(" ", "", paste(x, collapse="."), fixed=TRUE)})
            events <- events[!duplicated(events),]
            dfracf <- merge(dfracf, events, by='event', all.x=TRUE, all.y=FALSE)
            acfn <- do.call('rbind', lapply(names(acfn), 
                function(x){data.frame(n=acfn[[x]], event=x)}))
           
            acf_out <- list(acf = avgacf, acftable = allacf, 
                dataframe=dfracf, n=acfn, series = deparse(substitute(x)), FUN = fun)
            
        }
        invisible(acf_out)        
    } else {
        stop(sprintf("Not sufficient data to plot ACF: %s has %d elements.\n", deparse(substitute(x)), length(x)))
    }
}
 





#' Generate an ACF plot of model residuals. Works for lm, lmer, gam, bam, ....
#' 
#' @description Wrapper around \code{\link{acf_plot}} and 
#' \code{\link{acf_n_plots}} for regression models.
#' @export
#' @import mgcv
#' @import stats
#' @param model A regression model generated by \code{lm}, \code{glm},
#' \code{lmer}, \code{glmer}, \code{\link[mgcv]{gam}}, 
#' or \code{\link[mgcv]{bam}}.
#' (See examples for how to avoid errors due to missing values.)
#' @param split_pred Vector with names of model predictors that determine
#' the time series in the data, or should be used to split the ACF plot by.
#' Alternatively, \code{split_pred} can be a named list as being used by 
#' \code{\link{acf_plot}} and \code{\link{acf_n_plots}}.
#' Yet another option is to provide the text string "AR.start", for a model 
#' that includes an AR1 model. The events are derived from the AR.start column 
#' if that is provided.
#' @param n The number of plots to generate. If \code{n}=1 (default) then 
#' \code{\link{acf_plot}} is being called. If \code{n}>1 then 
#' \code{\link{acf_n_plots}} is being called.
#' @param plot Logical: whether or not to produce plot. Default is TRUE.
#' @param check.rho Numeric value: Generally leave at NULL. This value does 
#' not change anything, but it is used to check whether the model's AR1 
#' coefficient matches the expected value of rho. 
#' @param main Text string, title of plot.
#' @param ... Other arguments as input for \code{\link{acf_plot}} 
#' or \code{\link{acf_n_plots}}.
#' @return An aggregated ACF plot and / or optionally a list with the aggregated ACF values.
#' @author Jacolien van Rij
#' @seealso Use \code{\link[stats]{acf}} for the original ACF function, 
#' and \code{\link{acf_plot}}, or \code{\link{acf_n_plots}}.
#' @examples
#' data(simdat)
#'
#' # add missing values to simdat:
#' simdat[sample(nrow(simdat), 15),]$Y <- NA
#' 
#' \dontrun{
#' # Run GAMM model:
#' m1 <- bam(Y ~ te(Time, Trial)+s(Subject, bs='re'), data=simdat)
#'
#' # Using a list to split the data:
#' acf_resid(m1, split_pred=list(simdat$Subject, simdat$Trial))
#' # ...or using model predictors:
#' acf_resid(m1, split_pred=c("Subject", "Trial"))
#' 
#' # Calling acf_n_plots:
#' acf_resid(m1, split_pred=c("Subject", "Trial"), n=4)
#' # add some arguments:
#' acf_resid(m1, split_pred=c("Subject", "Trial"), n=4, max_lag=10)
#' 
#' # This does not work...
#' m2 <- lm(Y ~ Time, data=simdat)
#' acf_resid(m2, split_pred=c("Subject", "Trial"))
#' # ... but this is ok:
#' acf_resid(m2, split_pred=list(simdat$Subject, simdat$Trial))
#' 
#' # Using AR.start column:
#' simdat <- start_event(simdat, event=c("Subject", "Trial"))
#' r1 <- start_value_rho(m1)
#' m3 <- bam(Y ~ te(Time, Trial)+s(Subject, bs='re'), data=simdat, 
#'     rho=r1, AR.start=simdat$start.event)
#' acf_resid(m3, split_pred="AR.start")
#' # this is the same:
#' acf_resid(m3, split_pred=c("Subject", "Trial"))
#' # Note: use model comparison to find better value for rho
#' }
#' # see the vignette for examples:
#' vignette("acf", package="itsadug")
#' @family functions for model criticism
acf_resid <- function(model, split_pred=NULL, n=1, plot=TRUE, check.rho=NULL, main=NULL, ...){
	split_by=NULL
	res <- NULL
	if(!is.null(check.rho)){
		if( "AR1.rho" %in% names(model)) {
			mess <- "AR.start is not specified."
			if("(AR.start)" %in% names(model$model)){
				mess <- ""
			}
			if(model$AR1.rho != check.rho){
				message(sprintf("The model's value of rho %f does not match the expected value %f. %s", model$AR1.rho, check.rho, mess))
			}else{
				message(sprintf("The model's value of rho matches the expected value %f. %s", check.rho, mess))
			}
		}else{
			if(!inherits(model, "bam")){
				message("No value for rho specified in model. Use function bam() to specify rho.")
			}else if ("(AR.start)" %in% names(model$model)){
				message("AR.start is specified, but no value set for rho.")
			}else{
				message("Argument rho and AR.start not set in this model.")
			}
		}
	}
	if(!is.null(split_pred)){
		
		split_by=list()
		if(!is.list(split_pred)){
			# extract time series data from model:
			dat <- NULL
			if("lm" %in% class(model)){
				dat <- model$model
			}else if( "lmerMod" %in% class(model)){
				dat <- model@frame
			}
			if(length(split_pred)==1 & split_pred[1] %in% c("AR.start", "start.event")){
				if(split_pred %in% colnames(dat)){
					if(is.logical(dat[,split_pred])){
						split_by[[split_pred]] = derive_timeseries(dat[,split_pred])
					}else{
						split_by[[split_pred]] = dat[,split_pred]
					}
				}else{
					split_by[[split_pred]] = derive_timeseries(model)
				}
			}else if(!all(split_pred %in% colnames(dat))){
				notindata <- paste(split_pred[!split_pred %in% colnames(dat)], collapse=", ")
				stop(sprintf("split-pred value(s) %s is / are not included as predictor in the model.", 
					notindata))
			}else{
				for(i in split_pred){
					split_by[[i]] <- as.vector(dat[,i])
				}
			}
		}else{
			split_by <- split_pred
			# check for missing data
			me <- missing_est(model)
			if(!is.null(me)){
				split_by <- lapply(split_pred, function(x){
					return(x[-me])
				})
			}
		}
	}
	if(n > 1){
		if(!all(is.na(resid_gam(model, incl_na=TRUE)))){
			if(is.null(split_pred)){
				out <- acf_n_plots(resid_gam(model), split_by=split_by, n=n, plot=plot, ...)
				invisible(out)
			# }else if( !is.list(split_pred)){
			}else if( ("(AR.start)" %in% colnames(model$model)) || (( "AR1.rho" %in% names(model)) && model$AR1.rho > 0) ){
				out <- acf_n_plots(resid_gam(model, incl_na=TRUE), split_by=split_by, n=n, plot=plot, ...)
				invisible(out)
			}else{
				out <- acf_n_plots(resid_gam(model), split_by=split_by, n=n, plot=plot, ...)
				invisible(out)
			}
		}
	}else{
		if(!all(is.na(resid_gam(model, incl_na=TRUE)))){
			if(is.null(main)){
				main = ifelse(is.null(split_pred), sprintf("ACF resid_gam(%s)", deparse(substitute(model))), "ACF Average")
			}
			if(is.null(split_pred)){
				out <- acf_plot(resid_gam(model), split_by=split_by, plot=plot, main=main, ...)
				invisible(out)
			#}else if( !is.list(split_pred)){
			}else if( ("(AR.start)" %in% colnames(model$model)) || (( "AR1.rho" %in% names(model)) && model$AR1.rho > 0) ){
				out <- acf_plot(resid_gam(model, incl_na=TRUE), split_by=split_by, plot=plot, main=main, ...)
				invisible(out)
			}else{
				out <- acf_plot(resid_gam(model), split_by=split_by, plot=plot, main=main, ...)
				invisible(out)
			}
		}
	}
}





#' Derive the time series used in the AR1 model.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @param model GAMM model that includes an AR1 model.
#' @param AR.start Vector with AR.start information, 
#' necessary for the AR1 model. Optional, defaults to NULL.
#' @return A vector with time series indication based on the AR1 model.
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' 
#' # add missing values to simdat:
#' simdat[sample(nrow(simdat), 15),]$Y <- NA
#' simdat <- start_event(simdat, event=c("Subject", "Trial"))
#' 
#' \dontrun{
#' # Run GAMM model:
#' m1 <- bam(Y ~ te(Time, Trial)+s(Subject, bs='re'), data=simdat, 
#'     rho=.5, AR.start=simdat$start.event)
#' simdat$Event <- NA
#' simdat[!is.na(simdat$Y),]$Event <- derive_timeseries(m1)
#' acf_resid(m1, split_pred=list(Event=simdat$Event))
#' 
#' # And this works too:
#' simdat$Event <- derive_timeseries(simdat$start.event)
#' acf_resid(m1, split_pred=list(Event=simdat$Event))
#' 
#' # Note that acf_resid automatically makes use of derive_timeseries:
#' acf_resid(m1, split_pred="AR.start")
#' }
#' @family functions for model criticism
derive_timeseries <- function(model, AR.start=NULL){
	tmpdat <- c()
	na.el <- c()
	if(inherits(model, "lm")){
		tmpdat <- model$model
		na.el <- missing_est(model)
	}else if(inherits(model, "lmerMod")){
		tmpdat <- model@frame
		na.el <- missing_est(model)
	}else if(inherits(model, "logical")){
		tmpdat <- data.frame(val=model)
		names(tmpdat) <- "(AR.start)"
		na.el <- which(is.na(model))
	}
	if(! "(AR.start)" %in% names(tmpdat)){
		if(is.null(AR.start)){
			stop("Model does not contain an AR1 model. Please specify AR.start information.")
		}else{
			if(length(AR.start)==nrow(tmpdat)){
				tmpdat[,"(AR.start)"] <- AR.start
			}else if (length(AR.start) == (nrow(tmpdat)+length(na.el))){
				tmpdat[,"(AR.start)"] <- AR.start[-na.el]
			}else{
				stop(sprintf("AR.start has more values (%f) than observations used in the model (%f).", 
					length(AR.start), nrow(tmpdat)))
			}
		}
	}else{
		if(!is.null(AR.start)){
			warning(sprintf("Model %s includes AR1 model, AR.start will be ignored.", deparse(substitute(model))))
		}
	}
	starts <- which(tmpdat[,"(AR.start)"]==TRUE)
	ends <- move_n_point(starts - 1, n=-1)
	ends[length(ends)] <- nrow(tmpdat)
	return( as.factor(rep(1:length(starts), ends-starts+1)) )
}





#' Extract model residuals and remove the autocorrelation accounted for. 
#' 
#' @export
#' @import mgcv
#' @import stats
#' @aliases resid.gam
#' @param model A GAMM model build with \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param AR_start Optional: vector with logicals, indicating the start of 
#' events. 
#' Default is NULL, because generally the function can retrieve all necessary 
#' information from the model.
#' @param incl_na Whether or not to include missing values (NA)when returning 
#' the residuals. Defaults to FALSE.
#' @param return_all Default is FALSE. Returns a list with normal residuals, 
#' corrected residuals, and the value of rho.
#' @return Corrected residuals.
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' 
#' \dontrun{
#' # Add start event column:
#' simdat <- start_event(simdat, event=c("Subject", "Trial"))
#' head(simdat)
#' # bam model with AR1 model (toy example, not serious model):
#' m1 <- bam(Y ~ Group + te(Time, Trial, by=Group), 
#'    data=simdat, rho=.5, AR.start=simdat$start.event)
#' # Standard residuals:
#' res1 <- resid(m1)
#' # Corrected residuals:
#' res2 <- resid_gam(m1)
#' 
#' # Result in different ACF's:
#' par(mfrow=c(1,2))
#' acf(res1)
#' acf(res2)
#'
#' # Without AR.start included in the model:
#' m2 <- bam(Y ~ Group + te(Time, Trial, by=Group), 
#'    data=simdat)
#' acf(resid_gam(m2), plot=F)
#' # Same as resid(m2)!
#' acf(resid(m2), plot=F)
#'
#' ### MISSING VALUES ###
#' # Note that corrected residuals cannot be calculated for the last 
#' # point of each time series. These missing values are by default
#' # excluded.
#'
#' # Therefore, this will result in an error...
#' simdat$res <- resid_gam(m1)
#' # ... and this will give an error too:
#' simdat$res <- NA
#' simdat[!is.na(simdat$Y),] <- resid_gam(m1)
#' # ... but this works:
#' simdat$res <- resid_gam(m1, incl_na=TRUE)
#' 
#' # The parameter incl_na will NOT add missing values
#' # for missing values in the *data*. 
#' # Example:
#' simdat[sample(nrow(simdat), 15),]$Y <- NA
#' # Without AR.start included in the model:
#' m2 <- bam(Y ~ Group + te(Time, Trial, by=Group), 
#'    data=simdat)
#' # This works:
#' acf(resid_gam(m2))
#' # ...but this results in error, although no AR1 model specified:
#' simdat$res <- resid_gam(m2)
#' # ... for this type of missing data, this does not solve the problem:
#' simdat$res <- resid_gam(m2, incl_na=TRUE)
#' # instead try this:
#' simdat$res <- NA
#' simdat[-missing_est(m2),]$res <- resid_gam(m2)
#' 
#' # With AR.start included in the model:
#' m1 <- bam(Y ~ Group + te(Time, Trial, by=Group), 
#'    data=simdat, rho=.5, AR.start=simdat$start.event)
#' # This works (incl_na=TRUE):
#' simdat$res <- NA
#' simdat[-missing_est(m2),]$res <- resid_gam(m2, incl_na=TRUE)
#'
#' }
#' @seealso \code{\link[stats]{resid}}, \code{\link{missing_est}}
#' @family functions for model criticism
resid_gam <- function(model, AR_start = NULL, incl_na = FALSE, return_all = FALSE) {
    
    # message(sprintf('Version %s of package mgcv is loaded.', packageVersion('mgcv')))
    
    # help function
    next_point <- function(x) {
        x <- as.vector(x)
        
        if (length(x) > 1) 
            return(c(x[2:length(x)], NA)) else if (length(x) == 1) 
            return(NA) else return(NULL)
    }
    
    # extract time series data from model:
    tmpdat <- NULL
    n <- 1
    rho <- 0
    if("lm" %in% class(model)){
        tmpdat <- model$model
        # check AR_start:
        if (is.null(tmpdat$"(AR.start)")) {
            tmpdat$"(AR.start)" <- rep(FALSE, nrow(tmpdat))
            tmpdat[1, ]$"(AR.start)" <- TRUE
            
            if (!is.null(model$AR1.rho)) {
                if(model$AR1.rho > 0){
                    if (is.null(AR_start)) {
                        warning("No values for AR_start found in model specification. Please provide values for AR_start as argument to this function if you have run the model with an older version of mgcv (< 1.7.28).")
                    } else {
                        warning(sprintf("Values for argument AR_start may not be specified, although an AR1 model rho was included in model %s (rho = %f).", 
                          deparse(substitute(model)), model$AR1.rho))
                        tmpdat$"(AR.start)" <- AR_start
                    }
                }
            }
        }
        
        # additional check of AR.start:
        if ((!TRUE %in% unique(tmpdat$"(AR.start)")) | (!FALSE %in% unique(tmpdat$"(AR.start)"))) {
            stop("No event starts found in AR_start argument.\nImportant: check the values of the AR_start argument, add at least one event start (value TRUE), and run the model again.")
        }
        
        n <- which(tmpdat$"(AR.start)")
        if (is.null(model$AR1.rho)) {
            # warning("No rho specified in model. Assumed rho to equal 0.")
            model[["AR1.rho"]] <- 0
            rho <- 0
        }else{
            rho <- model$AR1.rho
        }
    }else if( "lmerMod" %in% class(model)){
        tmpdat <- model@frame
        rho <- 0
    }else{
        stop(sprintf('Function does not work for models of class %s.', class(model)[1]))
    }
        
    
    if (nrow(tmpdat) == length(resid(model))) {
        tmpdat$RES <- resid(model)
    } else {
        tmpdat$RES <- NA
        if(nrow(tmpdat)==length(resid(model))){
            tmpdat$RES <- resid(model)
        }else{
            tmpdat[!(1:nrow(tmpdat)) %in% model$na.action, ]$RES <- resid(model)
        } 
    }
    
    tmpdat$RES_next <- next_point(tmpdat$RES)
    if(length(n[(n - 1) > 0]) > 0){
        tmpdat[n[(n - 1) > 0], ]$RES_next <- rep(NA, length(n[(n - 1) > 0]))
    }
    
    res <- tmpdat$RES_next - rho * tmpdat$RES
    if (return_all) {
        return(list(res = tmpdat$RES, norm_res = res, AR1_rho = rho))
    }else{
        if(rho==0){
            res <- tmpdat$RES
        }
        if (incl_na) {
            return(res)
        } else {
            return(res[!is.na(res)])
        }
    }     
}
 





#' Determine the starting point for each time series.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @param data A data frame.
#' @param column Test string, name of the column that describes the order 
#' withing the time series. Default is "Time".
#' @param event A text string or vector indicating the columns that define the 
#' unique time series. Default is "Event".
#' @param label The name of the new column with the start point of each time 
#' series. Default is "start.event".
#' @param label.event In case \code{event} is not a single column, providing a 
#' text string will add a column with this name that defines unique time 
#' series. Default is NULL (no new column for time series is created).
#' @param order Logical: whether or not to order each time series. 
#' Default is TRUE, maybe set to FALSE with large data frames that are already ordered.
#' @return Data frame.
#' @author Jacolien van Rij
#' @examples 
#' data(simdat)
#' head(simdat)
#' test <- start_event(simdat, event=c("Subject", "Trial"), label.event="Event") 
#' head(test)
#' @family functions for model criticism
start_event <- function(data, column="Time", event="Event", label="start.event", label.event=NULL, order=TRUE){
	if(is.null(label)){
		stop("No output column specified in argument 'label'.")
	}
	if(!column %in% names(data)){
		stop(sprintf("No column '%s' in data frame %s.", column, deparse(substitute(data))))
	}
	if(!all(event %in% names(data))){
		el <- paste( which(event[!event %in% names(data)]), collapse=',' )
		stop(sprintf("Column name(s) '%s' not found in data frame %s.", el, deparse(substitute(data))))
	}
	if(length(column) > 1){
		warning(sprintf("Argument column has %d elements. Only first element is used.", length(column)))
		column <- column[1]
	}
	if(label %in% names(data)){
		warning(sprintf("Column %s already exists, will be overwritten.", label))
	}
	if(!is.null(label.event)){
		if(label.event %in% names(data)){
			warning(sprintf("Column %s already exists, will be overwritten.", label.event))
		}
	}
	tmp <- NULL
	if(!is.null(label.event) & length(event) > 1){
		data[, label.event] <- interaction(data[,event])
		tmp <- split(data, f=list(data[,label.event]), drop=TRUE)
	} else if (length(event) > 1){
		tmp <- split(data, f=as.list(data[,event]), drop=TRUE)
	} else {
		tmp <- split(data, f=list(data[,event]), drop=TRUE)
	}
	tmp <- lapply(tmp, function(x){
		min.x <- min(x[,column])
		x[,label] <- x[,column]==min.x
		if(order){
			x <- x[order(x[,column]),]
		}
		return(x)
	})
	if (requireNamespace("data.table", quietly = TRUE)) {
		tmp <- as.data.frame(data.table::rbindlist(tmp))
	}else{
		tmp <- do.call('rbind', tmp)
	}
	row.names(tmp) <- NULL
	return(tmp)
}





#' Extract the Lag 1 value from the ACF of the residuals of a gam, bam, lm, 
#' lmer model, ...
#' 
#' @description Wrapper around \code{\link{acf_plot}} for regression models.
#' @export
#' @import mgcv
#' @import stats
#' @param model A regression model generated by \code{lm}, \code{glm},
#' \code{lmer}, \code{glmer}, \code{\link[mgcv]{gam}}, 
#' or \code{\link[mgcv]{bam}}.
#' (See examples for how to avoid errors due to missing values.)
#' @param plot Logical: whether or not to produce plot. Default is TRUE.
#' @param lag Numeric value, indicating the lag. Default is 2. 
#' @param main Text string, title of ACF plot.
#' @param ... Other arguments for plotting the acf, see 
#' \code{\link[stats]{acf}}.
#' @return The autocorrelation value of data points with the data points 
#' at lag \code{lag}.
#' @author Jacolien van Rij
#' @seealso Use \code{\link[stats]{acf}} for the original ACF function, 
#' and \code{\link{acf_plot}}, or \code{\link{acf_resid}}.
#' @examples
#' data(simdat)
#'
#' # add missing values to simdat:
#' simdat[sample(nrow(simdat), 15),]$Y <- NA
#' 
#' \dontrun{
#' # Run GAMM model:
#' m1 <- bam(Y ~ te(Time, Trial)+s(Subject, bs='re'), data=simdat)
#'
#' # No plotting:
#' start_value_rho(m1)
#' # With plot:
#' rhom1 <- start_value_rho(m1, plot=TRUE)
#' 
#' }
#' # see the vignette for examples:
#' vignette("acf", package="itsadug")
#' @family functions for model criticism
#' 
start_value_rho <- function(model, plot=FALSE, lag=2, main=NULL,...){
	rval <- acf(resid(model), plot=plot, main=ifelse(is.null(main), sprintf("Series resid(%s)", deparse(substitute(model))), main), ...)$acf[lag]
	if(plot){
		par <- list(...)
		f <- par()$lwd
		if("lwd" %in% names(par)){
			f <- par[['lwd']]
		}
		points(lag-1, rval, col='red', type='h', lwd=1.5*f)
		abline(h=rval, col='red', lty=3)
	}
	return(rval)
}





