#' Get coefficients for the parametric terms (intercepts and random slopes).
#'
#' @description Wrapper around the function \code{\link[stats]{coef}}, 
#' and loosely based on \code{\link[mgcv]{summary.gam}}. This function 
#' provides a much faster alternative for \code{summary(model)$p.table}.
#' The function \code{\link[mgcv]{summary.gam}}) may take considerably 
#' more time for large models, because it additionally needs to calculate 
#' estimates for the smooth term table.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param se Logical: whether or not to return the standard errors.
#' @return The coefficients of the parametric terms.
#' @examples
#' data(simdat)
#'
#' # Condition as factor, to have a random intercept
#' # for illustration purposes:
#' simdat$Condition <- as.factor(simdat$Condition)
#'
#' # Model with random effect and interactions:
#' m1 <- bam(Y ~ Group * Condition + s(Time),
#'     data=simdat)
#'
#' # extract all parametric coefficients:
#' get_coefs(m1)
#' # calculate t-values:
#' test <- get_coefs(m1)
#' test <- cbind(test, test[,1] / test[,2] )
#' colnames(test)[3] <- 't-value'
#' test
#' 
#' # get_coefs returns the same numbers as shown in the parametric summary:
#' summary(m1)
#' # get_coefs is based on the function coef. This function returns 
#' # values of all coefficients, and does not provide SE:
#' coef(m1)
#' 
#' @author Jacolien van Rij
#' @family Model predictions
get_coefs <- function (model, se=TRUE) {
    # number of parametric terms including intercept:
    n <- model$nsdf
    coefs <- coef(model)[1:n]
    p.table <- NULL
    if(se==TRUE){
        covmat <- model$Vp[1:n, 1:n]
        name <- names(coefs)
        dimnames(covmat) <- list(name, name)
        se <- diag(covmat)^.5
        p.table <- cbind(coefs, se)
        dimnames(p.table) <- list(names(coefs), 
                c("Estimate", "Std. Error"))
    }else{
        p.table <- coefs[1:n]
    }
    return(p.table)
}





#' Get model predictions for differences between conditions.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param comp A named list with the two levels to compare.
#' @param cond A named list of the values to use for the other predictor 
#' terms. Variables omitted from this list will have the closest observed 
#' value to the median for continuous variables, or the reference level for 
#' factors. 
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is FALSE. Alternatively a vector of numbers with the 
#' mdoelterm number of the random effect(s) to remove. 
#' (See notes.)
#' @param se Logical: whether or not to return the confidence interval or 
#' standard error around the estimates.
#' @param f A number to scale the standard error. Defaults to 1.96, resulting 
#' in 95\% confidence intervals. For 99\% confidence intervals use a value of 
#' 2.58.
#' @param print.summary Logical: whether or not to print a summary of the 
#' values selected for each predictor. 
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @return Returns a data frame with the estimates of the difference and 
#' optionally the confidence intervals around that estimate.
#' @section Notes:
#' Other, not specified effects and random effects are generally canceled 
#' out, when calculating the difference. When the predictors that 
#' specify the conditions to compare are involved in other interactions
#' or included as random slopes, it may be useful to specify the values
#' of other predictors with \code{cond} or remove the random effects with
#' \code{rm.ranef}. 
#' @examples
#' data(simdat)
#' 
#' # first fit a simple model:
#' m1 <- bam(Y ~ Group+te(Time, Trial, by=Group), data=simdat)
#'
#' # get difference estimates:
#' diff <- get_difference(m1, comp=list(Group=c('Adults', 'Children')), 
#'     cond=list(Time=seq(0,500,length=100)))
#' head(diff)
#' @author Jacolien van Rij, Martijn Wieling
#' @family Model predictions
get_difference <- function(model, comp, cond=NULL,
	rm.ranef=NULL,
	se=TRUE, f=1.96, print.summary=getOption('itsadug_print')){
	if(!"lm" %in% class(model)){
		stop("This function does not work for class %s models.", class(model)[1])
	}else{
		newd <- NULL
		su <- model$var.summary
		dat <- model$model
		# check comp
		if(is.null(names(comp))){
			stop("Predictor specified in 'comp' unknown. Please provide a named list for 'comp', in the form of 'comp=list(Predictor=c('level1', 'level2'))'.")
		}
		if(all(names(comp) %in% colnames(dat))){
			for(i in 1:length(comp)){
				if(length(comp[[i]]) < 2){
					stop(sprintf('Provide two levels for %s to calculate difference.', names(comp)[i]))
				}else if(length(comp[[i]]) > 2){
					warning(sprintf('More than two levels provided for predictor %s. Only first two levels are being used.', names(comp)[i]))
				}
			}
		}else{
			errname <- paste(which(!names(comp)%in% colnames(dat)), collapse=", ")
			stop(sprintf('Grouping predictor(s) not found in model: %s.', errname))
		}
		if(any(names(cond) %in% names(comp))){
			for(i in names(cond)[names(cond) %in% names(comp)] ){
				cond[[i]] <- NULL
				warning(sprintf('Predictor %s specified in comp and cond. (The value in cond will be ignored.)', i))
			}
		}
		new.cond1 <- list()
		new.cond2 <- list()
		for(i in names(su)){
			if(i %in% names(comp)){
				new.cond1[[i]] <- comp[[i]][1]
				new.cond2[[i]] <- comp[[i]][2]
			}else if(i %in% names(cond)){
				new.cond1[[i]] <- new.cond2[[i]] <- cond[[i]]
			}else{
				if(class(su[[i]])=="factor"){
					new.cond1[[i]] <- as.character(su[[i]][1])
					new.cond2[[i]] <- as.character(su[[i]][1])
				}else if(class(su[[i]])=="numeric"){
					new.cond1[[i]] <- su[[i]][2]
					new.cond2[[i]] <- su[[i]][2]
				}
			}
		}
		newd1 <- expand.grid(new.cond1)
		newd2 <- expand.grid(new.cond2)
		p1 <- mgcv::predict.gam(model, newd1, type='lpmatrix') 
		p2 <- mgcv::predict.gam(model, newd2, type='lpmatrix')
		newd <- as.data.frame(newd1)
		newd.names <- colnames(newd)
		for(nn in newd.names){
			if(nn %in% names(comp)){
				newd[,nn] <- NULL
			}
			
		}
		mysummary <- summary_data(newd, print=FALSE)	
		# Check for random effects:
		if(class(rm.ranef)=="logical"){
			if(rm.ranef[1]==FALSE){
				rm.ranef <- NULL			
			}
		}
		if(!is.null(rm.ranef)){	
			# get random effects columns:
			smoothlabels.table <- as.data.frame( do.call('rbind', 
				lapply(model$smooth, 
					function(x){
						data.frame(Label=x[['label']], 
							Dim=x[['null.space.dim']], 
							Class = attr(x, "class")[1],
							stringsAsFactors=FALSE)
					} ) ) )
			# smoothlabels <- smoothlabels[smoothlabels$Dim==0,c("Label", "Class")]
			smoothlabels <- as.vector( smoothlabels.table[smoothlabels.table$Class %in% c("random.effect","fs.interaction"), "Label"] )
			if(class(rm.ranef)=="logical"){
				if(rm.ranef[1]==TRUE){
					rm.ranef <- smoothlabels			
				}else{
					rm.ranef <- ""
				}
			}else if (inherits(rm.ranef, c("numeric", "integer"))){
				smoothlabels.table <- smoothlabels.table[rm.ranef,]
				smoothlabels <- as.vector( smoothlabels.table[smoothlabels.table$Class %in% c("random.effect","fs.interaction"), "Label"] )
			}
			rm.col <- unlist(lapply(rm.ranef, 
				function(x){
					colnames(p1)[grepl(x, colnames(p1), fixed=TRUE)]
				}))
			rm.col <- unlist(lapply(smoothlabels,
				function(x){
					rm.col[grepl(x, rm.col, fixed=TRUE)]
				}))
			# cancel random effects
			p1[,rm.col] <- 0
			p2[,rm.col] <- 0
			# find terms that only occur in random effects:
			predictors <- do.call('rbind',
				lapply(model$smooth, 
					function(x){
						data.frame(Label=x[['label']],
							Terms=x[['term']])
					} )) 	
			test <- table(predictors$Terms) - table(predictors[predictors$Label %in% rm.ranef,]$Terms)
			for(pred in names(test[test==0])){
				if(pred %in% names(mysummary)){
					mysummary[[pred]] <- paste(mysummary[[pred]], "(Might be canceled as random effect, check below.)")
				}
			}
			if(length(rm.col)>0){
				mysummary[['NOTE']] =  sprintf("The following random effects columns are canceled: %s\n", 
				paste(smoothlabels, collapse=","))
			}else{
				warning("No random effects to cancel.\n")				
			}
		}
		# calculate the difference:
		p <- p1 - p2
		newd$difference <- as.vector(p %*% coef(model))
		if(se){
			newd$CI <- f*sqrt(rowSums((p%*%vcov(model))*p))
		}
		# print summary of chosen values
		if(print.summary==TRUE){
			print_summary(mysummary)
		}	
		
		return(newd)
	}
}





#' Get model all fitted values.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param se A number to scale the standard error. Defaults to 1.96, resulting 
#' in 95\% confidence intervals. For 99\% confidence intervals use a value of 
#' 2.58.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is FALSE. Alternatively a string (or vector of strings) with the 
#' name of the random effect(s) to remove.
#' @param as.data.frame Logical: return values as data frame or as vector. 
#' Default is FALSE (returning a vector).
#' @param print.summary Logical: whether or not to print a summary of the 
#' values selected for each predictor. 
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @return A data frame with estimates and optionally errors.
#' @examples
#' data(simdat)
#' \dontrun{
#' m1 <- bam(Y ~ Group + s(Time, by=Group)+ s(Subject, bs='re'), data=simdat)
#' 
#' # as.data.frame FALSE and rm.ranef=NULL results in fitted():
#' all( get_fitted(m1) == fitted(m1) )
#' 
#' # now fitted values without random effects:
#' all( get_fitted(m1, rm.ranef=TRUE) == fitted(m1) )
#' head(get_fitted(m1, rm.ranef=TRUE))
#' 
#' # without summary:
#' infoMessages("off")
#' head(get_fitted(m1, rm.ranef=TRUE))
#' infoMessages("on")
#' }
#' @author Jacolien van Rij
#' @family Model predictions
get_fitted <- function(model, se=1.96, rm.ranef=NULL, 
	as.data.frame=FALSE,
	print.summary=getOption('itsadug_print')){
	if(!"gam" %in% class(model)){
		stop("This function does not work for class %s models.", class(model)[1])
	}else{
		if(!is.null(rm.ranef)){
			if(class(rm.ranef)=="logical"){
				if(rm.ranef[1]==FALSE){
					rm.ranef <- NULL
				}
			}
		}
		newd <- model$model
		su <- model$var.summary
		mysummary <- summary_data(newd, print=FALSE)		
		p <- mgcv::predict.gam(model, newd, type='lpmatrix')
		if(!is.null(rm.ranef)){	
			# get random effects columns:
			smoothlabels.table <- as.data.frame( do.call('rbind', 
				lapply(model$smooth, 
					function(x){
						data.frame(Label=x[['label']], 
							Dim=x[['null.space.dim']], 
							Class = attr(x, "class")[1],
							stringsAsFactors=FALSE)
					} ) ) )
			# smoothlabels <- smoothlabels[smoothlabels$Dim==0,c("Label", "Class")]
			smoothlabels <- as.vector( smoothlabels.table[smoothlabels.table$Class %in% c("random.effect","fs.interaction"), "Label"] )
			if(class(rm.ranef)=="logical"){
				if(rm.ranef[1]==TRUE){
					rm.ranef <- smoothlabels			
				}else if(rm.ranef[1]==FALSE){
					rm.ranef <- ""
				}
			}else if (inherits(rm.ranef, c("numeric", "integer"))){
				smoothlabels.table <- smoothlabels.table[rm.ranef,]
				smoothlabels <- as.vector( smoothlabels.table[smoothlabels.table$Class %in% c("random.effect","fs.interaction"), "Label"] )
			}
			rm.col <- unlist(lapply(rm.ranef, 
				function(x){
					colnames(p)[grepl(x, colnames(p), fixed=TRUE)]
				}))
			rm.col <- unlist(lapply(smoothlabels,
				function(x){
					rm.col[grepl(x, rm.col, fixed=TRUE)]
				}))
			# cancel random effects
			p[,rm.col] <- 0
			# find terms that only occur in random effects:
			predictors <- do.call('rbind',
				lapply(model$smooth, 
					function(x){
						data.frame(Label=x[['label']],
							Terms=x[['term']])
					} )) 	
			test <- table(predictors$Terms) - table(predictors[predictors$Label %in% rm.ranef,]$Terms)
			for(pred in names(test[test==0])){
				mysummary[[pred]] <- paste(mysummary[[pred]], "(Might be canceled as random effect, check below.)")
			}
			if(length(rm.col)>0){
				mysummary[['NOTE']] =  sprintf("The following random effects columns are canceled: %s\n", 
				paste(smoothlabels, collapse=","))
			}else{
				warning("No random effects to cancel.\n")				
			}
		}
		if(print.summary){
			print_summary(mysummary)
		}
		
		newd$fit <- p %*% coef(model)
		newd$CI <- se*sqrt(rowSums((p%*%vcov(model))*p))
		if(!is.null(rm.ranef)){
			newd$rm.ranef <- paste(smoothlabels, collapse=",")
		}
		if(as.data.frame==FALSE){
			return(as.vector( newd$fit ))
		}else{
			return(newd)
		}
		
	}
}





#' Get estimated for selected model terms.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param select A number, indicating the model term to be selected. 
#' @param cond A named list of the values to restrict the estimates for the 
#' predictor terms. When NULL (default) values are automatically selected.
#' Only relevant for complex interactions, which involve more than two 
#' dimensions.
#' @param n.grid Number of data points estimated for each random smooth.
#' @return A data frame with estimates for the selected smooth term.
#' @param se Logical: whether or not to return the confidence interval or 
#' standard error around the estimates.
#' @param f A number to scale the standard error. Defaults to 1.96, resulting 
#' in 95\% confidence intervals. For 99\% confidence intervals use a value of 
#' 2.58.
#' @param as.data.frame Logical: whether or not to return a data.frame. 
#' Default is false, and a list will be returned.
#' @param print.summary Logical: whether or not to print a summary of the 
#' values selected for each predictor. 
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @return A list with two or more elements:
#' \itemize{
#' \item \code{fit}: Numeric vector with the fitted values;
#' \item \code{se.fit}: Optionally, only with \code{se=TRUE}.
#' Numeric vector with the error or confidence interval values (f*SE);
#' \item \code{f}: The multiplication factor for generating 
#' the confidence interval values;
#' \item \code{terms}: Numeric vector (for 1-dimensional smooth) 
#' or data frame (more 2- or more dimensional surfaces) 
#' with values of the modelterms.
#' \item \code{title}: String with name of the model term.
#' \item \code{xlab}, \code{ylab}, or \code{labels}:
#' Labels for x-axis and optionally y-axis. Precise structure depends 
#' on type of smooth term: for 1-dimensional smooth only x-label is provided, 
#' for 2-dimensional smooths x-label and y-label are provided, 
#' for more complex smooths a vector of of labels is provided.
#' }
#' @examples
#' data(simdat)
#'
#'\dontrun{ 
#' # Model with random effect and interactions:
#' m1 <- bam(Y ~ s(Time) + s(Trial) 
#' +ti(Time, Trial)
#' +s(Time, Subject, bs='fs', m=1),
#' data=simdat)
#' 
#' # Get list with predictions:
#' p <- get_modelterm(m1, select=1)
#' emptyPlot(range(p$terms), range(p$fit), h=0)
#' plot_error(p$terms, p$fit, p$se.fit, shade=TRUE, xpd=TRUE)
#' 
#' # Plot random effects in separate panels:
#' pp <- get_modelterm(m1, select=4, as.data.frame=TRUE)
#' require(lattice)
#' lattice::xyplot(fit~Time|Subject, 
#'     data=pp, type="l",
#'     xlab="Time", ylab="Partial effect")
#' 
#' # Plot selection of subjects:
#' pp <- get_modelterm(m1, select=4, 
#'     cond=list(Subject=c('a01', 'a03', 'c16')),
#'     as.data.frame=TRUE)
#' lattice::xyplot(fit~Time|Subject, 
#'     data=pp, type="l",
#'     xlab="Time", ylab="Partial effect")
#'
#' # Or using the package ggplot2:
#' require(ggplot2)
#' pp <- get_modelterm(m1, select=4, as.data.frame=TRUE)
#' pg <- ggplot2::qplot(Time, fit, data = pp, 
#'     geom = c("point", "line"), colour = Subject)
#' pg + ggplot2::guides(col = guide_legend(nrow = 18))
#' }
#' 
#' @author Jacolien van Rij
#' @family Model predictions
get_modelterm <- function(model, select, cond=NULL, n.grid=30, 
	se=TRUE, f=1.96, as.data.frame=FALSE, 
	print.summary=getOption('itsadug_print')){
	if(!"lm" %in% class(model)){
		stop("This function does not work for class %s models.", class(model)[1])
	}else{
		# find terms:
		smoothlabels <- model$smooth[[select[1]]][['label']]
		smoothterms <-  model$smooth[[select[1]]][['term']]
		# select right grouping predictor:
		if(model$smooth[[select[1]]][['by']] !="NA"){
			cond[[model$smooth[[select[1]]][['by']]]] <- model$smooth[[select[1]]][['by.level']]
		}
		su <- model$var.summary
		newd <- NULL
		new.cond <- list()
		for(i in names(su)){
			if((i %in% names(cond)) & any(grepl(i, smoothlabels)) ){
				new.cond[[i]] <- cond[[i]]
			}else{
				if(any(grepl(i, smoothlabels))){
					if(inherits(su[[i]],"factor")){
						new.cond[[i]] <- levels(su[[i]])
					}else if(inherits(su[[i]],"numeric")){
						new.cond[[i]] <- seq(su[[i]][1], su[[i]][3], length=n.grid)
					}
				}else{
					if(inherits(su[[i]],"factor")){
						new.cond[[i]] <- as.character(su[[i]][1])
					}else if(inherits(su[[i]],"numeric")){
						new.cond[[i]] <- su[[i]][2]
					}
				}
			}
		}
		newd <- expand.grid(new.cond)
		mysummary = summary_data(newd, print=FALSE)
		p <- mgcv::predict.gam(model, newd, type='terms', se.fit=se)
		# print summary of chosen values
		if(print.summary==TRUE){
			new_summary <- list()
			for(i in names(mysummary)){
				if(i %in% smoothterms){
					new_summary[[i]] <- mysummary[[i]]
				}
			}
			print_summary(new_summary)	
		}
		if(as.data.frame){
			fv <- c()
			if(length(smoothterms) > 1){
				fv <- newd[,smoothterms]
			}else{
				fv <- data.frame(st=newd[,smoothterms])
				names(fv) <- smoothterms
			}
			if(se){
				fv$fit <- p$fit[, smoothlabels]
				fv$se.fit <- f*p$se.fit[,smoothlabels]
			}else{
				fv$fit <- p[, smoothlabels]
			}
			return(fv)
		}else{
			fv <- list()
			
			if(se){
				fv[['fit']] <- p$fit[, smoothlabels]
				fv[['se.fit']] <- f*p$se.fit[,smoothlabels]
				fv[['f']] <- f
			}else{
				fv[['fit']] <- p[, smoothlabels]
			}
			fv[['terms']] <- newd[,smoothterms]
			fv[['title']] <- model$smooth[[select]]['label']
			if(length(smoothterms)==1){
				fv[['xlab']] <- smoothterms[1]
			}else if(length(smoothterms)==2){
				fv[['xlab']] <- smoothterms[1]
				fv[['ylab']] <- smoothterms[2]
			}else{
				fv[['labels']] <- smoothterms
			}	
			return(fv)
		}
	}
}





#' Get model predictions for specific conditions.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param cond A named list of the values to use for the predictor terms. 
#' Variables omitted from this list will have the closest observed value to 
#' the median for continuous variables, or the reference level for factors. 
#' @param se Logical: whether or not to return the confidence interval or 
#' standard error around the estimates.
#' @param f A number to scale the standard error. Defaults to 1.96, resulting 
#' in 95\% confidence intervals. For 99\% confidence intervals use a value of 
#' 2.58.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is FALSE. Alternatively a vector with numbers (modelterms) 
#' of the random effect(s) to remove.
#' @param print.summary Logical: whether or not to print a summary of the 
#' values selected for each predictor. 
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @return A data frame with estimates and optionally errors.
#' @examples
#' data(simdat)
#'
#' \dontrun{
#' m1 <- bam(Y ~ Group + s(Time, by=Group), data=simdat)
#' 
#' # Time value is automatically set:
#' pp <- get_predictions(m1, cond=list(Group="Adults"))
#' head(pp)
#' 
#' # Range of time values:
#' pp <- get_predictions(m1, 
#'     cond=list(Group="Adults", Time=seq(0,500,length=100)))
#' # plot:
#' emptyPlot(500, range(pp$fit), h=0)
#' plot_error(pp$Time, pp$fit, pp$CI, shade=TRUE, xpd=TRUE)
#'
#' # Warning: also unrealistical values are possible
#' range(simdat$Time)
#' pp <- get_predictions(m1, 
#'     cond=list(Group="Adults", Time=seq(-500,0,length=100)))
#' # plot of predictions that are not supported by data:
#' emptyPlot(c(-500,0), range(pp$fit), h=0)
#' plot_error(pp$Time, pp$fit, pp$CI, shade=TRUE, xpd=TRUE) 
#' }
#'
#' @author Jacolien van Rij
#' @family Model predictions
get_predictions <- function(model, cond=NULL, se=TRUE, f=1.96, rm.ranef=NULL, 
	print.summary=getOption('itsadug_print')){
	if(is.null(cond)){
		stop("Please specify values for at least one predictor in the parameter 'cond'.")
	}
	if(!"gam" %in% class(model)){
		stop("This function does not work for class %s models.", class(model)[1])
	}else{
		if(!is.null(rm.ranef)){
			if(class(rm.ranef)=="logical"){
				if(rm.ranef[1]==FALSE){
					rm.ranef <- NULL
				}
			}
		}
		newd <- NULL
		su <- model$var.summary
		new.cond <- list()
		for(i in names(su)){
			if(i %in% names(cond)){
				new.cond[[i]] <- cond[[i]]
			}else{
				if(class(su[[i]])=="factor"){
					new.cond[[i]] <- as.character(su[[i]][1])
				}else if(class(su[[i]])=="numeric"){
					new.cond[[i]] <- su[[i]][2]
				}
			}
		}
		newd <- expand.grid(new.cond)
		mysummary <- summary_data(newd, print=FALSE)		
		p <- mgcv::predict.gam(model, newd, type='lpmatrix')
		if(!is.null(rm.ranef)){	
			# get random effects columns:
			smoothlabels.table <- as.data.frame( do.call('rbind', 
				lapply(model$smooth, 
					function(x){
						data.frame(Label=x[['label']], 
							Dim=x[['null.space.dim']], 
							Class = attr(x, "class")[1],
							stringsAsFactors=FALSE)
					} ) ) )
			# smoothlabels <- smoothlabels[smoothlabels$Dim==0,c("Label", "Class")]
			smoothlabels <- as.vector( smoothlabels.table[smoothlabels.table$Class %in% c("random.effect","fs.interaction"), "Label"] )
			if(class(rm.ranef)=="logical"){
				if(rm.ranef[1]==TRUE){
					rm.ranef <- smoothlabels			
				}else if(rm.ranef[1]==FALSE){
					rm.ranef <- ""
				}
			}else if (inherits(rm.ranef, c("numeric", "integer"))){
				smoothlabels.table <- smoothlabels.table[rm.ranef,]
				smoothlabels <- as.vector( smoothlabels.table[smoothlabels.table$Class %in% c("random.effect","fs.interaction"), "Label"] )
			}
			rm.col <- unlist(lapply(rm.ranef, 
				function(x){
					colnames(p)[grepl(x, colnames(p), fixed=TRUE)]
				}))
			rm.col <- unlist(lapply(smoothlabels,
				function(x){
					rm.col[grepl(x, rm.col, fixed=TRUE)]
				}))
			# cancel random effects
			p[,rm.col] <- 0
			# find terms that only occur in random effects:
			predictors <- do.call('rbind',
				lapply(model$smooth, 
					function(x){
						data.frame(Label=x[['label']],
							Terms=x[['term']])
					} )) 	
			test <- table(predictors$Terms) - table(predictors[predictors$Label %in% rm.ranef,]$Terms)
			for(pred in names(test[test==0])){
				mysummary[[pred]] <- paste(mysummary[[pred]], "(Might be canceled as random effect, check below.)")
			}
			if(length(rm.col)>0){
				mysummary[['NOTE']] =  sprintf("The following random effects columns are canceled: %s\n", 
				paste(smoothlabels, collapse=","))
			}else{
				warning("No random effects to cancel.\n")				
			}
		}
		if(print.summary){
			print_summary(mysummary)
		}
		
		newd$fit <- p %*% coef(model)
		if(se){
			newd$CI <- f*sqrt(rowSums((p%*%vcov(model))*p))
		}
		if(!is.null(rm.ranef)){
			newd$rm.ranef <- paste(smoothlabels, collapse=",")
		}
		
		return(newd)
	}
}





#' Get coefficients for the random intercepts and random slopes.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param cond A named list of the values to restrict the estimates for the 
#' random predictor terms. When NULL (default) all levels are returned.
#' Only relevant for complex interactions, which involve more than two 
#' dimensions.
#' @param print.summary Logical: whether or not to print a summary of the 
#' values selected for each predictor. 
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @return The coefficients of the random intercepts 
#' and slopes.
#' @examples
#' data(simdat)
#'
#' \dontrun{
#' # Condition as factor, to have a random intercept
#' # for illustration purposes:
#' simdat$Condition <- as.factor(simdat$Condition)
#'
#' # Model with random effect and interactions:
#' m2 <- bam(Y ~ s(Time) + s(Trial)
#' + ti(Time, Trial)
#' + s(Condition, bs='re')
#' + s(Time, Subject, bs='re'),
#' data=simdat)
#'
#' # extract all random effects combined:
#' newd <- get_random(m2)
#' head(newd)
#' 
#' # extract coefficients for the random intercept for Condition:
#' # Make bar plot:
#' barplot(newd[[1]])
#' abline(h=0)
#' 
#' # or select:
#' get_random(m2, cond=list(Condition=c('2','3')))
#' }
#'
#' @author Jacolien van Rij
#' @family Model predictions
get_random <- function(model, cond=NULL,
	print.summary=getOption('itsadug_print')){
	if(!"lm" %in% class(model)){
		stop("This function does not work for class %s models.", class(model)[1])
	}else{
		# find random effects:
		smoothlabels <- as.data.frame( do.call('rbind', 
			lapply(model$smooth, 
				function(x){
					data.frame(Label=x[['label']], 
						Dim=x[['null.space.dim']], 
						Class = attr(x, "class")[1],
						stringsAsFactors=FALSE)
				} ) ) )
		randomeffects <- which(smoothlabels$Class %in% c("random.effect"))
		if(length(randomeffects) == 0){
			warning("No random effects specified in the model.")
			return(NULL)
		}
		coefs <- list()
		
		coeflabels <- as.vector( smoothlabels[smoothlabels$Class %in% c("random.effect"), "Label"] )
		for(i in coeflabels){
			coefs[[i]] = coef(model)[grepl(convertNonAlphanumeric(i), names(coef(model)))]
			components = model$smooth[[which(smoothlabels$Label == i)]]$term
			components = components[sapply(components, function(x){ class(model$model[,x])})=="factor"]
			if(length(components)==1){
				names(coefs[[i]]) <- levels(model$model[,components])
				if(!is.null(cond)){
					for(j in names(cond)){
						if(j %in% components){
							coefs[[i]] = coefs[[i]][cond[[j]]]
						}
					}
				}
			}else if(length(components)>1){
				model$model[,paste(components, collapse='.')] <- interaction(model$model[,components])
				names(coefs[[i]]) = levels( model$model[,paste(components, collapse='.')] )
				if(!is.null(cond)){
					for(j in names(cond)){
						if(j %in% components){
							incl.levels <- sort(unique(model$model[model$model[,j]==cond[[j]],paste(components, collapse='.')]))
							incl.levels[incl.levels %in% names(coefs[[i]])]
							coefs[[i]] = coefs[[i]][incl.levels]
						}
					}
				}
			}
		}
		return(coefs)
		
	}
}





