#' Return PCA predictions.
#'
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @description Produces perspective or contour plot views of gam model 
#' predictions of the additive effects interactions.
#' The code is based on the script for \code{\link[mgcv]{vis.gam}}, 
#' but allows to cancel random effects.
#'
#' @param x A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param pca.term Text string, name of model predictor that represents a 
#' principle component.
#' @param weights Named list with the predictors that are combined in the PC 
#' and their weights. See examples.
#' @param view A two-value vector containing the names of the two terms to 
#' plot. The two terms should be part of the PC. Note that 
#' variables coerced to factors in the model formula won't work as view 
#' variables.
#' @param cond A named list of the values to use for the other predictor 
#' terms (not in view). Used for choosing between smooths that share the 
#' same view predictors.
#' @param partial Logical value: whether or not to plot the partial effect 
#' (TRUE) or the summed effect (FALSE, default). 
#' @param select  A number, selecting a single model term for printing. e.g. 
#' if you want the plot for the second smooth term set select=2.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is TRUE.
#' @param se If less than or equal to zero then only the predicted surface 
#' is plotted, but if greater than zero, then 3 surfaces are plotted, one at 
#' the predicted values minus se standard errors, one at the predicted 
#' values and one at the predicted values plus se standard errors.
#' @param xlim A two item array giving the lower and upper limits for the x-
#' axis scale. NULL to choose automatically.
#' @param ylim A two item array giving the lower and upper limits for the y-
#' axis scale. NULL to choose automatically.
#' @param n.grid  The number of grid nodes in each direction used for 
#' calculating the plotted surface. 
#' @param print.summary Logical: whether or not to print a summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param as.data.frame Logical: whether the output is returned as data 
#' frame (TRUE, default) or as list (FALSE).
#' @author Jacolien van Rij
#' @seealso \code{plot_pca_surface}, \code{\link[stats]{prcomp}}
#' @examples
#' data(simdat)
#' # add hypothetical correlated term:
#' simdat$predictor <-  (simdat$Trial+10)^.75 + rnorm(nrow(simdat))
#' # principal components analysis:
#' pca <- prcomp(simdat[, c("Trial", "predictor")])
#' # only first PC term contributes:
#' summary(pca)
#' # get rotation (weights of predictors in PC):
#' pcar <- pca$rotation
#' # add PC1 to data:
#' simdat$PC1 <- pca$x[,1]
#' 
#' \dontrun{
#' # model:
#' m1 <- bam(Y ~ Group + te(Time, PC1, by=Group) 
#'     + s(Time, Subject, bs='fs', m=1, k=5), data=simdat)
#' # inspect surface:
#' fvisgam(m1, view=c("Time", "PC1"), cond=list(Group="Children"),
#'     rm.ranef=TRUE)
#' # how does Trial contribute?
#' p <- get_pca_predictions(m1, pca.term="PC1", weights=pcar[,"PC1"], 
#'     view=c("Time", "Trial"), cond=list(Group="Children"),
#'     rm.ranef=TRUE, partial=FALSE)
#' # Note that the range of Trial is estimated based on the values of PC1.
#' # A better solution is to specify the range:
#' p <- get_pca_predictions(m1, pca.term="PC1", weights=pcar[,"PC1"], 
#'     view=list(Time=range(simdat$Time), Trial=range(simdat$Trial)), 
#'     cond=list(Group="Children"),rm.ranef=TRUE, partial=FALSE)
#' # plotting of the surface:
#' plot_pca_surface(m1, pca.term="PC1", weights=pcar[,"PC1"], 
#'     view=c("Time", "Trial"), cond=list(Group="Children"),rm.ranef=TRUE)
#' }
#' @family Functions for PCA interpretation
get_pca_predictions <- function(x, pca.term=NULL, 
	weights = NULL, view = NULL, 
	cond = list(), select=NULL,
    n.grid = 30, se = 1.96,  
    xlim=NULL, ylim=NULL, 
    partial=TRUE, rm.ranef=NULL, 
    as.data.frame=TRUE, 
    print.summary=getOption('itsadug_print')){
    if(partial==TRUE){
    	if(is.null(select)){
    		stop("Specify a smooth term for which to extract the predictions with argument 'select'.")
    	}
    	if(is.null(pca.term)){
    		pca.term <- x$smooth[[1]]$term[1]
    		if(print.summary){
    			cat(sprintf("Predictor %s selected as PC term.\n", pca.term))
    		}
    	}
    }
    v.names <- names(x$var.summary)
    obs <- c()
    if(is.null(pca.term)){
    	stop("Please specify pca term.")
    }else{
    	if(!pca.term %in% v.names){
    		stop(paste(c("pca variable must be one of", v.names), collapse = ", "))
    	}
    }
    # check cond
    if(!is.null(cond)){
        cn <- names(cond)
        test <- sapply(cn, function(x){
            if(length(unique(cond[[x]]))>1){
                stop("Do not specify more than 1 value for conditions listed in the argument cond.")
            }else{
                TRUE
            }
        })
    }
    cond.0 <- cond
    if (is.null(view)) {
        stop("Specify view predictors for the x- and optionally y-axis.")
    } else {
       	if(!is.list(view)){
       		if(length(view) > 2){
            warning("View has more than two values. Only first two will be used.")
            view <- view[1:2]
          }
          view.list <- list()
          for(i in view){
            if(i %in% colnames(x$model)){
              view.list[[i]] <- range(x$model[,i], na.rm=TRUE)
            }else if(i %in% names(weights)){
              r.pc <- range(x$model[,pca.term], na.rm=TRUE)
              other.pca.components <- list()
              for(j in names(weights)){
                if(! j %in% view){
                  if(j %in% names(cond)){
                    other.pca.components[[j]] <- cond[[j]]
                  }
                }
              }
              other <- 0
              for(j in names(other.pca.components)){
                other <- other + weights[j]*other.pca.components[[j]]
              }
              other.pca.components <- NULL
              if(sum(view %in% names(weights))==1){
                # option 1: only this view predictor is in PC
                view.list[[i]] <- (r.pc -other)/ weights[i]
              }else if (sum(view %in% names(weights))==2){
                # option 2: two view predictors are in PC
                view.list[[i]] <- (r.pc - other) / sum(weights[view])
                warning("The ranges of the two view predictors are unknown, and will be assumed the same. It is in most cases better to specify the view predictors as a list, with (the range of) their values. See examples.")
              }             
            }
          }
          view = view.list
       	}
        # process view list
     		for(i in 1:length(view)){
     			if(length(view[[i]])==2){
     				view[[i]] = seq(view[[i]][1], view[[i]][2], length=n.grid)
     			}else if (length(view[[i]]) != nrow(x$model)){
     				el <- missing_est(x)
     				if(is.null(el)){
     					view[[i]] = seq(min(view[[i]], na.rm=TRUE), max(view[[i]], na.rm=TRUE), length=n.grid)
     				}else{
     					tmp = 1:length(view[[1]])
     					tmp = tmp[!tmp %in% el]
     					if(length(view[[i]][tmp])==nrow(x$model)){
     						if(!names(view)[i] %in% v.names){
     							x$model[,names(view)[i]] <- view[[i]][tmp]
     						}
     					}
     					view[[i]] = seq(min(view[[i]], na.rm=TRUE), max(view[[i]], na.rm=TRUE), length=n.grid)
     				}
     			}else{
     				if(!names(view)[i] %in% v.names){
     					x$model[,names(view)[i]] <- view[[i]]
     				}       				
     				view[[i]] = seq(min(view[[i]], na.rm=TRUE), max(view[[i]], na.rm=TRUE), length=n.grid)
     			}
     		}
       	
    }
    m1 <- view[[1]]
    m2 <- NULL
    if(length(view)>=2){
	   	m2 <- view[[2]]
	    if(length(view)>2){
	    	warning("Only first two view predictors are being used.")
	    }    	
    }
    view = names(view)[1:min(2, length(view))]
    if(!is.null(xlim)){
        if(length(xlim) != 2){
            warning("Invalid xlim values specified. Argument xlim is being ignored.")
        }else{ 
            m1 <- seq(xlim[1], xlim[2], length=n.grid)
        }
    }
    if(!is.null(ylim)){
        if(length(ylim) != 2){
            warning("Invalid ylim values specified. Argument ylim is being ignored.")
        }else if(!is.null(m2)){ 
            m2 <- seq(ylim[1], ylim[2], length=n.grid)
        }
    }
    # calculate PC1 values
    # PC = a*view1 + b*view2
    other.pca.components <- list()
    if(!any(view %in% names(weights))){
    	stop("None of the view predictors are found in weights.")
    }else if(! all(names(weights) %in% view)){
    	miss <- names(weights)[!names(weights) %in% view]
     	for(i in miss){
     		if(i %in% names(cond)){
     			other.pca.components[[i]] <- cond[[i]]
     			miss <- miss[miss != i]
     		}else{
          other.pca.components[[i]] <- 0
        }    		
    	}
    	if(length(miss)>0){
    		warning(sprintf("The following predictors are set to 0 when calculating the effect of %s: %s",
    		    			pca.term, paste(miss, collapse=", ")))
    	}
    }
    if(length(other.pca.components) > 0){
    	tmp <- 0
    	for(i in names(other.pca.components)){
    		tmp <- tmp+weights[i]*other.pca.components[[i]]
    	}
    	other.pca.components <- tmp
    }else{
    	other.pca.components <- 0
    }
    newd <- NULL
    obs <- "empty"
    if(!is.null(m2)){
  		tmp <- expand.grid(m1=m1, m2=m2)
  		names(tmp) <- view
      if(all(view %in% names(weights)) ){
          tmp[,pca.term] <- weights[view[1]]*tmp[,view[1]] + weights[view[2]]*tmp[,view[2]] + other.pca.components
          if( sum(view %in% colnames(x$model))==1){
            if( view[1] %in% colnames(x$model)){
              r.pc <- (x$model[,pca.term] - other.pca.components - weights[view[1]]*x$model[,view[1]]) / weights[view[2]]
              obs <- data.frame(x=x$model[,view[1]], y=r.pc)
              names(obs) <- view
            }else{
              r.pc <- (x$model[,pca.term] - other.pca.components - weights[view[2]]*x$model[,view[2]]) / weights[view[1]]
              obs <- data.frame(x=x$model[,view[2]], y=r.pc)
              names(obs) <- view
            }
          }else if( sum(view %in% colnames(x$model))==2 ){
            obs <- x$model[,view]
          }else{
            r.pc <- (x$model[,pca.term] - other.pca.components) / (weights[view[1]]+weights[view[2]])
            obs <- data.frame(x=r.pc, y=r.pc)
            names(obs) <- view
          }
      }else if(sum(view %in% names(weights))==1){
          tmp[,pca.term] <- weights[view[view %in% names(weights)]]*tmp[,view[view %in% names(weights)]] + other.pca.components
          if( sum(view %in% colnames(x$model))==1){
            if( view[1] %in% names(weights)){
              r.pc <- (x$model[,pca.term] - other.pca.components) / weights[view[1]]
              obs <- data.frame(x=r.pc, y=x$model[,view[2]])
              names(obs) <- view
            }else if( view[2] %in% names(weights)){
              r.pc <- (x$model[,pca.term] - other.pca.components) / weights[view[2]]
              obs <- data.frame(x=r.pc, y=x$model[,view[1]])
              names(obs) <- view
            }
          }else if( sum(view %in% colnames(x$model))==2 ){
            obs <- data.frame(x=x$model[, view[1]], y=x$model[,view[2]])
            names(obs) <- view
          }
      }		
  		cond[[pca.term]] <- unique(tmp[,pca.term])
      cond[[view[1]]] <- m1
      cond[[view[2]]] <- m2
  		if(partial==TRUE){
  			newd <- get_predictions(x, cond=cond, se=FALSE, 
  		        f=0, rm.ranef=FALSE,
  		        print.summary=FALSE)
  			newd$fit <- NULL
  			fv <- predict(x, newd, type="terms", se.fit=ifelse(se>0, TRUE, FALSE))
  			smooth.names <- x$smooth[[select]]$label
  			if(se>0){
  				newd <- cbind(newd, 
  					fit=fv$fit[,smooth.names], 
  					se.fit=fv$se.fit[,smooth.names]*se)
  			}else{
  				newd <- cbind(newd, 
  					fit=fv[,smooth.names])
  			}
  	
  		}else{
  			newd <- get_predictions(x, cond=cond, se=ifelse(se>0, TRUE, FALSE), 
  		        f=ifelse(se>0, se, 1.96), rm.ranef=rm.ranef,
  		        print.summary=print.summary)			
  		}
	    # add new predictors
      newd <- merge(newd, tmp, by=colnames(newd)[colnames(newd) %in% colnames(tmp)], all=TRUE)
	    newd <- newd[order(newd[,view[1]], newd[, view[2]]),]
	    tmp <- NULL
	}else{
		tmp <- data.frame(m1=m1)
		names(tmp) <- view
    
    tmp[,pca.term] <- weights[view[1]]*tmp[,view[1]] + other.pca.components 
    cond[[pca.term]] <- unique(tmp[,pca.term])
    cond[[view[1]]] <- m1
		if(partial==TRUE){
			newd <- get_predictions(x, cond=cond, se=FALSE, 
		        f=0, rm.ranef=FALSE,
		        print.summary=FALSE)
			newd$fit <- NULL
			fv <- predict(x, newd, type="terms", se.fit=ifelse(se>0, TRUE, FALSE))
			smooth.names <- x$smooth[[select]]$label
			if(se>0){
				newd <- cbind(newd, 
					fit=fv$fit[,smooth.names], 
					se.fit=fv$se.fit[,smooth.names]*se)
			}else{
				newd <- cbind(newd, 
					fit=fv[,smooth.names])
			}
	
		}else{
			newd <- get_predictions(x, cond=cond, se=ifelse(se>0, TRUE, FALSE), 
		        f=ifelse(se>0, se, 1.96), rm.ranef=rm.ranef,
		        print.summary=print.summary)			
		}
	    # add new predictors
	    newd[, view[1]] <- tmp[, view[1]]
	    ### newd[, view[2]] <- tmp[, view[2]] <- some way to extract these predictors
	    newd <- newd[order(newd[,view[1]]),]
	    tmp <- NULL
	}
	attr(newd, "partial") <- partial
	attr(newd, "select") <- select
	attr(newd, "rm.ranef") <- rm.ranef
	attr(newd, "se") <- se
	attr(newd, "weights") <- weights
	attr(newd, "view") <- view
	attr(newd, "obs") <- obs
	return(newd)
	
}





#' Visualization of the effect predictors in nonlinear interactions with 
#' principled components.
#'
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @description Produces perspective or contour plot views of gam model 
#' predictions of the additive effects interactions.
#' The code is based on the script for \code{\link[mgcv]{vis.gam}}, 
#' but allows to cancel random effects.
#'
#' @param x A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param pca.term Text string, name of model predictor that represents a 
#' principle component.
#' @param weights Named list with the predictors that are combined in the PC 
#' and their weights. See examples.
#' @param view A two-value vector containing the names of the two terms to 
#' plot. The two terms should be part of the PC. Note that 
#' variables coerced to factors in the model formula won't work as view 
#' variables.
#' @param cond A named list of the values to use for the other predictor 
#' terms (not in view). Used for choosing between smooths that share the 
#' same view predictors.
#' @param partial Logical value: whether or not to plot the partial effect 
#' (TRUE) or the summed effect (FALSE, default). 
#' @param select Numeric value, model term. In case \code{partial=TRUE} a 
#' model term needs to be selected.
#' @param n.grid  The number of grid nodes in each direction used for 
#' calculating the plotted surface. 
#' @param too.far Plot grid nodes that are too far from the points defined by 
#' the variables given in view can be excluded from the plot. too.far 
#' determines what is too far. The grid is scaled into the unit square along 
#' with the view variables and then grid nodes more than too.far from the 
#' predictor variables are excluded.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is TRUE.
#' @param col The colors for the facets of the plot.
#' @param color The color scheme to use for plots. One of "topo", "heat", 
#' "cm", "terrain", "gray" or "bw". 
#' @param contour.col sets the color of contours when using plot.
#' @param nCol The number of colors to use in color schemes.
#' @param plotCI Logical: whether or not to plot the confidence intervals. 
#' The value of \code{se} determines the size of the CI.
#' @param add.color.legend Logical: whether or not to add a color legend. 
#' Default is TRUE. If FALSE (omitted), one could use the function
#' \code{\link{gradientLegend}} to add a legend manually at any position.
#' @param se If less than or equal to zero then only the predicted surface 
#' is plotted, but if greater than zero, then 3 surfaces are plotted, one at 
#' the predicted values minus se standard errors, one at the predicted 
#' values and one at the predicted values plus se standard errors.
#' @param plot.type one of "contour" or "persp" (default is "contour").
#' @param zlim A two item array giving the lower and upper limits for the z-
#' axis scale. NULL to choose automatically.
#' @param xlim A two item array giving the lower and upper limits for the x-
#' axis scale. NULL to choose automatically.
#' @param ylim A two item array giving the lower and upper limits for the y-
#' axis scale. NULL to choose automatically.
#' @param print.summary Logical: whether or not to print a summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param transform Function for transforming the fitted values. 
#' Default is NULL.
#' @param transform.view List with two functions for transforming 
#' the values on the x- and y-axis respectively. If one of the axes 
#' need to be transformed, set the other to NULL (no transformation).
#' See examples below.
#' @param hide.label Logical: whether or not to hide the label 
#' (i.e., "fitted values"). Default is FALSE.
#' @param dec Numeric: number of decimals for rounding the color legend. 
#' When NULL, no rounding (default). If -1, automatically determined.  
#' Note: if value = -1, rounding will be applied also when 
#' \code{zlim} is provided.
#' @param ... other options to pass on to persp, image or contour. In 
#' particular ticktype="detailed" will add proper axes labeling to the plots.
#' @author Jacolien van Rij
#' data(simdat)
#' # add hypothetical correlated term:
#' simdat$predictor <-  (simdat$Trial+10)^.75 + rnorm(nrow(simdat))
#' # principal components analysis:
#' pca <- prcomp(simdat[, c("Trial", "predictor")])
#' # only first PC term contributes:
#' summary(pca)
#' # get rotation (weights of predictors in PC):
#' pcar <- pca$rotation
#' # add PC1 to data:
#' simdat$PC1 <- pca$x[,1]
#' 
#' # model:
#' m1 <- bam(Y ~ Group + te(Time, PC1, by=Group) 
#'     + s(Time, Subject, bs='fs', m=1, k=5), data=simdat)
#' # inspect surface:
#' fvisgam(m1, view=c("Time", "PC1"), cond=list(Group="Children"),
#'     rm.ranef=TRUE)
#' # how does Trial contribute?
#' plot_pca_surface(m1, pca.term="PC1", weights=pcar[,"PC1"], 
#'     view=c("Time", "Trial"), cond=list(Group="Children"),rm.ranef=TRUE)
#' # Note that the range of Trial is estimated based on the values of PC1.
#' # A better solution is to specify the range:
#' plot_pca_surface(m1, pca.term="PC1", weights=pcar[,"PC1"], 
#'     view=list(Time=range(simdat$Time), Trial=range(simdat$Trial)), 
#'     cond=list(Group="Children"),rm.ranef=TRUE)
#' 
#' # Partial effects:
#' pvisgam(m1, view=c("Time", "PC1"), cond=list(Group="Children"),
#'     select=1, rm.ranef=TRUE)
#' # PCA:
#' plot_pca_surface(m1, pca.term="PC1", weights=pcar[,"PC1"], 
#'     partial=TRUE, select=1,
#'     view=list(Time=range(simdat$Time), Trial=range(simdat$Trial)), 
#'     cond=list(Group="Children"))
#' 
#' @seealso \code{\link{fvisgam}}, \code{\link{pvisgam}}
#'
#' @family Functions for PCA interpretation
plot_pca_surface <- function(x, pca.term=NULL, 
	weights = NULL, view = NULL, cond = list(),
	partial=FALSE, select=NULL, se = -1, 
    n.grid = 30, too.far = 0, rm.ranef=NULL, 
    col = NA, color = "topo", contour.col = NULL,
    nCol = 50,  plotCI=FALSE,
    add.color.legend=TRUE, plot.type = "contour", 
    xlim=NULL, ylim=NULL, zlim = NULL,
    print.summary=getOption('itsadug_print'), 
    transform=NULL, transform.view=NULL, hide.label=FALSE, 
    dec=NULL, ...){
    dnm <- names(list(...))
    v.names <- names(x$var.summary)
    if(length(view) < 2){
    	stop("Specify 2 view predictors in a list; see examples.")
    }else{
    }
    newd <- get_pca_predictions(x=x, pca.term=pca.term, 
    	weights=weights, view=view,
    	cond=cond, partial=partial, select=select,
    	n.grid=n.grid, se=se, rm.ranef=rm.ranef,
    	xlim=xlim, ylim=ylim, print.summary=getOption('itsadug_print'))
    view = attr(newd, "view")
    m1 <- sort(unique(newd[,view[1]]))
    m2 <- sort(unique(newd[,view[2]]))
    x$model <- cbind(x$model, attr(newd, "obs"))
    # transform values x- and y-axes:
    errormessage <- function(name){
        return(sprintf("Error: the function specified in transformation.view cannot be applied to %s-values, because infinite or missing values are not allowed.", name))
    }
    if(!is.null(transform.view)){
        if(length(transform.view)==1){
            
            m1 <- sapply(m1, transform.view)
            m2 <- sapply(m2, transform.view)
            if(any(is.infinite(m1)) | any(is.nan(m1)) | any(is.na(m1))){
                stop(errormessage("x"))
            }
            if(any(is.infinite(m2)) | any(is.nan(m2)) | any(is.na(m2))){
                stop(errormessage("y"))
            }
            if(print.summary){
                cat("\t* Note: The same transformation is applied to values of x-axis and y-axis.\n")
            }
        }else if(length(transform.view) >= 2){
            if(!is.null(transform.view[[1]])){
                m1 <- sapply(m1, transform.view[[1]])
                if(any(is.infinite(m1)) | any(is.nan(m1)) | any(is.na(m1))){
                    stop(errormessage("x"))
                }
            }
            if(!is.null(transform.view[[2]])){
                m2 <- sapply(m2, transform.view[[2]])
                if(any(is.infinite(m2)) | any(is.nan(m2)) | any(is.na(m2))){
                    stop(errormessage("y"))
                }
            }
            if(print.summary){
                cat("\t* Note: Transformation function(s) applied to values of x-axis and / or y-axis.\n")
            }
        }          
    }
    too.far.raster <- rep(alpha('white', f=0), nrow(newd))
    newd.toofar <- newd
    ex.tf = NULL
    if (too.far > 0) {
        ex.tf <- mgcv::exclude.too.far(newd[,view[1]], newd[,view[2]], x$model[, view[1]], x$model[, view[2]], dist = too.far)
        newd.toofar$se.fit[ex.tf] <- newd.toofar$fit[ex.tf] <- NA
        too.far.raster[ex.tf] <- 'white'
    }
    # raster images are row-first, in contrast to images...
    too.far.raster <- matrix(too.far.raster, byrow=FALSE, n.grid, n.grid)
    too.far.raster <- as.raster(too.far.raster[nrow(too.far.raster):1,])
    z <- matrix(newd$fit, byrow=TRUE, n.grid, n.grid)
    z.toofar <- matrix(newd.toofar$fit, byrow=TRUE, n.grid, n.grid)
    zlab <- colnames(x$model)[!colnames(x$model) %in% names(cond)][1]
   
    if (plotCI==FALSE | se <= 0) {
        z.fit <- newd$fit
        z.fit.toofar <- newd.toofar$fit
        if(!is.null(transform)){
            z.fit <- sapply(z.fit, transform)
            z <- matrix(z.fit, byrow=TRUE, n.grid, n.grid)
            z.fit.toofar <- sapply(z.fit.toofar, transform)
        }
        old.warn <- options(warn = -1)
        av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), byrow=TRUE, n.grid, n.grid - 1)
        options(old.warn)
        max.z <- max(z, na.rm = TRUE)
        z[is.na(z)] <- max.z * 10000
        z <- matrix(z, byrow=TRUE, n.grid, n.grid)
        surf.col <- t(av) %*% z %*% av
        surf.col[surf.col > max.z * 2] <- NA
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
                stop("Something wrong with zlim")
            if(!is.null(dec)){
                if(dec == -1){
                    dec <- getDec(min(zlim, na.rm=TRUE))
                }
                zlim <- getRange(zlim, step=(.1^dec), n.seg=2)
            }
            min.z <- zlim[1]
            max.z <- zlim[2]
        } else {
            if(!is.null(dec)){
                if(dec == -1){
                    dec <- getDec(min(z.fit.toofar, na.rm = TRUE))
                }
                tmp <- getRange(range(z.fit.toofar, na.rm = TRUE), n.seg=2, step=(.1^dec))
            }else{
                tmp <- range(z.fit.toofar, na.rm = TRUE)
            }
            # min.z <- min(z.fit, na.rm = TRUE)
            # max.z <- max(z.fit, na.rm = TRUE)
            min.z <- tmp[1]
            max.z <- tmp[2]
        }
        surf.col <- surf.col - min.z
        surf.col <- surf.col/(max.z - min.z)
        surf.col <- round(surf.col * nCol)
        con.col <- 1
        if (color == "heat") {
            pal <- heat.colors(nCol)
            con.col <- 3
        } else if (color == "topo") {
            pal <- topo.colors(nCol)
            con.col <- 2
        } else if (color == "cm") {
            pal <- cm.colors(nCol)
            con.col <- 1
        } else if (color == "terrain") {
            pal <- terrain.colors(nCol)
            con.col <- 2
        } else if (color == "bpy") {
            if (requireNamespace("sp", quietly = TRUE)) {
                pal <- sp::bpy.colors(nCol)
                con.col <- 1
            } else {
                warning("Package 'sp' needed for bpy color palette. Using topo.colors instead (default).")
                color <- 'topo'
                pal <- topo.colors(nCol)
                con.col <- 2
            }
        } else if (color == "gray" || color == "bw") {
            pal <- gray(seq(0.1, 0.9, length = nCol))
            con.col <- 1
        } else stop("color scheme not recognized")
        if (is.null(contour.col)) 
            contour.col <- con.col
        surf.col[surf.col < 1] <- 1
        surf.col[surf.col > nCol] <- nCol
        if (is.na(col)) 
            col <- pal[as.array(surf.col)]
        z <- matrix(z, byrow=TRUE, n.grid, n.grid)
        if (plot.type == "contour") {
            stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("main" %in% 
                dnm, "", ",main=zlab"), ",...)", sep = "")
            if (color != "bw") {
                txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", stub, sep = "")
                eval(parse(text = txt))
                txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)", ifelse("add" %in% dnm, "", ",add=TRUE"), 
                  ",...)", sep = "")
                eval(parse(text = txt))
            } else {
                txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)", stub, sep = "")
                eval(parse(text = txt))
            }
            gfc <- getFigCoords('p')
            rasterImage(too.far.raster, xleft=gfc[1], xright=gfc[2], ybottom=gfc[3], ytop=gfc[4])
            if(add.color.legend){
                gradientLegend(c(min.z, max.z), n.seg=3, pos=.875, 
                    color=pal, dec=dec)
            }
	        if(hide.label==FALSE){
	        	addlabel = "fitted values"
	        	if(!is.null(rm.ranef)){
	        		if(rm.ranef !=FALSE){
	        			addlabel = paste(addlabel, "excl. random", sep=", ")
	        		}
	        	}
	        	mtext(addlabel, side=4, line=0, adj=0, 
	        		cex=.75, col='gray35', xpd=TRUE)
	        	if(!is.null(transform)){
	        		mtext("transformed", side=4, line=.75, adj=0, 
	        		cex=.75, col='gray35', xpd=TRUE)
	        	}
	        }
	            
        }else{
             stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("main" %in% 
                dnm, "", ",main=zlab"), ",...)", sep = "")
            if (color == "bw") {
                op <- par(bg = "white")
                txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", stub, sep = "")
                eval(parse(text = txt))
                par(op)
            }
            else {
                txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                  stub, sep = "")
                eval(parse(text = txt))
            } 
	        if(hide.label==FALSE){
	        	addlabel = "fitted values"
	        	if(!is.null(rm.ranef)){
	        		if(rm.ranef !=FALSE){
	        			addlabel = paste(addlabel, "excl. random", sep=", ")
	        		}
	        	}
	        	mtext(addlabel, side=4, line=0, adj=0, 
	        		cex=.75, col='gray35', xpd=TRUE)
	        	if(!is.null(transform)){
	        		mtext("transformed", side=4, line=.75, adj=0, 
	        		cex=.75, col='gray35', xpd=TRUE)
	        	}
	        }
        }
    } else {
        z.fit <- newd$fit
        z.cil <- newd$fit - newd$CI
        z.ciu <- newd$fit + newd$CI
        if(!is.null(transform)){
            z.fit <- sapply(z.fit, transform)
            z.cil <- sapply(z.cil, transform)
            z.ciu <- sapply(z.ciu, transform)
        }
        if (color == "bw" || color == "gray") {
            subs <- paste("grey are +/-", se, "s.e.")
            lo.col <- "gray"
            hi.col <- "gray"
        } else {
            subs <- paste("red/green are +/-", se, "s.e.")
            lo.col <- "green"
            hi.col <- "red"
        }
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
                stop("Something wrong with zlim")
            min.z <- zlim[1]
            max.z <- zlim[2]
        } else {
            z.max <- max(z.ciu, na.rm = TRUE)
            z.min <- min(z.cil, na.rm = TRUE)
        }
        zlim <- c(z.min, z.max)
        z <- matrix(z.cil, byrow=TRUE, n.grid, n.grid)
        if (plot.type == "contour") 
            warning("sorry no option for contouring with errors: try plot.gam")
        stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
            dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, "", ",sub=subs"), ",...)", sep = "")
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=lo.col"), stub, sep = "")
        eval(parse(text = txt))
        par(new = TRUE)
        z <- matrix(z.fit, byrow=TRUE, n.grid, n.grid)
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=\"black\""), stub, sep = "")
        eval(parse(text = txt))
        par(new = TRUE)
        z <- matrix(z.ciu, byrow=TRUE, n.grid, n.grid)
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=hi.col"), stub, sep = "")
        eval(parse(text = txt))
        if(hide.label==FALSE){
        	addlabel = "fitted values"
        	if(!is.null(rm.ranef)){
        		if(rm.ranef !=FALSE){
        			addlabel = paste(addlabel, "excl. random", sep=", ")
        		}
        	}
        	mtext(addlabel, side=4, line=0, adj=0, 
        		cex=.75, col='gray35', xpd=TRUE)
        	if(!is.null(transform)){
        		mtext("transformed", side=4, line=.75, adj=0, 
        		cex=.75, col='gray35', xpd=TRUE)
        	}
        }
    }
    invisible(list(fv = newd, m1 = m1, m2 = m2, zlim=c(min.z, max.z), too.far = ex.tf,
         note=paste(ifelse(is.null(transform), "type=lpmatrix, not on response scale", transform),
         	sprintf("Model predictor %s is transformed back to %s and %s", pca.term, view[1], view[2]),
         	sep=";")) )   
}
 





