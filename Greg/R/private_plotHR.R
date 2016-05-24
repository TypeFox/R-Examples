#' Gets the non-linear function's estimate
#' 
#' The function uses predict if not specified contrast
#' in order to attain the estimate, upper, and lower 
#' confidence interval
#' 
#' @param model The fit of the model to be plotted
#' @param term.label The name of the label
#' @param ylog If the outcome should be presented in the anti-log form, i.e.
#'  \code{exp()}.
#' @param cntrst A boolean that indicates if the \code{\link[rms]{contrast}}
#'  function is to be deployed for \pkg{rms} generated functions. The nice
#'  thing is that you get the median as a reference by default.
#' @param xlim The xlim if provided
#' @param new_data If not provided the function looks for the most common values
#'  i.e. median for continuous and mode for factors.
#' @importFrom rms contrast
#' @return \code{data.frame} with the columns xvalues, fit, ucl, lcl
#' @keywords internal
prPhEstimate <- function(model, 
                         term.label,
                         ylog,
                         cntrst,
                         xlim,
                         alpha,
                         new_data){
  if (missing(new_data)){
    new_data <- 
      prPhNewData(model = model, 
                  term.label = term.label, 
                  xlim = xlim)
  }else{
    # Unfortunately we cannot rely on all the
    # needed variables to be present among the
    # supplied data-frame and we need to see if 
    # we're missing something
    tmp <- 
      prPhNewData(model = model, 
                  term.label = term.label, 
                  xlim = xlim)
    
    # Add these variables
    add_vars <- colnames(tmp)[!colnames(tmp) %in% colnames(new_data) &
                                 colnames(tmp) != term.label]
    for (dv in add_vars){
      # Copy the element to the previous new_data
      # Note that the value is the same for the 
      # entire column and hence we only need the first value
      new_data[[dv]] <- tmp[1, dv]
    }

    # Remove these variables
    rm_vars <- colnames(new_data)[!colnames(new_data) %in% colnames(tmp) &
                                colnames(new_data) != term.label]
    for (r in rm_vars){
      new_data[[r]] <- NULL
    }
  }
  
  if(inherits(model, "cph")){
    if (cntrst){
      a <- list(new_data[,term.label])
      b <- list(attr(new_data, "median"))
      names(b) <- names(a) <- term.label

      for (n in names(new_data)){
        if (!n %in% names(a)){
          a[[n]] <- new_data[1,n]
          b[[n]] <- new_data[1,n]
        }
      }
      cntr <- contrast(model,
                       a=a,
                       b=b)
      df <- as.data.frame(cntr[c("Contrast", 
                                 "Lower",
                                 "Upper")])
    }else{
      pred <- predict (model, 
                       newdata=new_data, 
                       conf.int=1-alpha,
                       expand.na=FALSE, 
                       na.action=na.pass)
      df <- as.data.frame(pred[c("linear.predictors", "lower", "upper")])
    }
    
  }else{
    if (cntrst){
      stop("Contrast plotting is not defined for the models of class '", 
           paste(class(model), collapse="', '"), "'")
    }else{
      alt_label <- term.label
      if (!alt_label %in% names(model$assign)){
        # Assume that the term is contained within a function call
        # such as pspline() or similar
        tmp <- names(model$assign)[grep(alt_label, names(model$assign), fixed=TRUE)]
        if (length(tmp) != 1){
          stop("Could not identify the term",
               " '", alt_label,"'",
               " among the following model terms:",
               " '", paste(names(model$assign), collapse="', '"),"'")
        }
        alt_label <- tmp
        rm(tmp)
      }
      pred <- predict(model, newdata=new_data, type="terms", 
                      se.fit = TRUE, terms = alt_label)
      pred$upper <- as.double(pred$fit + qnorm(1 -alpha/2) * pred$se.fit)
      pred$lower <- as.double(pred$fit - qnorm(1 -alpha/2) * pred$se.fit)
      df <- as.data.frame(pred[c("fit", "lower", "upper")])
    }
  }

  colnames(df) <- c("estimate",
                    "lower",
                    "upper")
  
  df <- cbind(xvalues= new_data[,term.label],
              df)
  
  # Change to exponential form
  if (ylog == FALSE){
    for (n in names(df)){
      if (n != "xvalues")
        df[,n] <- exp(df[,n])
    }
  }
    
  attr(df, "new_data") <- new_data
  return(df)    
}

#' A function for retrieving new_data argument for predict
#' 
#' @param model The model fit from \code{\link[survival]{coxph}}
#'  or \code{\link[rms]{cph}}
#' @param term.label The label that is the one that \code{\link{plotHR}}
#'  intends to plot.
#' @param xlim The x-limits for the plot if any
#' @return \code{data.frame}
#' @keywords internal
prPhNewData <- function(model, term.label, xlim){
  # Get new data to use as basis for the prediction
  new_data <- prGetModelData(model)
  
  # Remove any Surv class variable as this is not 
  # part of the predictors
  new_data <-
    new_data[,sapply(new_data, function(x) !inherits(x, "Surv")),
             drop=FALSE]
  
  getMode <- function(x){
    tbl <- table(x)
    ret <- names(tbl)[which.max(tbl)]
    factor(ret, levels=names(tbl))
  }
  
  # Set all other but the variable of interest to the
  # mode or the median
  for (variable in colnames(new_data)){
    if (variable != term.label){
      if (is.numeric(new_data[,variable])){
        new_data[, variable] <- median(new_data[, variable])
      }else{
        new_data[, variable] <- getMode(new_data[, variable])
      }
    }
  }
  if (is.numeric(new_data[, term.label])){
    nd_median <- median(new_data[, term.label], na.rm=TRUE)
  }else{
    nd_median <- getMode(new_data[, term.label])
  }
  
  new_data <- new_data[!duplicated(new_data[, term.label]), ,drop=FALSE]
  
  if (!missing(xlim)){
    if (NCOL(new_data) == 1){
      new_data <- as.data.frame(
        matrix(
          new_data[new_data >= min(xlim) & 
                     new_data <= max(xlim)], 
          ncol=1, 
          dimnames=list(NULL, c(term.label))))
      
    }else{
      new_data <- new_data[new_data[, term.label] >= min(xlim) &
                             new_data[, term.label] <= max(xlim), ,
                           drop=FALSE]
      
    }
  }
  
  attr(new_data, "median") <- nd_median
  return(new_data)
}


#' Plots the confidence intervals
#' 
#' Uses \code{\link[graphics]{polygon}} or 
#' \code{\link[graphics]{lines}} to plot confidence
#' intervals.
#' 
#' @param model_data A data frame with 'xvalues', 'upper', and 'lower'
#'  columns.
#' @param color The color of the line/polygon
#' @param polygon Boolean indicating polygon or line
#' @param lwd Line width - see \code{\link[grid]{gpar}}
#' @param lty Line type - see \code{\link[grid]{gpar}}
#' @keywords internal
#' @return \code{void} The function performs the print
prPhConfIntPlot <- function(model_data, color, polygon, lwd, lty){
  current_i.backw <- order(model_data$xvalues , decreasing = TRUE)
  current_i.forw <- order(model_data$xvalues)
  
  if (polygon){
    # The x-axel is always the same
    x.poly <- c(model_data$xvalues[current_i.forw] , model_data$xvalues[current_i.backw])
    # The y axel is based upin the current model
    y.poly <- c(model_data$upper[current_i.forw] , model_data$lower[current_i.backw])
    polygon(x.poly , y.poly , col = color, border = NA)
  }else{
    lines(model_data$xvalues[current_i.forw] , model_data$upper[current_i.forw], 
          col = color, 
          lwd = lwd, lty= lty)
    lines(model_data$xvalues[current_i.forw], 
          model_data$lower[current_i.forw], 
          col = color, 
          lwd = lwd, lty= lty)
  }
}

#' Plot a rug on the datapoints
#' 
#' @param xvalues The xvalues that are used for the density
#' @return \code{void}
#' @keywords internal
prPhRugPlot <- function(xvalues){
  # rugs at datapoints
  axis(side = 1, 
       line = -1.2, 
       at = jitter(xvalues), 
       labels = FALSE, 
       tick = TRUE, 
       tcl = 0.8, 
       lwd.ticks = 0.1, 
       lwd = 0)
  
  # rugs and labels at 1Q, median and 3Q
  axis(side = 1, 
       line = -1.0, 
       at = fivenum(xvalues)[2:4], 
       lwd = 0, 
       tick = TRUE, 
       tcl = 1.2, 
       lwd.ticks = 1, 
       col.ticks = "black", 
       labels = c("Quartile 1","Median","Quartile 3"), 
       cex.axis = 0.7, 
       col.axis = "black", 
       padj = -2.8)
  axis(side = 1,
       line = 0.0, 
       at = fivenum(xvalues)[2:4], 
       lwd = 0, 
       tick = TRUE, 
       tcl = 0.2, 
       lwd.ticks = 1, 
       col.ticks = "black", 
       labels = FALSE)
}

#' Plot a density on the datapoints
#' 
#' @param xvalues The xvalues that are used for the density
#' @param color The color of the density polygon
#' @return \code{void}
#' @keywords internal
prPhDensityPlot <- function(xvalues, color){
  # calculate the coordinates of the density function
  density <- density( xvalues )
  # the height of the densityity curve
  max.density <- max(density$y)
  
  # Get the boundaries of the plot to
  # put the density polygon at the x-line
  plot_coordinates <- par("usr")
  
  # get the "length" and range of the y-axis
  y.scale <- plot_coordinates[4] - plot_coordinates[3]
  
  # transform the y-coordinates of the density
  # to the lower 10% of the plotting panel
  density$y <- (0.1 * y.scale / max.density) * density$y + plot_coordinates[3]
  
  ## plot the polygon
  polygon(density$x , density$y, 
          border = FALSE, col = color)
}
