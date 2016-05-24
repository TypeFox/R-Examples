################################################################################
# 
# plot.predict.downscale.R
# Version 1.1
# 13/02/2015
#
# Updates:
#   13/03/2015: if 0's predicted don't plot them
#
# Plot the observed and predicted area of occupancy against grain size on
# log-log axes.
#
# Args:
#   predict.object: an object of class 'predict.downscale' containing observed
#                   and predicted data 
#   ...: arguments, including graphical parameters, passed to other methods.
#
# Returns:
#   no object returned.
#
################################################################################

plot.predict.downscale <- function(x,
                                   xlim = NULL,
                                   ylim  = NULL,
                                   xlab = NULL,
                                   ylab = NULL,
                                   main = NULL,
                                   lwd.obs = NULL,
                                   lwd.pred = NULL,
                                   col.obs = NULL,
                                   col.pred = NULL,
                                   ...) {
  predict.object <- x
  # error checking
  if (class(predict.object) != "predict.downscale"){
    stop("Input data not of class 'predict.downscale'")
  }
  
  model.run <- predict.object$model
  observed <- predict.object$observed
  predicted <- predict.object$predicted
  predicted[predicted == 0] <- NA

  ### plotting pars
  if(is.null(xlim)) {
    xlim <- c(min(c(observed[, "Cell.area"], 
                    predicted[, "Cell.area"]), na.rm = TRUE),
              max(c(observed[, "Cell.area"], 
                    predicted[, "Cell.area"]), na.rm = TRUE))
  }
  if(is.null(ylim)) {
    ylim <- c(min(c(observed[, "Occupancy"], 
                    predicted[, "Occupancy"]), na.rm = TRUE), 1)
  }
  if(is.null(xlab)) {
    xlab <- "Log cell area"
  }
  if(is.null(ylab)) {
    ylab <- "Log occupancy"
  }
  if(is.null(main)) {
    main <- paste(model.run, "model")
  }
  if(is.null(col.obs)) {
    col.obs <- "black"
  }
  if(is.null(col.pred)) {
    col.pred <- "red"
  }
  if(is.null(lwd.obs)) {
    lwd.obs <- 2
  }
  if(is.null(lwd.pred)) {
    lwd.pred <- 2
  }

  plot(observed[, "Occupancy"] ~ observed[, "Cell.area"],
       type = "n",
       log = "xy",
       xlim = xlim,
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       main = main,
       ...)
  
  points(observed[, "Occupancy"] ~ observed[, "Cell.area"],
         type = "b",
         col = col.obs,
         lwd = lwd.obs,
         ...)
  
  points(predicted[, "Occupancy"] ~ predicted[, "Cell.area"],
         type = "b",
         col = col.pred,
         lwd = lwd.pred,
         ...)
}
