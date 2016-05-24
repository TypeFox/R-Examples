#' Multiple comparison of the temperature dependent variance of the refMFI
#' 
#' MFIerror is used for a fast multiple comparison of the temperature dependent
#' variance of the refMFI. MFIerror returns an object of the class data.frame
#' with columns ``Temperature'', ``Location'' (Mean, Median), ``Deviation''
#' (Standard Deviation, Median Absolute Deviation) and ``Coefficient of
#' Variation''.
#' 
#' 
#' @param x is the column of a data frame for the temperature.
#' @param y are multiple columns of fluorescence values from a
#' \code{data.frame} (e.g., [, c(1:n)]).
#' @param CV If \code{CV} is true the coefficient of variation (RSD, CV) is
#' plotted. If set to FLASE the deviation as Standard Deviation or Median
#' Absolute Deviation is plotted.
#' @param RSD Setting the option \code{RSD=TRUE} shows the relative standard
#' deviation (RSD) in percent.
#' @param rob Using the option \code{rob} as TRUE the median and the median
#' absolute deviation (MAD) is plotted instead of the mean and standard
#' deviation.
#' @param errplot sets \code{MFIerror()} to plot the results (default). In the
#' default setting (\code{CV=FALSE}) the mean with the standard deviations is
#' plotted.
#' @param type is a graphical parameter setting the plot use lines, points or
#' both (see \code{\link{plot}}).
#' @param pch is a graphical parameter used to define the symbol used in the
#' plot.
#' @param length \code{length} is a graphical parameter used to define the
#' length of the error bar used in the plot.
#' @param col \code{col} is a graphical parameter used to define the color of
#' the error bar used in the plot.
#' @return \item{res }{ returns a \code{data.frame} containing the
#' "Temperature", "Location" (mean, median), "Deviation" (standard deviation,
#' median absolute deviation), "Coefficient of Variance" (CV, RSD) sequential
#' in the columns.}
#' @author Stefan Roediger
#' @seealso \code{\link{mcaSmoother}}
#' @keywords deviation
#' @examples
#' 
#' # First Example
#' # Temperature dependent variance of the refMFI using standard measures 
#' # (Mean, Standard Deviation (SD)).
#' # Use Standard Deviation (SD) in the plot
#' 
#' data(MultiMelt)
#' MFIerror(MultiMelt[, 1], MultiMelt[, c(2L:13)])
#' 
#' # Second Example
#' # Temperature dependent relative variance of the refMFI using robust 
#' # measures (Median, Median Absolute Deviation (MAD)). The parameter 
#' # errplot is set to FALSE in order to prevent the plot of the 
#' # coefficient of variation versus the temperature.
#' 
#' MFIerror(MultiMelt[, 1], MultiMelt[, c(2L:13)], errplot = FALSE, 
#' 	  RSD = TRUE, rob = TRUE)
#' 
#' # Third Example
#' # Temperature dependent relative variance of the refMFI using 
#' # robust measures (Median, Median Absolute Deviation (MAD)).
#' MFIerror(MultiMelt[, 1], MultiMelt[, c(2L:13)], RSD = TRUE, 
#'   rob = TRUE)
#' 
#' @export MFIerror
MFIerror <- function(x, y, CV = FALSE, RSD = FALSE, rob = FALSE, 
		       errplot = TRUE, type = "p", pch = 19, length = 0.05, 
		       col = "black") {
  #Define if "robust" or standard function should be used as measures
  old.warn <- options("warn")[[1]]
  options(warn = -1)
  #Test if x and y exist.
  if (is.null(x)) 
      stop("Enter temperature")
  if (is.null(y)) 
      stop("Enter fluorescence data")

  if (rob) {
      loc.fct <- median
      dev.fct <- mad
  } else {
      loc.fct <- mean
      dev.fct <- sd
    }
	
  if (ncol(data.frame(y)) > 1) {
      y.m <- apply(y, 1, loc.fct)
      y.sd <- apply(y, 1, dev.fct)
  } else {
      y.m <- y
      y.sd <- rep(0, length(y))
  }

  if (RSD) {
      y.cv <- (y.sd / y.m) * 100
  } else {
      y.cv <- y.sd / y.m
    }
  
  res <- data.frame(x, y.m, y.sd, y.cv)

  if (rob == TRUE && RSD == FALSE) {
      names(res) <- c("Temperature", "Location (Median)", "Deviation (MAD)", 
		      "Coefficient of Variance (RSD [%])")
  }
  
  if (rob == FALSE && RSD == FALSE) {
      names(res) <- c("Temperature", "Location (Mean)", "Deviation (SD)", 
		      "Coefficient of Variance (RSD [%])")
  }
  
  if (rob == TRUE && RSD == TRUE) {
      names(res) <- c("Temperature", "Location (Median)", "Deviation (MAD)", 
		      "Coefficient of Variance (RSD)")
  }
  
  if (rob == FALSE && RSD == TRUE) {
      names(res) <- c("Temperature", "Location (Mean)", "Deviation (SD)", 
		      "Coefficient of Variance (RSD)")
  }

  #Plot the Coefficient of Variance
  if (errplot) {
    if (CV) {
	plot(res[, 1], res[,4], xlab = "Temperature", ylab = "CV", 
	     col = col, pch = pch)
  #Plot the location with error bars.
    } else {
	plot(res[, 1], res[, 2], ylim = c(min(res[, 2] - res[, 3]), 
	     max(res[, 2] + res[, 3])), xlab = "Temperature", 
	     ylab = "MFI", type = type)
		
	    arrows(res[, 1], res[, 2] + res[, 3], res[, 1], 
		   res[, 2] - res[, 3], angle = 90, code = 3, 
		   length = length, col = col)
      }
  }
  #restore old warning value
  options(warn = old.warn)
  
  # res is the an object of the type data.frame containing the 
  # temperature, location, deviation and coefficient of variance.
  res
}
