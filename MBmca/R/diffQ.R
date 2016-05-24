#' Calculation of the melting temperature (Tm) from the first derivative
#' 
#' \code{diffQ} is used to calculate the melting temperature (Tm) but also for
#' elementary graphical operations (e.g., show the Tm or the derivative). It
#' does not require smoothed data for the MCA. The parameter \code{rsm} can be
#' used to double the temperature resolution by calculation of the mean
#' temperature and mean fluorescence. Note: mcaSmoother has the \code{n}
#' parameter with a similar functionality. First the approximate Tm is
#' determined as the \code{min()} and/or \code{max()} from the first
#' derivative. The first numeric derivative (Forward Difference) is estimated
#' from the values of the function values obtained during an experiment since
#' the exact function of the melting curve is unknown. The method used in
#' \code{diffQ} is suitable for independent variables that are equally and
#' unequally spaced. Alternatives for the numerical differentiation include
#' Backward Differences, Central Differences or Three-Point (Forward or
#' Backward) Difference based on Lagrange Estimation (currently not implemented
#' in \code{diffQ}). The approximate peak value is the starting-point for a
#' function based calculation. The function takes a defined number n (maximum
#' 8) of the left and the right neighbor values and fits a quadratic
#' polynomial. The quadratic regression of the X (temperature) against the Y
#' (fluorescence) range gives the coefficients. The optimal quadratic
#' polynomial is chosen based on the highest adjusted R-squared value.
#' \code{diffQ} returns an objects of the class \code{list}. To accessing
#' components of lists is done as described elsewhere either by name or by
#' number. \code{diffQ} has a simple plot function. However, for sophisticated
#' analysis and plots its recommended to use \code{diffQ} as presented in the
#' examples as part of algorithms.
#' 
#' 
#' @param xy is a \code{data.frame} containing in the first column the
#' temperature and in the second column the fluorescence values. Preferably the
#' output from \code{mcaSmoother} is used.
#' @param fct accepts \code{min} or \code{max} as option and is used to define
#' whether to find a local minimum (``negative peak'') or local maximum
#' (``positive peak'').
#' @param fws defines the number (n) of left and right neighbors to use for the
#' calculation of the quadratic polynomial.
#' @param col is a graphical parameter used to define the length of the line
#' used in the plot.
#' @param plot shows a plot of a single melting curve. To draw multiple curves
#' in a single plot set \code{plot = FALSE} and create and empty plot instead
#' (see examples).
#' @param verbose shows additional information (e.g., approximate derivative,
#' ranges used for calculation, approximate Tm) of the calculation.
#' @param warn diffQ tries to keep the user as informed as possible about the
#' quality of the analysis. However, in some scenarios are the warning and
#' message about analysis not needed or disturbing.  \code{warn} can be used to
#' stop the flodding of the output.
#' @param peak shows the peak in the plot (see examples).
#' @param negderiv uses the positive first derivative instead of the negative.
#' @param deriv shows the first derivative with the color assigned to
#' \code{col} (see examples).
#' @param derivlimits shows the neighbors (fws) used to calculate the Tm as
#' points in the plot (see examples).
#' @param derivlimitsline shows the neighbors (fws) used to calculate the Tm as
#' line in the plot (see examples).
#' @param vertiline draws a vertical line at the Tms (see examples).
#' @param rsm performs a doubling of the temperature resolution by calculation
#' of the mean temperature and mean fluorescence between successive temperature
#' steps. Note: \code{mcaSmoother} has the "n" parameter with a similar but
#' advanced functionality.
#' @param inder Interpolates first derivatives using the five-point stencil.
#' See \code{chipPCR} package for details.
#' @return \item{diffQ() }{returns a comprehensive list (if parameter verbose
#' is TRUE) with results from the first derivative. The list includes a
#' \code{data.frame} of the derivative ("xy"). The temperature range
#' ("limits.xQ") and fluorescence range ("limits.diffQ") to calculate the peak
#' value. "fluo.x" is the approximate fluorescence at the approximate melting
#' temperature. The calculated melting temperature ("Tm") with the
#' corresponding fluorescence intensity ("fluoTm"). The number of neighbors
#' ("fws"), the adjusted R-squared ("adj.r.squared") and the
#' normalized-root-mean-squared-error ("NRMSE") to fit. The quality of the
#' calculated melting temperature ("Tm") can be checked with devsum which
#' reports the relative deviation (in percent) between the approximate melting
#' temperature and the calculated melting temperature, if NRMSE is less than
#' 0.08 and the adjusted R-squared is less than 0.85. A relative deviation
#' larger than 10 percent will result in a warning. Reducing fws might improve
#' the result. }
#' 
#' \item{Tm }{returns the calculated melting temperature ("Tm").}
#' 
#' \item{fluoTm }{returns the calculated fluorescence at the calculated melting
#' temperature.}
#' 
#' \item{Tm.approx }{returns the approximate melting temperature.}
#' 
#' \item{fluo.x }{returns the approximate fluorescence at the calculated
#' melting temperature.}
#' 
#' \item{xy }{returns the approximate derivative value used for the calculation
#' of the melting peak.}
#' 
#' \item{limits.xQ }{returns a data range of temperature values used to
#' calculate the melting temperature.}
#' 
#' \item{limits.diffQ }{returns a data range of fluorescence values used to
#' calculate the melting temperature.}
#' 
#' \item{adj.r.squared }{returns the adjusted R-squared from the quadratic
#' model fitting function (see also \code{fit}).}
#' 
#' \item{NRMSE }{returns the normalized root-mean-squared-error (NRMSE) from
#' the quadratic model fitting function (see also \code{fit}).}
#' 
#' \item{fws }{returns the number of points used for the calculation of the
#' melting temperature.}
#' 
#' \item{devsum }{returns measures to show the difference between the
#' approximate and calculated melting temperature.}
#' 
#' \item{temperature }{returns measures to investigate the temperature
#' resolution of the melting curve. Raw fluorescence measurements at irregular
#' temperature resolutions (intervals) can introduce artifacts and thus lead to
#' wrong melting point estimations.}
#' 
#' \item{temperature$T.delta }{returns the difference between two successive
#' temperature steps.}
#' 
#' \item{temperature$mean.T.delta }{returns the mean difference between two
#' temperature steps.}
#' 
#' \item{temperature$sd.T.delta }{returns the standard deviation of the
#' temperature.}
#' 
#' \item{temperature$RSD.T.delta }{returns the relative standard deviation
#' (RSD) of the temperature in percent.}
#' 
#' \item{fit }{returns the summary of the results of the quadratic model
#' fitting function.}
#' @author Stefan Roediger
#' @seealso \code{\link{diffQ2}}, \code{\link{mcaSmoother}}
#' @references A Highly Versatile Microscope Imaging Technology Platform for
#' the Multiplex Real-Time Detection of Biomolecules and Autoimmune Antibodies.
#' S. Roediger, P. Schierack, A. Boehm, J. Nitschke, I. Berger, U. Froemmel, C.
#' Schmidt, M. Ruhland, I. Schimke, D. Roggenbuck, W. Lehmann and C.
#' Schroeder.  \emph{Advances in Biochemical Bioengineering/Biotechnology}.
#' 133:33--74, 2013. \url{http://www.ncbi.nlm.nih.gov/pubmed/22437246}
#' 
#' Nucleic acid detection based on the use of microbeads: a review. S.
#' Roediger, C. Liebsch, C. Schmidt, W. Lehmann, U. Resch-Genger, U. Schedler,
#' P. Schierack. \emph{Microchim Acta} 2014:1--18. DOI:
#' 10.1007/s00604-014-1243-4
#' 
#' Roediger S, Boehm A, Schimke I. Surface Melting Curve Analysis with R.
#' \emph{The R Journal} 2013;5:37--53.
#' @keywords melting Tm
#' @examples
#' 
#' # First Example
#' # Plot the first derivative of different samples for single melting curve
#' # data. Note that the argument "plot" is TRUE.
#' 
#' data(MultiMelt)
#' par(mfrow = c(1,2))
#' sapply(2L:14, function(i) {
#'         tmp <- mcaSmoother(MultiMelt[, 1], MultiMelt[, i])
#'         diffQ(tmp, plot = TRUE)
#'   }
#' )
#' par(mfrow = c(1,1))
#' # Second example
#' # Plot the first derivative of different samples from MultiMelt
#' # in a single plot.
#' data(MultiMelt)
#' 
#' # First create an empty plot
#' plot(NA, NA, xlab = "Temperature", ylab ="-d(refMFI)/d(T)",
#'         main = "Multiple melting peaks in a single plot", xlim = c(65,85),
#'         ylim = c(-0.4,0.01), pch = 19, cex = 1.8)
#' # Prepossess the selected melting curve data (2,6,12) with mcaSmoother 
#' # and apply them to diffQ. Note that the argument "plot" is FALSE
#' # while other arguments like derivlimitsline or peak are TRUE. 
#' sapply(c(2,6,12), function(i) {
#' 	tmp <- mcaSmoother(MultiMelt[, 1], MultiMelt[, i], 
#' 			    bg = c(41,61), bgadj = TRUE)
#' 	diffQ(tmp, plot = FALSE, derivlimitsline = TRUE, deriv = TRUE, 
#' 	      peak = TRUE, derivlimits = TRUE, col = i, vertiline = TRUE)
#'   }
#' )
#' legend(65, -0.1, colnames(MultiMelt[, c(2,6,12)]), pch = c(15,15,15), 
#' 	col = c(2,6,12))
#' 
#' # Third example
#' # First create an empty plot
#' plot(NA, NA, xlim = c(50,85), ylim = c(-0.4,2.5), 
#'      xlab = "Temperature", 
#'      ylab ="-refMFI(T) | refMFI'(T) | refMFI''(T)",
#'      main = "1st and 2nd Derivatives", 
#'      pch = 19, cex = 1.8)
#' 
#' # Prepossess the selected melting curve data with mcaSmoother 
#' # and apply them to diffQ and diffQ2. Note that 
#' # the argument "plot" is FALSE while other 
#' # arguments like derivlimitsline or peak are TRUE.
#' 
#' tmp <- mcaSmoother(MultiMelt[, 1], MultiMelt[, 2], 
#' 		    bg = c(41,61), bgadj = TRUE)
#' lines(tmp, col= 1, lwd = 2)
#' 
#' # Note the different use of the argument derivlimits in diffQ and diffQ2
#' diffQ(tmp, fct = min, derivlimitsline = TRUE, deriv = TRUE, 
#' 	    peak = TRUE, derivlimits = FALSE, col = 2, vertiline = TRUE)
#' diffQ2(tmp, fct = min, derivlimitsline = TRUE, deriv = TRUE, 
#' 	    peak = TRUE, derivlimits = TRUE, col = 3, vertiline = TRUE)
#' 
#' # Add a legend to the plot
#' legend(65, 1.5, c("Melting curve",
#' 		  "1st Derivative", 
#' 		  "2nd Derivative"), 
#' 		  pch = c(19,19,19), col = c(1,2,3))
#' 
#' # Fourth example
#' # Different curves may potentially have different quality in practice. 
#' # For example, using data from MultiMelt should return a 
#' # valid result and plot.
#' data(MultiMelt)
#' 
#' diffQ(cbind(MultiMelt[, 1], MultiMelt[, 2]), plot = TRUE)$Tm
#' # limits_xQ
#' #  77.88139
#' 
#' # Imagine an experiment that went terribly wrong. You would 
#' # still get an estimate for the Tm. The output from diffQ, 
#' # with an error attached, lets the user know that this Tm 
#' # is potentially meaningless. diffQ() will give several 
#' # warning messages.
#' 
#' set.seed(1)
#' y = rnorm(55,1.5,.8)
#' diffQ(cbind(MultiMelt[, 1],y), plot = TRUE)$Tm
#' 
#' # The distribution of the curve data indicates noise.
#' # The data should be visually inspected with a plot 
#' # (see examples of diffQ). The Tm calculation (fit, 
#' # adj. R squared ~ 0.157, NRMSE ~ 0.279) is not optimal 
#' # presumably due to noisy data. Check raw melting 
#' # curve (see examples of diffQ).
#' # Calculated Tm 
#' #      56.16755
#' 
#' 
#' # Sixth example
#' # Curves may potentially have a low temperature resolution. The rsm 
#' # parameter can be used to double the temperature resolution.
#' # Use data from MultiMelt column 15 (MLC2v2).
#' data(MultiMelt)
#' tmp <- cbind(MultiMelt[, 1], MultiMelt[, 15])
#' 
#' # Use diffQ without and with the rsm parameter and plot
#' # the results in a single row
#' par(mfrow = c(1,2))
#' 
#' diffQ(tmp, plot = TRUE)$Tm
#'   text(60, -0.15, "without rsm parameter")
#' 
#' diffQ(tmp, plot = TRUE, rsm = TRUE)$Tm
#'   text(60, -0.15, "with rsm parameter")
#' par(mfrow = c(1,1))
#' 
#' @export diffQ
"diffQ" <- function(xy, fct = min, fws = 8, col = 2, plot = FALSE, 
                    verbose =  FALSE, warn = TRUE, 
                    peak = FALSE, negderiv = TRUE, deriv = FALSE, 
                    derivlimits = FALSE, derivlimitsline = FALSE, 
                    vertiline = FALSE, rsm = FALSE, inder = FALSE) { 
  # Test if fws (number of neighbors) is within a meaningful range.
  old.warn <- options("warn")[["warn"]]
  options(warn = -1)
  fws <- round(fws)
  if (fws < 2 || fws > 8) 
    stop("Fit window size must be within 2 and 8.")
  
  list.res <- list()
  
  # Take the x and y values from the object of type data.frame.
  x <- xy[, 1]
  y <- xy[, 2]
  
  # Test if x and y exist.
  if (is.null(x)) 
    stop("Enter temperature")
  if (is.null(y)) 
    stop("Enter fluorescence data")
  
  # Determine the temperature resolution of the melting curve data
  deltaT <- vector()
  for (i in 1L:(length(x) - 1)){
    tmp <- abs(x[i] - x[i + 1])
    deltaT <- c(deltaT,tmp)
  }
  
  deltaT.mean 	<- mean(deltaT)
  deltaT.sd 	<- sd(deltaT)
  deltaT.RSD	<- deltaT.sd/deltaT.mean * 100
  temperature <- list(T.delta = deltaT, mean.T.delta = deltaT.mean, 
                      sd.T.delta = deltaT.sd, 
                      RSD.T.delta = deltaT.RSD)
  
  # Function to increase the resolution of the melting curve by 
  # calculation of the mean temperature and mean fluorescence
  if (rsm == TRUE) {
    ind.low <- seq(1, length(x) - 1, 1)
    ind.up <- seq(2, length(x), 1)
    xt <- (x[ind.up] + x[ind.low]) / 2
    yt <- (y[ind.up] + y[ind.low]) / 2
    xy <- data.frame(c(x, xt), c(y, yt))
    xy <- xy[order (xy[, 1]), ]
    x <- xy[, 1]
    y <- xy[, 2]
  }
  
  # Function to calculate the approximate derivative.
  if (negderiv) {y <- y * -1} 	# Change sign of curve data
  
  N <- length(x); low <- 1:2; up <- (N - 1):N
  
  if (inder) {	inder.tmp <- inder(x, y, Nip = 1)
               xQ <- inder.tmp[, "x"]
               diffQ <- inder.tmp[, "d1y"]
               
  } else {  xQ <- x[3L:length(x) - 1] 
            h <- (x[-low] - x[-up]) 	# Step size
            diffQ <- (y[-low] - y[-up]) / h
  }
  
  out <- data.frame(xQ, diffQ)
  names(out) <- c("Temperature", "d(F) / dT")
  
  # First: Find approximate range of neighbors for minimum or 
  # maximum of temperature peak.Second: Calculate quadratic 
  # polynomial for ranges of minimum or maximum temperature peaks.
  # Return the results in an object of the type list.
  for (i in 2L:9) {
    suppressMessages(RANGE <- which(diffQ == fct(diffQ[(i + 1):(N - i)]))) 
    if (class(try(na.omit(xQ[(RANGE - i):(RANGE + i)]), 
                  silent = T)) == "try-error") {
      list.res[[i-1]] <- NA
    } else {
      limits.xQ <- na.omit(xQ[(RANGE - i):(RANGE + i)])
      limits.diffQ <- na.omit(diffQ[(RANGE - i):(RANGE + i)])
      fluo.x <- limits.diffQ[length(limits.diffQ) - i]
      lm2 <- lm(limits.diffQ ~ limits.xQ + I(limits.xQ^2))
      coeflm2 <- data.frame(lm2[1])[c(1:3), ]
      coeflm2.y <- coeflm2[1] + coeflm2[2] * limits.xQ + coeflm2[3] * limits.xQ^2
      # lm2sum <- summary(lm(coeflm2.y ~ limits.diffQ))
      lm2sum <- summary(lm2)
      list.res[[i-1]] <- list(i, limits.xQ, limits.diffQ, fluo.x, 
                              lm2, coeflm2, coeflm2.y, lm2sum
      )
    }
  }
  
  # Determine the optimal fitted quadratic polynomial for ranges 
  # of minimum or maximum temperature peaksbases on the adjusted 
  # R squared.
  Rsq <- matrix(NA, nrow = fws, ncol = 2)
  colnames(Rsq) <- c("fw", "Rsqr")
  for (i in 1L:fws) {
    Rsq[i, 1] <- list.res[[i]][[1]]
    Rsq[i, 2] <- list.res[[i]][[8]][["adj.r.squared"]]
  }
  list.res <- list.res[[which(Rsq[, 2] == max(na.omit(Rsq[, 2])))]]
  names(list.res) <- list("fw", "limits.xQ", "limits.diffQ", 
                          "fluo.x", "lm2", "coeflm2", 
                          "coeflm2.y", "lm2sum")
  
  limits.xQ <- list.res[["limits.xQ"]]
  limits.diffQ <- list.res[["limits.diffQ"]]
  
  fluo.x <- list.res[["fluo.x"]]
  names(fluo.x) <- c("Signal height at approximate Tm")
  
  lm2 <- list.res[["lm2"]]
  coeflm2 <- list.res[["coeflm2"]]
  coeflm2.y <- list.res[["coeflm2.y"]]
  lm2sum <- list.res[["lm2sum"]]
  # fw <- list.res$fw
  
  # Polynom to fit the area of the calculated Tm
  poly.fct <- function(xi) coeflm2[1] + coeflm2[2] * xi + coeflm2[3] * xi^2
  # Calculate the Tm and assign meaningful names to variables
  abl <- -lm2[["coefficients"]][2] / (2 * lm2[["coefficients"]][3])
  names(abl) <- "Calculated Tm"
  y <-	coeflm2[1] + coeflm2[2] * abl + coeflm2[3] * abl^2
  names(y) <- c("Signal height at calculated Tm")
  
  # Optional draw line of calculated approximate derivative.
  if ((!plot) && (.Device != "null device")) {
    if (vertiline) {
      abline(v = abl, col = "grey")
    }
    if (deriv) {
      lines(xQ, diffQ, col = col)  
    } 
    if (derivlimits) {
      points(limits.xQ, limits.diffQ, cex = 1, pch = 19, col = col)
    }
    if (derivlimitsline) {
      lines(spline(limits.xQ[1L:length(limits.xQ)], 
                   lm2[["fitted.values"]][1L:length(lm2[["fitted.values"]])]), 
            col = "orange", lwd = 2)
    }
    if (peak) {
      points(abl, y, cex = 1, pch = 19, col = col)
    }
  }
  
  # Test if the calculated Tm and the approximate Tm differ strongly
  dev.dat <- data.frame(limits.xQ, limits.diffQ)
  dev <- dev.dat[which(dev.dat[, 2] == max(dev.dat[, 2])), ]
  dev.var <- sd(c(dev[1, 1], abl)) / mean(c(dev[1, 1], abl)) * 100
  dev.sum <- data.frame(dev.var, dev[1, 1], abl)
  colnames(dev.sum) <- c("Relative Deviation (%)", "Approximate Tm", 
                         "Calculated Tm")
  rownames(dev.sum) <- NULL
  
  warn.approx.calc <- c("Approximate and calculated Tm varri. This is an expected behaviour \n
			 but the calculation should be confirmed with a plot (see examples of diffQ).")
  
  if (warn && dev.sum[1] > 5) {
    message(warn.approx.calc)
    #TO DO: maybe incorporate print into message?
    message(dev.sum)
  } else {warn.approx.calc = "ok"}
  
  # Calculates the Root Mean Squared Error
  NRMSE <- function(model = model, mes = mes) {
    RMSE <- sqrt(mean(residuals(model)^2))
    NRMSE <- RMSE / (max(mes) - min(mes))
    NRMSE.warning <- ifelse(NRMSE > 0.08, "NRMSE bad", "NRMSE ok")
    list(NRMSE = NRMSE, RMSE = RMSE, 
         NRMSE.warning = NRMSE.warning)
  }
  
  NRMSE.res <- NRMSE(model = lm2, mes = limits.diffQ)
  
  message.shapiro.test <- c("The distribution of the curve data indicates noise. The data should be visually \ninspected with a plot (see examples of diffQ).")
  message.polynomial.fit <- paste0("The Tm calculation (fit, adj. R squared ~ ", 
                   round(max(na.omit(Rsq[, 2])), 3), 
                   ", NRMSE ~ ", round(NRMSE.res[["NRMSE"]], 3), 
                   ") is not optimal presumably\ndue to noisy data. Check raw melting curve (see examples of diffQ).")
  
  # Simple test if data come from noise or presumably a melting curve
  if (warn && shapiro.test(xy[, 2])[["p.value"]] >= 10e-8) {
    message(message.shapiro.test)
    } else {message.shapiro.test = "ok"}
  
  # Simple test if polynomial fit performed accaptable
  if (warn && (max(na.omit(Rsq[, 2])) < 0.85)) {
    message(message.polynomial.fit)
    } else {message.polynomial.fit <- "ok"}
  
  if (plot) {
    plot(xQ, diffQ, xlab = "Temperature", ylab = "-d(F) / dT", 
         type = "b", col = col)
    points(limits.xQ, limits.diffQ, col = "orange", pch = 19)
    curve(poly.fct, limits.xQ[1], limits.xQ[length(limits.xQ)], 
          col = col, add = TRUE)
    points(abl, y, pch = 19, col = 2)
    if (vertiline) {
      abline(v = abl, col = "grey")
    }
    if (derivlimits) {
      points(limits.xQ, limits.diffQ, cex = 1, pch = 19, col = col)
    }
    if (derivlimitsline) {
      lines(spline(limits.xQ[1L:length(limits.xQ)], 
                   lm2[["fitted.values"]][1L:length(lm2[["fitted.values"]])]), 
            col = "orange", lwd = 2)
    }
  }
  
  #restore old warning value
  options(warn = old.warn)
  
  # Returns an object of the type list containing the data and 
  # data.frames from above including the approximate difference 
  # quotient values, melting temperatures, intensities and used 
  #neighbors.
  if (verbose) {
    res <- list(Tm = abl, fluoTm = y, 
         Tm.approx = dev[1], fluo.x = fluo.x, 
         xy = out, limits.xQ = limits.xQ, 
         limits.diffQ = limits.diffQ,
         adj.r.squared = lm2sum$adj.r.squared, 
         NRMSE = NRMSE.res$NRMSE,
         fws = list.res$fw, devsum=dev.sum, 
         temperature = temperature, 
         fit = summary(list.res$lm2),
         approx.calc = warn.approx.calc,
         shapiro.test = message.shapiro.test,
         polynomial.fit = message.polynomial.fit
         )
  } else {
    res <- list(Tm = abl, fluoTm = y)
  }
  class(res) <- "diffQobject"
  res
}
