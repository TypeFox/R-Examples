#' Function to pre-process melting curve data.
#' 
#' The function \code{mcaSmoother()} is used for data pre-processing.
#' Measurements from experimental systems may occasionally include missing
#' values (NA). \code{mcaSmoother()} uses \code{approx()} to fill up NAs under
#' the assumption that all measurements were equidistant. The original data
#' remain unchanged and only the NAs are substituted. Following it calls
#' \code{smooth.spline()} to smooth the curve. Different strengths can be set
#' using the option \code{df.fact} (f default~0.95). Internally it takes the
#' degree of freedom value from the spline and multiplies it with a factor
#' between 0.6 and 1.1. Values lower than 1 result in stronger smoothed curves.
#' The outcome of the differentiation depends on the temperature resolution of
#' the melting curve. It is recommended to use a temperature resolution of at
#' least 0.5 degree Celsius. In addition equal distances (e.g., 50 -> 50.5 ->
#' 60 degree Celsius) rather than unequal distances (e.g., 50 -> 50.4 -> 60.1
#' degree Celsius) for the temperature steps are recommended. The parameter
#' \code{n} can be used to increase the temperature resolution of the melting
#' curve data. \code{mcaSmoother} uses the spline function for this purpose. A
#' temperature range for a simple linear background correction. The linear
#' trend is estimated by a robust linear regression using \code{lmrob()}. In
#' case criteria for a robust linear regression are violated \code{lm()} is
#' automatically used. The parameter \code{n} can be combined with the
#' parameter \code{Trange} to make transform all melting curves of question to
#' have the #same range and similar resolution. Optionally a Min-Max
#' normalization between 0 and 1 can be used by setting the option
#' \code{minmax} to \code{TRUE}. This is useful in many situations. For
#' example, if the fluorescence values between samples vary considerably (e.g.,
#' due to high background, different reporter dyes,...), particularly in
#' solution or for better comparison of results.
#' 
#' 
#' @param x is the column of a data frame for the temperature.
#' @param y is the column of a data frame for the fluorescence values.
#' @param bgadj is used to adjust the background signal. This causes
#' \code{mcaSmoother} to use the data subset defined by \code{bg} for the
#' linear regression and background correction.
#' @param bg is used to define the range for the background reduction (e.g.,
#' \code{bg = c(50,55)}, between 50 and 55 degree Celsius)).
#' @param Trange is used to define the temperature range (e.g., \code{Trange =
#' c(50,95)}, between 50 and 95 degree Celsius) for melting curve analysis.
#' @param minmax is used to scale the fluorescence a Min-Max normalization
#' between 0 and 1 can be used by setting the option \code{minmax} to
#' \code{TRUE}.
#' @param df.fact is a factor to smooth the curve. Different strengths can be
#' set using the option \code{df.fact} (f default ~ 0.95). Internally it takes
#' the degree of freedom value from the spline and multiplies it with a factor
#' between 0.6 and 1.1. Values lower than 1 result in stronger smoothed curves.
#' @param n is number of interpolations to take place. This parameter uses the
#' spline function and increases the temperature resolution of the melting
#' curve data.
#' @return \item{xy }{returns a \code{data.frame} with the temperature ("x") in
#' the first and the pre-processed fluorescence values ("y.sp") in the second
#' column.}
#' @author Stefan Roediger
#' @seealso \code{\link{MFIerror}}, \code{\link{lmrob}},
#' \code{\link{smooth.spline}}, \code{\link{spline}}, \code{\link{lm}},
#' \code{\link{approx}}
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
#' @keywords smooth background normalization
#' @examples
#' 
#' # First Example
#' # Use mcaSmoother with different n to increase the temperature 
#' # resolution of the melting curve artificially. Compare the 
#' # influence of the n on the Tm and fluoTm values
#' data(MultiMelt)
#' 
#' Tm	<- vector()
#' fluo	<- vector()
#' for (i in seq(1,3.5,0.5)) {
#'   res.smooth <- mcaSmoother(MultiMelt[, 1], MultiMelt[, 14], n = i)
#'   res <- diffQ(res.smooth)
#'   Tm <- c(Tm, res$Tm)
#'   fluo <- c(fluo, res$fluoTm)
#' }
#' plot(fluo, Tm, ylim = c(76,76.2))
#' abline(h = mean(Tm))
#' text(fluo, seq(76.1,76.05,-0.02), 
#'       paste("n:", seq(3.5,1,-0.5), sep = " "), col = 2)
#' abline(h = c(mean(Tm) + sd(Tm), mean(Tm) - sd(Tm)), col = 2)
#' 
#' legend(-0.22, 76.2, c("mean Tm", "mean Tm +/- SD Tm"), 
#' 	col = c(1,2), lwd = 2)
#' 
#' # Second Example
#' # Use mcaSmoother with different strengths of smoothing 
#' # (f, 0.6 = strongest, 1 = weakest). 
#' data(DMP)
#' plot(DMP[, 1], DMP[,6], 
#'       xlim = c(20,95), xlab = "Temperature",
#'       ylab = "refMFI", pch = 19, col = 8)
#' f <- c(0.6, 0.8, 1.0)
#' for (i in c(1:3)) { 
#'  	lines(mcaSmoother(DMP[, 1],
#'          DMP[,6], df.fact = f[i]),
#'          col = i, lwd = 2)
#' }
#' legend(20, 1.5, paste("f", f, sep = ": "),
#'       cex = 1.2, col =  1:3, bty = "n",
#'       lty = 1, lwd = 4)
#' 
#' # Third Example
#' # Plot the smoothed and trimmed melting curve
#' data(MultiMelt)
#' tmp <- mcaSmoother(MultiMelt[, 1], MultiMelt[, 14])
#' tmp.trimmed <- mcaSmoother(MultiMelt[, 1], MultiMelt[, 14], 
#' 			    Trange = c(49,85))
#' plot(tmp, pch = 19, xlab = "Temperature", ylab = "refMFI", 
#' 	  main = "MLC-2v, mcaSmoother using Trange")
#'   points(tmp.trimmed, col = 2, type = "b", pch = 19)
#'   legend(50, 1, c("smoothed values",
#' 		  "trimmed smoothed values"), 
#' 		   pch = c(19,19), col = c(1,2))
#' 
#' # Fourth Example
#' # Use mcaSmoother with different n to increase the temperature 
#' # resolution of the melting curve. Caution, this operation may 
#' # affect your data negatively if the resolution is set to high. 
#' # Higher resolutions will just give the impression of better 
#' # data quality. res.st uses the default resolution (no 
#' # alteration)
#' # res.high uses the double resolution.
#' data(MultiMelt)
#' res.st <- mcaSmoother(MultiMelt[, 1], MultiMelt[, 14])
#' res.high <- mcaSmoother(MultiMelt[, 1], MultiMelt[, 14], n = 2)
#' 
#' par(fig = c(0,1,0.5,1))
#' plot(res.st, xlab = "Temperature", ylab = "F", 
#' 	main = "Effect of n parameter on the temperature 
#' 	resolution")
#' points(res.high, col = 2, pch = 2)
#'   legend(50, 1, c(paste("default resolution.", nrow(res.st), 
#' 		  "Temperature steps", sep = " "), 
#' 		  paste("double resolution.", nrow(res.high),
#' 		  "Temperature steps", sep = " ")), 
#' 		  pch = c(1,2), col = c(1,2))
#' par(fig = c(0,0.5,0,0.5), new  =  TRUE)
#' diffQ(res.st, plot = TRUE)
#'   text(65, 0.025, paste("default resolution.", nrow(res.st), 
#' 			"Temperature steps", sep = " "))
#' par(fig = c(0.5,1,0,0.5), new = TRUE)
#' diffQ(res.high, plot = TRUE)
#'   text(65, 0.025, paste("double resolution.", nrow(res.high), 
#' 			"Temperature steps", sep = " "))
#' 
#' # Fifth example
#' # Different experiments may have different temperature 
#' # resolutions and temperature ranges. The example uses a 
#' # simulated melting curve with a temperature resolution of 
#' # 0.5 and 1 degree Celsius and a temperature range of 
#' # 35 to 95 degree Celsius.
#' #
#' # Coefficients of a 3 parameter sigmoid model.  Note: 
#' # The off-set, temperature range and temperature resolution 
#' # differ between both simulations. However, the melting 
#' # temperatures should be very 
#' # similar finally.
#' b <- -0.5; e <- 77
#' 
#' # Simulate first melting curve with a temperature 
#' # between 35 - 95 degree Celsius and 1 degree Celsius 
#' # per step temperature resolution.
#' 
#' t1 <- seq(35, 95, 1)
#' f1 <- 0.3 + 4 / (1 + exp(b * (t1 - e)))
#' 
#' # Simulate second melting curve with a temperature 
#' # between 41.5 - 92.1 degree Celsius and 0.5 degree Celsius 
#' # per step temperature resolution.
#' t2 <- seq(41.5, 92.1, 0.5)
#' f2 <- 0.2 + 2 / (1 + exp(b * (t2 - e)))
#' 
#' # Plot both simulated melting curves
#' plot(t1, f1, pch = 15, ylab = "MFI", 
#'      main = "Simulated Melting Curves", 
#'      xlab = "Temperature", col = 1)
#' points(t2, f2, pch = 19, col = 2)
#' legend(50, 1, 
#'        c("35 - 95 degree Celsius, 1 degree Celsius per step", 
#'        "41.5 - 92.1 degree Celsius, 0.5 degree Celsius per step", 
#'        sep = " "), pch = c(15,19), col = c(1,2))
#' 
#' # Use mcaSmoother with n = 2 to increase the temperature 
#' # resolution of the first simulated melting curve. The minmax 
#' # parameter is used to make the peak heights compareable. The 
#' # temperature range was limited between 45 to 90 degree Celsius for 
#' # both simulations
#' 
#' t1f1 <- mcaSmoother(t1, f1, Trange= c(45, 90), minmax = TRUE, n = 2)
#' t2f2 <- mcaSmoother(t2, f2, Trange= c(45, 90), minmax = TRUE, n = 1)
#' 
#' # Perform a MCA on both altered simulations. As expected, the melting
#' # temperature are almost identical.
#' par(mfrow = c(2,1))
#' # Tm 77.00263, fluoTm -0.1245848
#' diffQ(t1f1, plot = TRUE)
#' text(60, -0.08, 
#'      "Raw data: 35 - 95 degree Celsius,\n 1 degree Celsius per step")
#' 
#' # Tm 77.00069, fluoTm -0.1245394
#' diffQ(t2f2, plot = TRUE)
#' text(60, -0.08, "Raw data: 41.5 - 92.1 degree Celsius,
#'       \n 0.5 degree Celsius per step")
#' par(mfrow = c(1,1))
#' 
#' @export mcaSmoother
mcaSmoother <- function(x, y, bgadj = FALSE, bg = NULL, Trange = NULL, 
				minmax = FALSE, df.fact = 0.95, n = NULL) {
  
  old.warn <- options("warn")[["warn"]]
  options(warn = -1)
  # Test if df.fact is within a meaningful range.
  if (df.fact < 0.6 || df.fact > 1.1) 
      stop("df.fact size must be within 0.6 and 1.1.")
  # Test if x and y exist and have identical lengths.
  if (is.null(x)) 
      stop("Enter temperature")
  if (is.null(y)) 
      stop("Enter fluorescence data")
  if (length(x) != length(y)) 
      stop("Use temperature and fluorescence data with same number of elements")
	
  # Test if bg has only two values
  if (!is.null(bg) && length(bg) != 2)
      stop("Use only two temperatures (e.g., bg = c(45,55)) to set the range for the background correction")
  bg <- sort(bg)
  # Test if background adjustment was set and background vector is not empty.
  if ((is.null(bg)) && (bgadj == TRUE))
      stop("Enter temperature background range (e.g., bg = c(45,55)).")
  # Test if Trange has only two values
  if (!is.null(Trange) && length(Trange) != 2)
      stop("Use only two temperatures (e.g., Trange = c(40,70)) to set the range for the melting curve analysis")
  # Test if bg range is in the range of Trange
  Trange <- sort(Trange)
  if (!is.null(Trange) && !is.null(bg) && (bgadj == TRUE)) {
      if ((bg[1] < Trange[1]) || (bg[1] > Trange[2])) {stop("Trange and bg overlapp wrongly")}
      if ((bg[2] < Trange[1]) || (bg[2] > Trange[2])) {stop("Trange and bg overlapp wrongly")}
  }
  # Test if temperature background range is unchanged. If not change to new values.
  if (!is.null(bg) && (bgadj == TRUE)) {
      tmp.data <- data.frame(x,y)
      bg <- c(head(which(tmp.data[, 1] >= bg[1]))[1]:tail(which(tmp.data[, 1] <= bg[2]))[1])
  }
  # Test if temperature range for melting curve analysis is unchanged.
  if (!is.null(Trange)) {
      tmp.data <- data.frame(x,y)
      range <- c(head(which(tmp.data[, 1] >= Trange[1]))[1]:tail(which(tmp.data[, 1] <= Trange[2]))[5])
      x <- tmp.data[range, 1]
      y <- tmp.data[range, 2]
  }
  # Test if y contains missing values. In case of missing values a regression is 
  # used to estimate the missing value.
  if (length(which(is.na(y) == TRUE)) > 0) { 
      y[which(is.na(y))] <- approx(x, y, n = length(x))[["y"]][c(which(is.na(y)))]
  }

  # Smooth the curve with a cubic spline. Takes first the degree of freedom from the cubic spline.
  # The degree of freedom is than used to smooth the curve by a user defined factor.
  df.tmp <- data.frame(smooth.spline(x,y)[["df"]])
  y.sp <- smooth.spline(x, y, df = (df.tmp * df.fact))[["y"]]
  
  if (!is.null(n)) {
      if (n < 0.1 || n > 10) 
	  stop("n must be a number between 0.1 and 10")
	  tmp.xy <- spline(x, y.sp, n = n * length(x))
	  x <- tmp.xy$x
	  y.sp <- tmp.xy[["y"]]
  }

  # If the argument bgadj is set TRUE, bg must be  used to define a temperature range for a linear 
  # background correction. The linear trend is estimated by a robust linear regression using lmrob().
  # In case criteria for a robust linear regression are violated lm() is automatically used.
  if (bgadj) {
      if (class(try(lmrob(y.sp[bg] ~ x[bg]), silent = T)) == "try-error") { 
	  coefficients <- data.frame(lm(y.sp[bg] ~ x[bg])[1]) 
	  } else {
	      lmrob.control <- suppressWarnings(lmrob(y.sp[bg] ~ x[bg])) 
	      if ((class(lmrob.control) != "try-error") && (lmrob.control$converged == TRUE)) { 
			      coefficients <- data.frame(lmrob(y.sp[bg] ~ x[bg])[1]) 
	      } else { 
		  coefficients <- data.frame(lm(y.sp[bg] ~ x[bg])[1]) 
		} 
	    } 
	    y.norm <- y.sp - (coefficients[2, 1] * x + coefficients[1, 1]) # Subtracts the linear trend from the smoothed values.
  } else {
      y.norm <- data.frame(y.sp)
    }
  # Performs a "Min-Max Normalization" between 0 and 1.	  
  if (minmax) {
      y.norm <- (y.norm - min(y.norm)) / (max(y.norm) - min(y.norm))
  }
  
  #restore old warning value
  options(warn = old.warn)
  
  # Returns an object of the type data.frame containing the temperature in the first column 
  # and the pre-processed fluorescence data in the second column.
  data.frame(x,y.norm)
}
