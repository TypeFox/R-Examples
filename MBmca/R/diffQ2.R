#' Calculation of the melting temperatures (Tm, Tm1D2 and Tm2D2) from the first
#' and the second derivative
#' 
#' \code{diffQ2()} calls instances of \code{diffQ()} to calculate the Tm1D2 and
#' Tm2D2. The options are similar to \code{diffQ()}. Both \code{diffQ()} and
#' \code{diffQ2()} return objects of the class \code{list}.  To accessing
#' components of lists is done as described elsewhere either be name or by
#' number. \code{diffQ2} has no standalone plot function. For sophisticated
#' analysis and plots its recommended to use \code{diffQ2} as presented in the
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
#' @param peak shows the peak in the plot.
#' @param plot shows a plot of a single melting curve with (Tm) as vertical
#' line and the second derivatives (Tm1D2 and Tm2D2). To draw multiple curves
#' in a single plot set \code{plot = FALSE} and create and empty plot instead
#' (see examples).
#' @param verbose shows additional information (e.g., first and second
#' approximate derivatives, ranges used for calculation, approximate Tm, Tm1D2,
#' Tm2D2) of the calculation.
#' @param negderiv calculates the negative derivative (default). If
#' \code{FALSE} the positive first negative is calculated.
#' @param deriv shows the first derivative with the color assigned to
#' \code{col} (see examples).
#' @param derivlimits shows the number (n) used to calculate the Tm as points
#' in the plot (see examples).
#' @param derivlimitsline shows the number (n) used to calculate the Tm as line
#' in the plot (see examples).
#' @param vertiline draws a vertical line at the Tms (see examples).
#' @param rsm performs a doubling of the temperature resolution by calculation
#' of the mean temperature and mean fluorescence between successive temperature
#' steps. Note: mcaSmoother has the "n" parameter with a similar but advanced
#' functionality.
#' @param inder Interpolates derivatives using the five-point stencil.  See
#' \code{chipPCR} package for details.
#' @param warn diffQ tries to keep the user as informed as possible about the
#' quality of the analysis. However, in some scenarios are the warning and
#' message about analysis not needed or disturbing.  \code{warn} can be used to
#' stop the flodding of the output.
#' @return
#' 
#' \item{$TmD1}{\code{TmD1} returns a comprehensive list (if parameter verbose
#' is TRUE) with results from the first derivative. The list includes a
#' \code{data.frame} of the derivative ("xy"). The temperature range
#' ("limits.xQ") and fluorescence range ("limits.diffQ") to calculate the peak
#' value. "fluo.x" is the approximate fluorescence at the approximate melting
#' temperature. The calculated melting temperature ("Tm") with the
#' corresponding fluorescence intensity ("fluoTm"). The number of points
#' ("fws") and the adjusted R-squared ("adj.r.squared") to fit.}
#' 
#' \item{$TmD1$Tm}{returns the calculated melting temperature ("Tm") from the
#' first derivative.}
#' 
#' \item{$TmD1$fluoTm}{returns the calculated fluorescence at the calculated
#' melting temperature ("Tm").}
#' 
#' \item{$TmD1$Tm.approx}{returns the approximate melting temperature ("Tm")
#' from the first derivative.}
#' 
#' \item{$TmD1$fluo.x}{returns the approximate fluorescence at the calculated
#' melting temperature ("Tm").}
#' 
#' \item{$TmD1$xy}{is a \code{data.frame} containing in the first column the
#' temperature and in the second column the fluorescence values. Preferably the
#' output from \code{mcaSmoother} is used.}
#' 
#' \item{$TmD1$limits.xQ}{returns a data range of temperature values used to
#' calculate the melting temperature.}
#' 
#' \item{$TmD1$limits.diffQ}{returns a data range of fluorescence values used
#' to calculate the melting temperature.}
#' 
#' \item{$TmD1$adj.r.squared}{returns the adjusted R-squared from the quadratic
#' model fitting function (see also \code{fit}) of the first derivative.}
#' 
#' \item{$TmD1$NRMSE}{returns the normalized root-mean-squared-error (NRMSE)
#' from the quadratic model fitting function (see also \code{fit}) of the first
#' derivative.}
#' 
#' \item{$TmD1$fws}{returns the number of points used for the calculation of
#' the melting temperature of the first derivative.}
#' 
#' \item{$TmD1$devsum}{returns measures to show the difference between the
#' approximate and calculated melting temperature of the first derivative.}
#' 
#' \item{$TmD1$fit}{returns the summary of the results of the quadratic model
#' fitting function of the first derivative.}
#' 
#' \item{$Tm1D2}{returns the "left" melting temperature ("Tm1D2 ") values from
#' the second derivative.}
#' 
#' \item{$Tm1D2$Tm}{returns the "left" calculated melting temperature ("Tm1D2")
#' from the second derivative.}
#' 
#' \item{$Tm1D2$fluoTm}{returns the "left" calculated fluorescence at the
#' calculated melting temperature ("Tm1D2") from the second derivative.}
#' 
#' \item{$Tm1D2$Tm.approx}{returns the "left" approximate melting temperature
#' ("Tm1D2") from the second derivative.}
#' 
#' \item{$Tm1D2$fluo.x}{returns the "left" approximate fluorescence at the
#' calculated melting temperature ("Tm1D2") from the second derivative.}
#' 
#' \item{$Tm1D2$xy}{is a \code{data.frame} containing in the first column the
#' temperature and in the second column the fluorescence values of the "left"
#' melting temperature ("Tm1D2") from the second derivative. Preferably the
#' output from \code{mcaSmoother} is used.}
#' 
#' \item{$Tm1D2$limits.xQ}{returns a data range of temperature values used to
#' calculate the melting temperature of the "left" melting temperature
#' ("Tm1D2") from the second derivative.}
#' 
#' \item{$Tm1D2$limits.diffQ}{returns a data range of fluorescence values used
#' to calculate the melting temperature of the "left" melting temperature
#' ("Tm1D2") from the second derivative.}
#' 
#' \item{$Tm1D2$adj.r.squared}{returns the adjusted R-squared from the
#' quadratic model fitting function (see also \code{fit}) of the "left" melting
#' temperature ("Tm1D2") from the second derivative.}
#' 
#' \item{$Tm1D2$NRMSE}{returns normalized root-mean-squared-error (NRMSE) from
#' the quadratic model fitting function (see also \code{fit}) of the "left"
#' melting temperature ("Tm1D2") from the second derivative.}
#' 
#' \item{$Tm1D2$fws}{ returns the number of points used for the calculation of
#' the melting temperature of the "left" melting temperature ("Tm1D2") from the
#' second derivative.}
#' 
#' \item{$Tm1D2$devsum}{ returns measures to show the difference between the
#' approximate and alculated melting temperature of the "left" melting
#' temperature ("Tm1D2") from the second derivative.}
#' 
#' \item{$Tm1D2$fit}{returns the summary of the results of the quadratic model
#' fitting function of the "left" melting temperature ("Tm1D2") from the second
#' derivative.}
#' 
#' \item{$Tm2D2}{returns the "right" melting temperature ("Tm2D2 ") values from
#' the second derivative.}
#' 
#' \item{$Tm2D2$Tm}{returns the "right" calculated melting temperature
#' ("Tm2D2") from the second derivative.}
#' 
#' \item{$Tm2D2$fluoTm}{returns the "right" calculated fluorescence at the
#' calculated melting temperature ("Tm2D2") from the second derivative.}
#' 
#' \item{$Tm2D2$Tm.approx}{returns the "right" approximate melting temperature
#' ("Tm1D2") from the second derivative.}
#' 
#' \item{$Tm2D2$fluo.x}{returns the "left" approximate fluorescence at the
#' calculated melting temperature ("Tm2D2") from the second derivative.}
#' 
#' \item{$Tm2D2$xy}{is a \code{data.frame} containing in the first column the
#' temperature and in the second column the fluorescence values of the "right"
#' melting temperature ("Tm2D2") from the second derivative.  Preferably the
#' output from \code{mcaSmoother} is used.}
#' 
#' \item{$Tm2D2$limits.xQ}{returns a data range of temperature values used to
#' calculate the melting temperature of the "right" melting temperature
#' ("Tm2D2") from the second derivative.}
#' 
#' \item{$Tm2D2$limits.diffQ}{returns a data range of fluorescence values used
#' to calculate the melting temperature of the "right" melting temperature
#' ("Tm"D2") from the second derivative.}
#' 
#' \item{$Tm2D2$adj.r.squared}{returns the adjusted R-squared from the
#' quadratic model fitting function (see also \code{fit}) of the "right"
#' melting temperature ("Tm2D2") from the second derivative.}
#' 
#' \item{$Tm2D2$NRMSE}{returns normalized root-mean-squared-error (NRMSE) from
#' the quadratic model fitting function (see also \code{fit}) of the "right"
#' melting temperature ("Tm2D2") from the second derivative.}
#' 
#' \item{$Tm2D2$fws}{returns the number of points used for the calculation of
#' the melting temperature of the "right" melting temperature ("Tm2D2") from
#' the second derivative.}
#' 
#' \item{$Tm2D2$devsum}{returns measures to show the difference between the
#' approximate and calculated melting temperature of the "right" melting
#' temperature ("Tm2D2") from the second derivative.}
#' 
#' \item{$Tm2D2$fit}{returns the summary of the results of the quadratic model
#' fitting function of the "right" melting temperature ("Tm2D2") from the
#' second derivative.}
#' 
#' \item{$xTm1.2.D2}{returns only the "left" and right calculated melting
#' temperature ("Tm1D2, Tm2D2") from the second derivative.}
#' 
#' \item{$yTm1.2.D2}{returns only the "left" and right calculated fluorescence
#' ("Tm1D2, Tm2D2") from the second derivative.}
#' 
#' \item{$temperature}{returns measures to investigate the temperature
#' resolution of the melting curve. Raw fluorescence measurements at irregular
#' temperature resolutions (intervals) can introduce artifacts and thus lead to
#' wrong melting point estimations.}
#' 
#' \item{$temperature$T.delta}{returns the difference between two successive
#' temperature steps.}
#' 
#' \item{$temperature$mean.T.delta}{returns the mean difference between two
#' temperature steps.}
#' 
#' \item{$temperature$sd.T.delta}{returns the standard deviation of the
#' temperature.}
#' 
#' \item{$temperature$RSD.T.delta}{returns the relative standard deviation
#' (RSD) of the temperature in percent.}
#' @author Stefan Roediger
#' @seealso \code{\link{diffQ}}, \code{\link{mcaSmoother}}
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
#' @keywords Tm
#' @examples
#' 
#' # First Example
#' # Plot the first and the second derivative melting curves of MLC-2v
#' # for a single melting curve. Should give a warning message but the graph 
#' # will show you that the calculation is ok
#' data(MultiMelt)
#' tmp <- mcaSmoother(MultiMelt[, 1], MultiMelt[, 14])
#' diffQ2(tmp, fct = min, verbose = FALSE, plot = TRUE)
#' 
#' # Second Example
#' # Calculate the maximum fluorescence of a melting curve, Tm, 
#' # Tm1D2 and Tm2D2 of HPRT1 for 12 microbead populations and assign the 
#' # values to the matrix HPRT1
#' data(MultiMelt)
#' HPRT1 <- matrix(NA,12,4,
#' 	dimnames = list(colnames(MultiMelt[, 2L:13]),
#' 	    c("Fluo", "Tm", "Tm1D2", "Tm2D2")))
#' for (i in 2L:13) {
#'        tmp <- mcaSmoother(MultiMelt[, 1], 
#'                      MultiMelt[, i])
#'        tmpTM <- diffQ2(tmp, fct = min, verbose = TRUE)
#'        HPRT1[i-1, 1] <- max(tmp$y)
#'        HPRT1[i-1, 2] <- tmpTM$TmD1$Tm
#'        HPRT1[i-1, 3] <- tmpTM$Tm1D2$Tm
#'        HPRT1[i-1, 4] <- tmpTM$Tm2D2$Tm
#' }
#' HPRT1
#' 
#' # Third Example
#' # Use diffQ2 to determine the second derivative.
#' 
#' data(MultiMelt)
#' HPRT1 <- matrix(NA,12,4,
#' 	dimnames = list(colnames(MultiMelt[, 2L:13]),
#' 	    c("Fluo", "Tm", "Tm1D2", "Tm2D2")))
#' for (i in 2L:13) {
#'        tmp <- mcaSmoother(MultiMelt[, 1], 
#'                      MultiMelt[, i])
#'        tmpTM <- diffQ2(tmp, fct = min, verbose = TRUE)
#'        HPRT1[i-1, 1] <- max(tmp[["y.sp"]])
#'        HPRT1[i-1, 2] <- tmpTM[["TmD1"]][["Tm"]]
#'        HPRT1[i-1, 3] <- tmpTM[["Tm1D2"]][["Tm"]]
#'        HPRT1[i-1, 4] <- tmpTM[["Tm2D2"]][["Tm"]]
#' }
#' plot(HPRT1[, 1], HPRT1[, 2], 
#'        xlab = "refMFI", ylab = "Temperature", 
#'        main = "HPRT1", xlim = c(2.1,2.55),
#'        ylim = c(72,82), pch = 19,
#'        col = 1:12, cex = 1.8)
#' points(HPRT1[, 1], HPRT1[, 3], pch = 15)
#' points(HPRT1[, 1], HPRT1[, 4], pch = 15)
#' abline(lm(HPRT1[, 2] ~ HPRT1[, 1]))
#' abline(lm(HPRT1[, 3] ~ HPRT1[, 1]))
#' abline(lm(HPRT1[, 4] ~ HPRT1[, 1]))
#' 
#' # Fourth Example
#' # Use diffQ2 with inder parameter to determine the second derivative.
#' data(MultiMelt)
#' 
#' tmp <- mcaSmoother(MultiMelt[, 1], MultiMelt[, 14])
#' diffQ2(tmp, fct = min, verbose = FALSE, plot = TRUE, inder = FALSE)
#' diffQ2(tmp, fct = min, verbose = FALSE, plot = TRUE, inder = TRUE)
#' par(mfrow = c(1,1))
#' 
#' @export diffQ2
diffQ2 <- function(xy, fct = max, fws = 8, col = 2, plot = FALSE, 
                   verbose = FALSE, peak = FALSE, 
                   deriv = FALSE, negderiv = TRUE, derivlimits = FALSE, 
                   derivlimitsline = FALSE, vertiline = FALSE, 
                   rsm = FALSE, inder = FALSE, warn = TRUE) {
  # Test if fws (number of neighbors) is within a meaningful range.
  options(warn = -1)
  fws <- round(fws)
  if (fws < 2 || fws > 8) 
    stop("Fit window size must be within 2 and 8.")
  
  # Calls instances of diffQ to calculate the Tm of the first and the 
  # second derivatives.The arguments are similar to diffQ.
  # Calculates the first derivative and provides the starting information 
  # for the second derivative.
  TmD1 <- diffQ(xy, fct = fct, fws = fws, negderiv = negderiv, 
                verbose = TRUE, rsm = rsm, inder = inder, 
                warn = warn)
  
  if (TmD1[["temperature"]][["mean.T.delta"]] >= 0.5) {
    tmp.xy <- predict(smooth.spline(TmD1[["xy"]][, 1], TmD1[["xy"]][, 2]), 
                      seq(min(TmD1[["xy"]][, 1]), max(TmD1[["xy"]][, 1]), 
                          TmD1[["temperature"]][["mean.T.delta"]] / 1.3))
    xy.smoothed <- data.frame(tmp.xy[["x"]], tmp.xy[["y"]])
  } else (xy.smoothed <-  data.frame(TmD1[["xy"]][, 1], 
                                     smooth.spline(TmD1[["xy"]][, 1], 
                                                   TmD1[["xy"]][, 2])[["y"]])
  )
  
  # Calculates the second derivative and from the first derivative for the 
  # minimum of the melting curve.
  Tm1D2 <- diffQ(xy.smoothed, fct = min, fws = fws, verbose =  TRUE, 
                 col = col, peak = peak, deriv = deriv, negderiv = FALSE, 
                 derivlimits = derivlimits, derivlimitsline = derivlimitsline, 
                 vertiline = vertiline, inder = inder, 
                 warn = warn)
  # Calculates the second derivative and from the first derivative for the 
  # maximum of the melting curve.
  Tm2D2 <- diffQ(xy.smoothed, fct = max, fws = fws, verbose =  TRUE, 
                 col = col, peak = peak, deriv = deriv, negderiv = FALSE, 
                 derivlimits = derivlimits, derivlimitsline = derivlimitsline, 
                 vertiline = vertiline, inder = inder, 
                 warn = warn)
  
  # Vectors of the two melting temperatures of the second derivative.
  x <- c(Tm1D2[["Tm"]], Tm2D2[["Tm"]])
  # Vectors of the two intensities at the melting temperatures of the 
  # second derivative.
  # y <- c(coeflm2.y.1, coeflm2.y.2)
  y <- c(Tm1D2[["fluoTm"]], Tm2D2[["fluoTm"]])
  
  # Polynom to fit the area of the calculated Tm
  poly.fct.TmD1 <- function(xi) TmD1[["fit"]][[4]][1] + TmD1[["fit"]][[4]][2] * xi + TmD1[["fit"]][[4]][3] * xi^2
  poly.fct.Tm1D2 <- function(xi) Tm1D2[["fit"]][[4]][1] + Tm1D2[["fit"]][[4]][2] * xi + Tm1D2[["fit"]][[4]][3] * xi^2
  poly.fct.Tm2D2 <- function(xi) Tm2D2[["fit"]][[4]][1] + Tm2D2[["fit"]][[4]][2] * xi + Tm2D2[["fit"]][[4]][3] * xi^2
  
  
  if (plot) {
    # Plot the first derivative
    par(fig = c(0,1,0.475,1))
    plot(TmD1[["xy"]], xlab = "Temperature", ylab = "-d(F) / dT", type = "b")
    abline(v = (TmD1[["Tm"]]), col = "grey", lwd = 1.25)
    abline(v = (Tm1D2[["Tm"]]), col = "grey")
    abline(v = (Tm2D2[["Tm"]]), col = "grey")
    points(TmD1[["limits.xQ"]], TmD1[["limits.diffQ"]], col = "orange", pch = 19)
    points(TmD1[["Tm"]], TmD1[["fluoTm"]], pch = 19, col = 2)
    curve(poly.fct.TmD1, TmD1[["limits.xQ"]][1], TmD1[["limits.xQ"]][length(TmD1[["limits.xQ"]])], col = "red", add = TRUE)
    if (derivlimits) {
      points(TmD1[["limits.xQ"]], TmD1[["limits.diffQ"]], cex = 1, pch = 19, col = col)
    }
    
    # Plot the second derivative
    par(fig = c(0,1,0,0.525), new = TRUE)
    plot(Tm1D2[["xy"]], xlim = c(min(TmD1[["xy"]][, 1]), max(TmD1[["xy"]][, 1])), 
         xlab = "Temperature", ylab = "-d^2(F) / dT^2", type = "b")
    abline(v = (TmD1[["Tm"]]), col = "grey")
    points(Tm1D2[["limits.xQ"]], Tm1D2[["limits.diffQ"]], col = "green", pch = 19)
    points(Tm2D2[["limits.xQ"]], Tm2D2[["limits.diffQ"]], col = "blue", pch = 1)
    points(Tm1D2[["Tm"]], Tm1D2[["fluo.x"]], pch = 19, col = 2)
    points(Tm2D2[["Tm"]], Tm2D2[["fluo.x"]], pch = 19, col = 2)
    curve(poly.fct.Tm1D2, Tm1D2[["limits.xQ"]][1], Tm1D2[["limits.xQ"]][length(Tm1D2[["limits.xQ"]])], col = "red", add = TRUE)
    curve(poly.fct.Tm2D2, Tm2D2[["limits.xQ"]][1], Tm2D2[["limits.xQ"]][length(Tm2D2[["limits.xQ"]])], col = "red", add = TRUE)
  }
  
  # Returns an object of the type list containing the data and data.frames from above including the approximate 
  # difference quotient values, melting temperatures of the first derivative and the second derivative, intensities and used neighbors.
  if (verbose) {
    res <- list(TmD1 = TmD1, Tm1D2 = Tm1D2, Tm2D2 = Tm2D2, xTm1.2.D2 = x, 
         yTm1.2.D2 = y, temperature = TmD1[["temperature"]])
  } else {
    res <- list(Tm = TmD1[["Tm"]], fluoTm = TmD1[["fluoTm"]], 
         xTm1.2.D2 = x, yTm1.2.D2 = y)
  }
  class(res) <- "diffQ2object"
  res
  
}
