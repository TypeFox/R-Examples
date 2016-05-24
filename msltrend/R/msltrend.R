    # ----------------------------------------------------------------
    # Mean Sea Level Analysis using annual average time series data files
    # Developed by P.J.Watson (27 October 2015)
    # ----------------------------------------------------------------
    #' Projected sea level rise integrated with historical record.
    #'
    #' @param object of class \dQuote{msl.trend} (see \code{\link{msl.trend}} and
    #' \code{\link{s}}).
    #' @param slr numeric, enables a user defined amount of projected sea level
    #' rise in millimetres. The user range is [200 to 1500] where 800 is the
    #' default setting.
    #' @param plot logical, if \sQuote{TRUE} then the original time series is
    #' plotted to the screen along with the trend component, the result of gap
    #' filling (where necessary) and the added quantum of sea level rise
    #' selected. 95\% confidence intervals have also been applied. Default = TRUE.
    #'
    #' @details This routine adds a user specified quantum of sea level rise
    #' from the end of the deconstructed historical record to the year 2100.
    #' All internal parameters captured in the \code{\link{msl.trend}} object
    #' are passed directly to \code{\link{msl.forecast}}.
    #'
    #' @return An object of class \dQuote{msl.forecast} is returned with the
    #' following elements:
    #' \describe{
    #'  \item{\strong{$Station.Name: }}{the name of the data record.}
    #'  \item{\strong{$Summary: }}{a summary data frame of relevant attributes
    #'  relating to the trend and the inputted annual average data set extended
    #'  to 2100 with projected sea level rise, including:}
    #'  \itemize{
    #'   \item{$Year: input data; }
    #'   \item{$MSL: input data; }
    #'   \item{$Trend: mean sea level trend; }
    #'   \item{$TrendSD: standard deviation of the determined mean sea level
    #'   trend; }
    #'   \item{$Vel: velocity (or first derivative) of mean sea level trend
    #'   (mm/year); }
    #'   \item{$VelSD: standard deviation of the velocity of the mean sea level
    #'   trend; }
    #'   \item{$Acc: acceleration (or second derivative) of mean sea level trend
    #'   (mm/year/year); }
    #'   \item{$AccSD: standard deviation of the acceleration of the mean sea
    #'   level trend; and }
    #'   \item{$FilledTS: gap-filled time series (where necessary). }
    #'    }
    #' }
    #'
    #' \describe{
    #'  \item{\strong{$Velocity: }}{outputs the peak velocity and the year in
    #'  which it occurs.}
    #'  \item{\strong{$Acceleration: }}{outputs the peak acceleration and the
    #'  year in which it occurs.}
    #'  \item{\strong{$Historical.Record: }}{outputs details of the start, end
    #'  and length of the input data set.}
    #'  \item{\strong{$Historical.Fillgaps: }}{outputs the extent of missing data
    #'  (years) in the original record and the gap filling method used (where
    #'  necessary).}
    #'  \item{\strong{$Projected.SLR: }}{details the amount of sea level rise
    #'  applied between the end of the historical record and the year 2100.}
    #'  \item{\strong{$Bootstrapping.Iterations: }}{outputs the number of
    #'  iterations used to generate the respective standard deviations for error
    #'  margins.}
    #'  }
    #'
    #' @seealso \code{\link{msl.trend}}, \code{\link{msl.plot}}, \code{\link{msl.pdf}},
    #' \code{\link{summary}}, \code{\link{Balt}}, \code{\link{s}}, \code{\link{t}}.
    #'
    #' @examples
    #' # -------------------------------------------------------------------------
    #' # Isolate trend from Baltimore record, filling gaps with spline interpolation,
    #' # 500 iterations and adding 1000 mm of slr to 2100. Use raw 'Balt.csv' data file.
    #' # Note: ordinarily user would call 'File.csv' direct from working directory
    #' # using the following sample code:
    #' # s <- msl.trend('Balt.csv', fillgaps = 3, iter = 500, 'BALTIMORE, USA')
    #' # t <- msl.forecast(s, slr = 1000)
    #' # -------------------------------------------------------------------------
    #'
    #' data(s) # msl.trend object from above-mentioned example
    #' data(t) # msl.forecast object from above-mentioned example
    #' str(t) # check structure of msl.forecast object
    #' msl.plot(s, type=2) # check screen output of gapfilling and trend estimate
    #' msl.plot(t, type=2) # check screen output of adding 1000 mm of sea level rise
    #'
    #' @export
msl.forecast <- function(object, slr = 800, plot = TRUE) {
    # -----------------------------------------------------------------
    grDevices::graphics.off()  # close active graphics windows

    summ <- object$Summary  # 'msl.trend' object
    n <- length(summ[, 1])
    # -----------------------------------------------------------------
    # requires 'msl.trend' class object
    # sea level rise range 200 - 1500 mm to 2100
    # -----------------------------------------------------------------
    if (object$Station.Name == "Station Name not entered") {
        object$Station.Name <- NULL
    } else {
        object$Station.Name <- object$Station.Name
    }
    # -----------------------------------------------------------------
    # if not a msl.trend object
    if (class(object) != "msl.trend") {
        stop("input object is not an msl.trend object: forecasting terminated")
    }
    if (slr >= 200 && slr <= 1500) {
        # within designated range of slr
        slr <- slr
    } else {
        # If slr not entered or entered outside range
        print("default slr setting (800mm to 2100) applied")
        slr <- 800
    }
    # -----------------------------------------------------------------
    # quadratic from last point in historical record to 2100
    # s = ut + 0.5t^2 , s = slr, t = time,
    # u = last velocity in historical record
    len2 <- 2100 - summ[n, 1]  # last year of historical record
    u <- summ[n, 5]
    a <- 2 * (slr - u * len2)/len2 ^ 2
    t2 <- c(1:(len2 + 30))
    predslr <- u * t2 + 0.5 * a * t2 ^ 2
    predslr <- predslr + summ[n, 3]  # normalise to trend of historical ts
    msl <- c(summ[, 3], predslr)  # complete trend to 2130
    len3 <- c(summ[1, 1]:2130)
    n2 <- length(len3)
    ts1 <- stats::ts(msl, start = summ[1, 1])  # consolidated ts of trend + proj slr
    trend <- ts1
    # -----------------------------------------------------------------
    # fitted smooting spline to trend
    # to determine velocity and acceleration
    k <- round((length(ts1)/8), 0)  # DoF for spline fitting (1 every 8 years)
    #
    ss <- stats::smooth.spline(trend, df = k)
    ssVEL <- stats::predict(ss, type = "response", deriv = 1)
    ssVEL <- ssVEL$y  # Velocity vector
    ssACC <- stats::predict(ss, type = "response", deriv = 2)
    ssACC <- ssACC$y  # Acceleration vector
    # -----------------------------------------------------------------
    # standard errors bootsrapping routines
    if (requireNamespace("Rssa", quietly = TRUE)) {
      Rssa::grouping.auto
      Rssa::igapfill
      Rssa::reconstruct
      Rssa::ssa
    }
    reps <- object$Bootstrapping.Iterations
    TS <- vector(mode = "list", length = reps)
    VEL <- vector(mode = "list", length = reps)
    ACC <- vector(mode = "list", length = reps)
    # -----------------------------------------------------------------
    for (i in 1:reps) {
        samp <- sample(summ$Resids, n2, replace = TRUE)
        newTS <- trend + samp
        model2.ssa <- Rssa::ssa(newTS, kind = "1d-ssa", svd.method = "auto")
        ssa2.auto <- Rssa::grouping.auto(model2.ssa, base = "series",
                                   groups = list(1:30), freq.bins = list(0.01),
                                   threshold = 0.75)
        recon2 <- Rssa::reconstruct(model2.ssa, groups = ssa2.auto)
        trend2 <- recon2$F1
        TS[[i]] <- trend2
        ss2 <- stats::smooth.spline(trend2, df = k)
        ssVEL2 <- stats::predict(ss2, type = "response", deriv = 1)
        VEL[[i]] <- ssVEL2$y  ## Isolate Velocity Vector
        ssACC2 <- stats::predict(ss2, type = "response", deriv = 2)
        ACC[[i]] <- ssACC2$y  ## Isolate acceleration Vector
    }
    # -----------------------------------------------------------------
    # standard deviations
    TSboot <- as.data.frame(do.call(cbind, TS))
    TSsd <- apply(TSboot[, c(1:reps)], 1, sd)
    VELboot <- as.data.frame(do.call(cbind, VEL))
    VELsd <- apply(VELboot[, c(1:reps)], 1, sd)
    ACCboot <- as.data.frame(do.call(cbind, ACC))
    ACCsd <- apply(ACCboot[, c(1:reps)], 1, sd)
    # -----------------------------------------------------------------
    # add vectors for initial time series out to 2130 to
    # match msl.trend data frame output
    vec2 <- c(summ[1, 1]:2130)  # year
    vec3 <- vector(mode = "numeric", length = length(t2))
    vec3[vec3 == 0] <- NA
    vec3 <- c(summ$MSL, vec3)  ## vector to expand MSL record to 2130
    vec4 <- vector(mode = "numeric", length = length(t2))
    vec4[vec4 == 0] <- NA
    vec4 <- c(summ$FilledTS, vec4)  ## vector to expand FilledTS record to 2130
    # -----------------------------------------------------------------
    # summary dataframe for forecast (summ2)
    summ2 <- as.data.frame(cbind(vec2, vec3, as.vector(trend), TSsd, ssVEL,
                                 VELsd, ssACC, ACCsd, vec4))
    colnames(summ2) <- c("Year", "MSL", "Trend", "TrendSD", "Vel", "VelSD",
                         "Acc", "AccSD", "FilledTS")
    # NA's at first and last 3 acceleration values (not accurate with spline)
    summ2$Acc[1:3] <- NA
    summ2$AccSD[1:3] <- NA
    summ2 <- summ2[summ2$Year < 2101, ]
    # -----------------------------------------------------------------
    # summary object
    n1 <- length(summ2[, 1])
    n2 <- n1 - n
    summ3 <- summ2[4:n1, ]  # non NA portion of summ
    object2 <- NULL
    object2$Station.Name <- object$Station.Name
    object2$Summary <- summ2
    object2$Velocity <- list(Peak.mm.yr = max(summ2$Vel),
                             Year = with(summ2, Year[Vel == max(Vel)]))
    object2$Acceleration <- list(Peak.mm.yr.yr = max(summ3$Acc),
                                 Year = with(summ3,Year[Acc == max(summ3$Acc)]))
    object2$Historical.Record <- object$Record.Length
    object2$Historical.Fillgaps <- object$Fillgaps
    object2$Projected.SLR <- list(Slr.mm = slr,
                                  Start.period = object$Record.Length$End,
                                  End = 2100, Years = n2)
    object2$Bootstrapping.Iterations <- object$Bootstrapping.Iterations
    class(object2) <- "msl.forecast"
    # -----------------------------------------------------------------
    # plot routines
    slrport <- summ2$Trend
    slrport[1:n] <- NA
    m1 <- object$Record.Length$End  # last year of historical record
    msl.sub <- summ2$MSL[summ2$Year <= m1]
    # -----------------------------------------------------------------
    # conditioning parameter for plot routines
    if (any(is.na(msl.sub)) == FALSE) {
        # If time series doesn't contain missing data
        p <- 1  # conditioning parameter
    } else {
        p <- 0
    }
    # -----------------------------------------------------------------
    # If plot option not entered or entered outside range
        if (plot == TRUE | plot == FALSE) {
        plot <- plot
    } else {
        print("default TYPE setting applied")
        plot <- TRUE
    }
    if (summ$Trend[n] - summ$Trend[1] < 0) {
        # check slope of graph for locating chart 1 legends
        laba <- paste("bottomleft")  # Chart key to bottomleft
        labb <- paste("topright")  # SLR notation to topright
    } else {
        # default setting
        laba <- paste("bottomright")  # Chart key to bottomright
        labb <- paste("topleft")  # SLR notation to topleft
    }
    labc <- paste0("Projected Sea Level Rise = ", slr, " mm")  # legend for slr
    # -----------------------------------------------------------------
    if (plot == TRUE) {
      if (requireNamespace("plyr", quietly = TRUE)) {
      plyr::round_any
    }
    # -----------------------------------------------------------------
    # setting up x-axis scale
      if (n < 100) {
      xtic = 10  # year ticks on x-axis
      xlo <- plyr::round_any(min(summ2[, 1]), 10, floor)
      xhi <- plyr::round_any(max(summ2[, 1]), 10, ceiling)
    } else {
      # default
      xtic = 20
      xlo <- plyr::round_any(min(summ2[, 1]), 20, floor)
      xhi <- plyr::round_any(max(summ2[, 1]), 20, ceiling)
    }
    # -----------------------------------------------------------------
    # setting up y-axis scale
    if (p == 0) {
      # gap filling routine required
      M1 <- max(max(summ2$MSL, na.rm = TRUE),
                max(summ2$Trend + 1.96 * summ2$TrendSD),
                max(summ2$FilledTS, na.rm = TRUE))
      M2 <- min(min(summ2$MSL, na.rm = TRUE),
                min(summ2$Trend - 1.96 * summ2$TrendSD),
                min(summ2$FilledTS, na.rm = TRUE))
    }
    if (p == 1) {
      # no gaps in record
      M1 <- max(max(summ2$MSL, na.rm = TRUE), max(summ2$Trend + 1.96 * summ2$TrendSD))
      M2 <- min(min(summ2$MSL, na.rm = TRUE), min(summ2$Trend - 1.96 * summ2$TrendSD))
    }
    ylen <- M1 - M2
    if (ylen <= 500) {
      # setting up y-axis plotting parameters
      ytic = 50  # year ticks on y-axis
      ylo <- plyr::round_any(M2, 50, floor)
      yhi <- plyr::round_any(M1, 50, ceiling)
    }
    if (ylen > 500 & ylen <= 1000) {
      # setting up y-axis plotting parameters
      ytic = 100  # year ticks on y-axis
      ylo <- plyr::round_any(M2, 100, floor)
      yhi <- plyr::round_any(M1, 100, ceiling)
    }
    if (ylen > 1000) {
      # setting up y-axis plotting parameters
      ytic = 200  # year ticks on y-axis
      ylo <- plyr::round_any(M2, 200, floor)
      yhi <- plyr::round_any(M1, 200, ceiling)
    }
    ylim = c(ylo, yhi)
    # -----------------------------------------------------------------
    if (p == 0) {
        # gap filling routine required
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2, 0.5), las = 1)
        graphics::plot(summ2$Year, summ2$MSL, type = "l", lty = 0, xlab = " ",
             ylim = ylim, xaxt = "n", yaxt = "n", ylab = " ",
             main = object$Station.Name)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::legend(laba, bg = "white",
               legend = c("KEY", "Annual Average Data", "MSL Trend",
                          "95% Confidence Intervals", "Gap Filling",
                          "Sea Level Rise"), text.font = c(2, 1, 1, 1, 1, 1),
               lty = c(0, 1, 1, 2, 1, 1), lwd = c(1, 1, 3, 1, 1, 3),
               col = c("black", "darkgrey", "black", "black", "red", "blue"),
               cex = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8))
        graphics::legend(labb, legend = labc, bty = "n")  # dispaly quantum of slr
        graphics::lines(summ2$Year, summ2$FilledTS, col = "red")
        graphics::lines(summ2$Year, summ2$MSL, col = "darkgrey")
        graphics::lines(summ2$Year, summ2$Trend, lwd = 2)
        graphics::lines(summ2$Year, slrport, lwd = 2, col = "blue")
        graphics::lines(summ2$Year, summ2$Trend + 1.96 * summ2$TrendSD, lty = 2)
        graphics::lines(summ2$Year, summ2$Trend - 1.96 * summ2$TrendSD, lty = 2)
        graphics::title(ylab = "Relative Mean Sea Level (mm)", font.lab = 2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
    }
    if (p == 1) {
        # no gaps in record
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2, 0.5), las = 1)
        graphics::plot(summ2$Year, summ2$MSL, type = "l", xlab = " ", ylim = ylim,
             xaxt = "n", yaxt = "n", ylab = " ", main = object$Station.Name)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::legend(laba, bg = "white",
               legend = c("KEY", "Annual Average Data", "MSL Trend",
                          "95% Confidence Intervals", "Sea Level Rise"),
               text.font = c(2, 1, 1, 1, 1, 1), lty = c(0, 1, 1, 2, 1),
               lwd = c(1, 1, 3, 1, 3), col = c("black", "darkgrey",
                                               "black", "black", "blue"),
               cex = c(0.8, 0.8, 0.8, 0.8, 0.8))
        graphics::legend(labb, legend = labc, bty = "n")  # dispaly quantum of slr
        graphics::lines(summ2$Year, summ2$MSL, col = "darkgrey")
        graphics::lines(summ2$Year, summ2$Trend, lwd = 2)
        graphics::lines(summ2$Year, slrport, lwd = 2, col = "blue")
        graphics::lines(summ2$Year, summ2$Trend + 1.96 * summ2$TrendSD, lty = 2)
        graphics::lines(summ2$Year, summ2$Trend - 1.96 * summ2$TrendSD, lty = 2)
        graphics::title(ylab = "Relative Mean Sea Level (mm)", font.lab = 2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
     }
    }
    # -----------------------------------------------------------------
    # no plotting
    if (plot == FALSE) {
        print("no plotting required")
    }
    return(object2)
}

#' Isolate trend component from mean sea level records.
#'
#' @param file csv format input file \strong{with no header row} of annual
#' average water levels. This file must contain 2 columns with the first column
#' the time period (in years) and the second column annual average ocean water
#' levels (in millimetres). Missing data must be denoted by \dQuote{NA}. Missing
#' data and maximum missing data gap are limited to 15\% and 5\%, respectively,
#' of the data record. The minimum length data record processed by the package
#' is 80 years.
#'
#' \strong{Warning: } If input data files do not conform to these pre-conditions,
#' the analysis will be terminated. It should be further noted that the existence
#' of quasi 60 year oscillations in global mean sea level have been well recognised
#' in the literature. Therefore, in order to be effective for climate change and
#' sea level research, only input files with a minimum length exceeding 80 years
#' have been considered in order that the package can identify and isloate such
#' signals.
#'
#' @param station_name character string, providing the name of the data record.
#'
#' \strong{Note: }This field can be left blank, however, it is retained for use
#' in banner labelling of all plotting and pdf outputs.
#'
#' @param fillgaps numeric, provides 3 alternative gap filling procedures for
#' missing data. The default procedure (fillgaps = 1) is based on iterative gap
#' filling using Singular Spectrum Analysis (refer \code{\link[Rssa]{igapfill}})
#' . The alternatives (2 and 3) are based on linear interpolation and cubic
#' spline interpolation, respectively (refer \code{\link[zoo]{na.approx}}).
#'
#' \strong{Note: }Gap filled portions of the time series are denoted in red on
#' the default screen plot. This is done specifically to provide ready visual
#' observation to discern if the selected gap filling method provides an
#' appropriate estimate within the gaps in keeping with the remainder of the
#' historical record. Depending on the nature of the record and extent of gaps,
#' some trial and error between alternatives might be necessary to optimise gap
#' filling.
#'
#' @param iter numeric, enables a user defined number of iterations for
#' bootstrapping to determine error margins. The user range is [500 to 10000]
#' where 10000 is the default setting.
#'
#' \strong{Warning: }Although the default setting provides a more accurate basis
#' for estimating error margins, the degree of iterations slows the analysis and
#' can take several minutes to run.
#'
#' @param plot logical, if \sQuote{TRUE} then the original time series is
#' plotted to the screen along with the trend component and the result of gap
#' filling (where necessary). 95\% confidence intervals have also been applied.
#' Default = TRUE.
#'
#' @details This is the key entry point to the package. This function deconstructs
#' annual average time series data into a trend and associated real-time velocities
#' and accelerations, filling necessary internal structures to facilitate all
#' other functions in this package. The trend is isloated using Singular Spectrum
#' Analysis, in particular, aggregating components whose low frequency band [0 to
#' 0.01] exceed a threshold contribution of 75\%. Associated velocities and
#' accelerations are determined through the fitting of a cubic smoothing spline
#' to the trend with 1 degree of freedom per every 8 years of record length.
#' Refer Watson (2016a,b) for more detail.
#'
#' @return An object of class \dQuote{msl.trend} is returned with the following
#' elements:
#' \describe{
#'  \item{\strong{$Station.Name: }}{the name of the data record.}
#'  \item{\strong{$Summary: }}{a summary data frame of the relevant attributes
#'  relating to the trend and
#'   the inputted annual average data set, including:}
#'  \itemize{
#'   \item{$Year: input data; }
#'   \item{$MSL: input data; }
#'   \item{$Trend: mean sea level trend; }
#'   \item{$TrendSD: standard deviation of the determined mean sea level
#'   trend; }
#'   \item{$Vel: velocity (or first derivative) of mean sea level trend
#'   (mm/year); }
#'   \item{$VelSD: standard deviation of the velocity of the mean sea level
#'   trend; }
#'   \item{$Acc: acceleration (or second derivative) of mean sea level trend
#'   (mm/year/year); }
#'   \item{$AccSD: standard deviation of the acceleration of the mean sea level
#'   trend; }
#'   \item{$Resids: time series of uncorrelated residuals; and }
#'   \item{$FilledTS: gap-filled time series (where necessary). }
#'    }
#' }
#'
#' \describe{
#'  \item{\strong{$Velocity: }}{outputs of the peak velocity and the year in
#'  which it occurred.}
#'  \item{\strong{$Acceleration: }}{outputs of the peak acceleration and the
#'  year in which it occurred.}
#'  \item{\strong{$Record.Length: }}{outputs details of the start, end and
#'  length of the input data set.}
#'  \item{\strong{$Fillgaps: }}{outputs the extent of missing data (years) in
#'  the original record and the gap filling method used (where necessary).}
#'  \item{\strong{$Bootstrapping.Iterations: }}{outputs the number of iterations
#'  used to generate the respective standard deviations for error margins.}
#'  \item{\strong{$Changepoints: }}{outputs the number and time at which
#'  changepoints in the variance of the uncorrelated residuals occur (if any).
#'  Where changepoints are identified, block bootstrapping procedures are used
#'  with residuals quarantined between changepoints.}
#'    }
#'
#' @references Watson, P.J., 2016a. Identifying the best performing time series
#' analytics for sea-level research. In: \emph{Time Series Analysis and
#' Forecasting, Contributions to Statistics}, ISBN 978-3-319-28725-6, Springer
#' International Publishing (in press).
#'
#' Watson, P.J., 2016b. How to improve estimates of real-time acceleration in
#' the mean sea level signal. In: Vila-Concejo, A., Bruce, E., Kennedy, D.M.,
#' and McCarroll, R.J. (eds.), Proceedings of the 14th International Coastal
#' Symposium (Sydney, Australia). \emph{Journal of Coastal Research},
#' Special Issue, No. 75. Coconut Creek (Florida), ISSN 0749-0208 (in press).
#'
#' @seealso \code{\link{msl.forecast}}, \code{\link{msl.plot}},
#' \code{\link{msl.pdf}}, \code{\link{summary}}, \code{\link{Balt}}, \code{\link{s}}.
#'
#' @examples
#' # -------------------------------------------------------------------------
#' # Isolate trend from Baltimore record, filling gaps with spline interpolation and
#' # 500 iterations. Use raw 'Balt.csv' data file. Note: ordinarily user would call
#' # 'File.csv' direct from working directory using the following sample code:
#' # s <- msl.trend('Balt.csv', fillgaps = 3, iter = 500, 'BALTIMORE, USA') # DONT RUN
#' # -------------------------------------------------------------------------
#'
#' data(s) # msl.trend object from above-mentioned example
#' str(s) # check structure of msl.trend object
#' msl.plot(s, type=2) # check screen output of gapfilling and trend estimate
#'
#'
#' @export
msl.trend <- function(file, station_name = " ", fillgaps = 1, iter = 10000,
                      plot = TRUE) {
  # -----------------------------------------------------------------
  grDevices::graphics.off()  # close active graphics windows
  y <- utils::read.csv(file, header = FALSE)
  if (class(y[1, 1]) != "integer" | class(y[1, 2]) != "integer") {
    # If non numeric header row included in input file
    stop("non-numeric header row in input file : program terminated. Check
         manual for input file format requirements.")
  }
  ts1 <- stats::ts(y[2], start = y[1, 1])
  n <- length(ts1)
  n1 <- sum(is.na(ts1))
  L1 <- ceiling(0.2 * n)  # round up to full number window for ssa shape
  p <- 0  # conditioning parameter (assumes missing data to start)
  # -----------------------------------------------------------------
  # error handling
  if (n < 80) {
    # If too short
    stop("length of time series less than 80 years: program terminated. Check
         manual for input file format requirements.")
  }
  if (min(diff(y[[1]]), na.rm = TRUE) < 1) {
    # If time steps less than annual
    stop("dataset must be an annual time series: program terminated. Check
         manual for input file format requirements.")
  }
  # if gaps larger than 15% of record length, terminate analysis
  if (n1/n > 0.15) {
    # If missing values > 15% stop
    stop("missing values larger than 15% of time series: program terminated. Check
         manual for input file format requirements.")
  }
  w <- y[, 2]  # check that max sequence of NA's is not greater than 5%
  w[is.na(w)] <- 0  # replace NA's with 0
  w <- rle(w)  # sequences (runs of equal values)
  w1 <- max(w$lengths[w$values == 0])  # max sequence of NA's (now 0)
  if (w1/n > 0.05) {
  # Continuous sequence of missing values >5% - gaps too large
    stop("maximum gap lengths too large (> 5% of time series length): program
         terminated. Check manual for input file format requirements.")
  }
  if (is.na(y[1, 2]) == TRUE | is.na(y[n, 2]) == TRUE) {
    # If time series starts or ends with missing values terminate
    stop("cannot fill missing values at start or end of time series : program
         terminated.")
  }
  # -----------------------------------------------------------------
    # if 'fillgaps' not entered set default
  if (fillgaps == 1 | fillgaps == 2 | fillgaps == 3) {
    fillgaps <- fillgaps
  } else {
    # If 'fillgaps' not entered or entered outside range
    print("default fillgaps setting applied")
    fillgaps = 1
  }
  # -----------------------------------------------------------------
  # If station name not entered
  if (station_name == " ") {
    station_name <- " "  # empty
    lab5 <- paste("Station Name not entered")
    print("no station name entered, field left blank")
  } else {
    station_name <- station_name
    lab5 <- station_name
  }
  # -----------------------------------------------------------------
  # where there are no gaps in the record
  # or if fillgaps argument entered when not needed
  if ((any(is.na(y[, 2])) == FALSE) && (fillgaps == 1 | fillgaps == 2 |
                                        fillgaps == 3)) {
    ts1 <- ts1
    lab1 <- paste("zero")  # number of filled gaps
    lab2 <- paste("NA")  # gap filling method
    print("no missing data in time series, gap-filling not required")
    TSfill <- vector(mode = "numeric", length = n)
    TSfill[TSfill == 0] <- NA
    p <- 1  # conditioning parameter (no gaps in record)
  }
  # -----------------------------------------------------------------
  # fill gaps using SSA (Package 'Rssa')
  if ((fillgaps == 1) && (p == 0)) {
    print("gapfilling using iterative SSA setting applied")
    if (requireNamespace("Rssa", quietly = TRUE)) {
      Rssa::grouping.auto
      Rssa::igapfill
      Rssa::reconstruct
      Rssa::ssa
    }
    if (requireNamespace("tseries", quietly = TRUE)) {
      tseries::na.remove
    }
    lab1 <- n1  # number of filled gaps
    lab2 <- paste("SSA")  # gap filling method
    shape <- Rssa::ssa(ts1, L = L1, kind = "1d-ssa", svd.method = "auto")  # SSA shape
    kk <- sum(shape$sigma > 0.1)
    recon <- Rssa::reconstruct(shape, groups = c(1:kk))
    lt <- NULL
    reps <- kk
    for (i in 1:reps) {
      f <- tseries::na.remove(recon[[i]])
      ss <- stats::spec.pgram(f, spans = 3, detrend = FALSE, log = "no", plot = FALSE)
      freq1 <- ss$freq[ss$spec == max(ss$spec)]  # max freq from max spectogram
      if (freq1 <= 0.2) {
        newel <- i
        lt <- c(lt, newel)
      }
    }
    newTS <- Rssa::igapfill(shape, groups = list(lt))
    ts1 <- newTS
    if (requireNamespace("zoo", quietly = TRUE)) {
      zoo::na.approx
      zoo::na.spline
    }
    TSfill <- newTS[, 1]  # extract TSfill values
  }
  # -----------------------------------------------------------------
  # fill gaps using linear interpolation from Package 'zoo')
  if ((fillgaps == 2) && (p == 0)) {
    print("gapfilling using linear interpolation setting applied")
    if (requireNamespace("zoo", quietly = TRUE)) {
      zoo::na.approx
      zoo::na.spline
    }
    lab1 <- n1  # number of filled gaps
    lab2 <- paste("Linear Interpolation")  # gap filling method
    TSfill <- zoo::na.approx(ts1, na.rm = FALSE)
    ts1 <- TSfill
    TSfill <- TSfill[, 1]  # extract TSfill values
  }
  # -----------------------------------------------------------------
  # fill gaps using cubic spline interpolation from Package 'zoo')
  if ((fillgaps == 3) && (p == 0)) {
    print("gapfilling using spline interpolation setting applied")
    if (requireNamespace("zoo", quietly = TRUE)) {
      zoo::na.approx
      zoo::na.spline
    }
    lab1 <- n1  # number of filled gaps
    lab2 <- paste("Spline Interpolation")  # gap filling method
    TSfill <- zoo::na.spline(ts1, na.rm = FALSE)
    ts1 <- TSfill
    TSfill <- TSfill[, 1]  # extract TSfill values
  }
  # -----------------------------------------------------------------
  # trend isolation using SSA in package (Rssa)
  if (requireNamespace("Rssa", quietly = TRUE)) {
    Rssa::grouping.auto
    Rssa::igapfill
    Rssa::reconstruct
    Rssa::ssa
  }
  model.ssa <- Rssa::ssa(ts1, kind = "1d-ssa", svd.method = "auto")  # default 1D-ssa
  # Auto-detection of trend components from ssa decomposition low frequency interval
  # [0, 0.01], threshold contribution = 0.75, 30 leading components
  ssa.auto <- Rssa::grouping.auto(model.ssa, base = "series", groups = list(1:30),
                            freq.bins = list(0.01), threshold = 0.75)
  recon <- Rssa::reconstruct(model.ssa, groups = ssa.auto)
  trend <- recon$F1
  # -----------------------------------------------------------------
  # fit cubic smooting spline to trend to determine velocity and acceleration
  k <- round((length(ts1)/8), 0)  # DoF for smooth spline (1 every 8 years)
  ss <- stats::smooth.spline(trend, df = k)
  ssVEL <- stats::predict(ss, type = "response", deriv = 1)
  ssVEL <- ssVEL$y  # Velocity vector
  ssACC <- stats::predict(ss, type = "response", deriv = 2)
  ssACC <- ssACC$y  # Acceleration vector
  # -----------------------------------------------------------------
  # determine standard errors
  # remove auto-correlation from residuals
  res <- trend - ts1  # correlated residuals
  if (requireNamespace("forecast", quietly = TRUE)) {
    forecast::auto.arima
  }
  aa <- forecast::auto.arima(res, seasonal = FALSE)  # auto.arima fit non-season
  noncor <- aa$residuals  # time series of non-correlated residuals
  # -----------------------------------------------------------------
  # detect changepoints in variance of non-correlated residuals
  if (requireNamespace("changepoint", quietly = TRUE)) {
    changepoint::cpt.var
    changepoint::cpts
  }
  cpt <- changepoint::cpt.var(noncor, minseglen = 15)  # min segment 15 years (default)
  ff <- changepoint::cpts(cpt)  # point in time series at which changepoint is detected
  # -----------------------------------------------------------------
  # bootsrapping routines
  TS <- vector(mode = "list", length = iter)
  VEL <- vector(mode = "list", length = iter)
  ACC <- vector(mode = "list", length = iter)
  # -----------------------------------------------------------------
  # no changepoints in the variance of residuals
  if (length(ff) == 0) {
    lab3 <- paste("zero")  # number of changepoints
    lab4 <- paste("NA")  # changepoint year
    for (i in 1:iter) {
      samp <- sample(noncor, replace = TRUE)
      newTS <- trend + samp
      model2.ssa <- Rssa::ssa(newTS, kind = "1d-ssa", svd.method = "auto")  # 1D-ssa
      ssa2.auto <- Rssa::grouping.auto(model2.ssa, base = "series",
                                       groups = list(1:30), freq.bins = list(0.01),
                                       threshold = 0.75)
      recon2 <- Rssa::reconstruct(model2.ssa, groups = ssa2.auto)
      trend2 <- recon2$F1
      TS[[i]] <- trend2
      ss2 <- stats::smooth.spline(trend2, df = k)
      ssVEL2 <- stats::predict(ss2, type = "response", deriv = 1)
      VEL[[i]] <- ssVEL2$y  ## Isolate Velocity Vector
      ssACC2 <- stats::predict(ss2, type = "response", deriv = 2)
      ACC[[i]] <- ssACC2$y  ## Isolate acceleration Vector
    }
  }
  # -----------------------------------------------------------------
  # 1 changepoint in variance of residuals detected
  if (length(ff) == 1) {
    lab3 <- 1  # number of changepoints
    lab4 <- y[ff, 1]  # changepoint year
    for (i in 1:iter) {
      q1 <- noncor[1:ff - 1]
      r1 <- sample(q1, replace = TRUE)
      q2 <- noncor[ff:n]
      r2 <- sample(q2, replace = TRUE)
      samp <- as.numeric(samp <- paste(c(r1, r2)))
      newTS <- trend + samp
      model2.ssa <- Rssa::ssa(newTS, kind = "1d-ssa", svd.method = "auto")  # 1D-ssa
      ssa2.auto <- Rssa::grouping.auto(model2.ssa, base = "series",
                                       groups = list(1:30), freq.bins = list(0.01),
                                       threshold = 0.75)
      recon2 <- Rssa::reconstruct(model2.ssa, groups = ssa2.auto)
      trend2 <- recon2$F1
      TS[[i]] <- trend2
      ss2 <- stats::smooth.spline(trend2, df = k)
      ssVEL2 <- stats::predict(ss2, type = "response", deriv = 1)
      VEL[[i]] <- ssVEL2$y  ## Isolate Velocity Vector
      ssACC2 <- stats::predict(ss2, type = "response", deriv = 2)
      ACC[[i]] <- ssACC2$y  ## Isolate acceleration Vector
    }
  }
  # -----------------------------------------------------------------
  # standard deviations
  TSboot <- as.data.frame(do.call(cbind, TS))
  TSsd <- apply(TSboot[, c(1:iter)], 1, sd)
  VELboot <- as.data.frame(do.call(cbind, VEL))
  VELsd <- apply(VELboot[, c(1:iter)], 1, sd)
  ACCboot <- as.data.frame(do.call(cbind, ACC))
  ACCsd <- apply(ACCboot[, c(1:iter)], 1, sd)
  # -----------------------------------------------------------------
  # summary dataframe
  summ <- as.data.frame(cbind(y[, 1], y[, 2], as.vector(trend), TSsd, ssVEL, VELsd,
                              ssACC, ACCsd, noncor, TSfill))
  colnames(summ) <- c("Year", "MSL", "Trend", "TrendSD", "Vel", "VelSD", "Acc",
                      "AccSD", "Resids", "FilledTS")
  st <- c(1:3, n - 2, n - 1, n)
  summ$Acc[st] <- NA
  # NA's at first and last 3 acceleration values (not accurate with spline)
  summ$AccSD[st] <- NA
  # -----------------------------------------------------------------
  # summary object
  summ2 <- summ[4:(n - 3), ]  # non NA portion of summ
  object <- NULL
  object$Station.Name <- lab5
  object$Summary <- summ
  object$Velocity <- list(Peak.mm.yr = max(summ$Vel),
                          Year = with(summ, Year[Vel == max(Vel)]))
  object$Acceleration <- list(Peak.mm.yr.yr = max(summ2$Acc),
                              Year = with(summ2, Year[Acc == max(summ2$Acc)]))
  object$Record.Length <- list(Start = y[1, 1], End = y[n, 1], Years = n)
  object$Fillgaps <- list(Gaps = lab1, Method = lab2)
  object$Bootstrapping.Iterations <- iter
  object$Changepoints <- list(Number = lab3, Year = lab4)
  class(object) <- "msl.trend"
  # -----------------------------------------------------------------
  # plot routines
  # plot option not entered or entered outside range
  if (plot == TRUE | plot == FALSE) {
    plot <- plot
  } else {
    print("default TYPE setting applied")
    plot <- TRUE
  }
  if (summ$Trend[n] - summ$Trend[1] < 0) {
    # check slope of graph for locating chart 1 legends
    laba <- paste("bottomleft")  # Chart key to bottomleft
  } else {
    # default setting
    laba <- paste("bottomright")  # Chart key to bottomright
  }
  # -----------------------------------------------------------------
  if (plot == TRUE) {
    # set plotting limits
    if (requireNamespace("plyr", quietly = TRUE)) {
      plyr::round_any
    }
    if (n < 100) {
      # setting up x-axis plotting parameters
      xtic = 10  # year ticks on x-axis
      xlo <- plyr::round_any(min(summ[, 1]), 10, floor)
      xhi <- plyr::round_any(max(summ[, 1]), 10, ceiling)
    } else {
      # default
      xtic = 20
      xlo <- plyr::round_any(min(summ[, 1]), 20, floor)
      xhi <- plyr::round_any(max(summ[, 1]), 20, ceiling)
    }
    xlim = c(xlo, xhi)
    # -----------------------------------------------------------------
    # set y-axis scale
    if (p == 0) {
      # gap filling routine required
      M1 <- max(max(summ[, 2], na.rm = TRUE), max(summ[, 10]))
      M2 <- min(min(summ[, 2], na.rm = TRUE), min(summ[, 10]))
    }
    if (p == 1) {
      # no gaps in record
      M1 <- max(summ[, 2])
      M2 <- min(summ[, 2])
    }
    ylen <- M1 - M2
    if (ylen < 200) {
      # setting up y-axis plotting parameters
      ytic = 20  # year ticks on y-axis
      ylo <- plyr::round_any(M2, 20, floor)
      yhi <- plyr::round_any(M1, 20, ceiling)
    } else {
      # default
      ytic = 50
      ylo <- plyr::round_any(M2, 50, floor)
      yhi <- plyr::round_any(M1, 50, ceiling)
    }
    ylim = c(ylo, yhi)
    # -----------------------------------------------------------------
    opar <- graphics::par(no.readonly = TRUE)  # capture current settings
    graphics::par(mar = c(3.6, 5.1, 2, 0.5), las = 1)
    graphics::plot(summ$Year, summ$MSL, type = "l", xaxt = "n", yaxt = "n", lty = 0,
         xlim = xlim, ylim = ylim, xlab = " ", ylab = " ", main = station_name)
    graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
         col = "lightcyan1")
    graphics::title(ylab = "Relative Mean Sea Level (mm)", font.lab = 2, line = 3.5)
    graphics::title(xlab = "Year", font.lab = 2, line = 2.5)
    graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
    graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
    # -----------------------------------------------------------------
    if (p == 0) {
      # gap filling routine required
      graphics::legend(laba, bg = "white", legend = c("KEY", "Annual Average Data",
                                            "MSL Trend", "95% Confidence Intervals",
                                            "Gap Filling"),
             text.font = c(2, 1, 1, 1, 1), lty = c(0, 1, 1, 2, 1),
             lwd = c(1, 1, 3, 1, 1), col = c("black", "darkgrey", "black",
                                             "black", "red"),
             cex = c(0.8, 0.8, 0.8, 0.8, 0.8))
      graphics::lines(summ$Year, summ$FilledTS, col = "red")
      graphics::lines(summ$Year, summ$MSL, col = "darkgrey")
      graphics::lines(summ$Year, summ$Trend, lwd = 2)
      graphics::lines(summ$Year, summ$Trend + 1.96 * summ$TrendSD, lty = 2)
      graphics::lines(summ$Year, summ$Trend - 1.96 * summ$TrendSD, lty = 2)
    }
    # -----------------------------------------------------------------
    if (p == 1) {
      # no gaps in record
      graphics::legend(laba, bg = "white", legend = c("KEY", "Annual Average Data",
                                            "MSL Trend", "95% Confidence Intervals"),
             text.font = c(2, 1, 1, 1, 1), lty = c(0, 1, 1, 2), lwd = c(1, 1, 3, 1),
             col = c("black", "darkgrey", "black", "black"),
             cex = c(0.8, 0.8, 0.8, 0.8))
      graphics::lines(summ$Year, summ$MSL, col = "darkgrey")
      graphics::lines(summ$Year, summ$Trend, lwd = 2)
      graphics::lines(summ$Year, summ$Trend + 1.96 * summ$TrendSD, lty = 2)
      graphics::lines(summ$Year, summ$Trend - 1.96 * summ$TrendSD, lty = 2)
    }
  }
  # -----------------------------------------------------------------
  # no plotting
  if (plot == FALSE) {
    print("no plotting required")
  }
  return(object)
}

#' Screen plotting options.
#'
#' @param x object of class \dQuote{msl.trend} (see \code{\link{msl.trend}} and
#' \code{\link{s}}) or \dQuote{msl.forecast} (see \code{\link{msl.forecast}} and
#' \code{\link{t}}).
#'
#' @param type numeric, enables a user defined input to select the type of chart
#' to be plotted. The default setting (type = 1) provides 3 charts in the same
#' plot area with the time series in the top panel, instantaneous velocity in
#' the middle panel and instantaneous acceleration in the bottom panel. The
#' alternatives (2, 3 and 4) are single panel plots of time series, instantaneous
#' velocity and instantaneous acceleration, respectively.
#'
#' @param ci numeric, enables a user defined input to select the type of
#' confidence interval to be displayed on the plots. The default setting (ci = 1)
#' corresponds to a 95\% confidence interval whilst ci=2 provides a 99\%
#' confidence interval.
#'
#' @details This routine provides a range of screen plotting options for both
#' \dQuote{msl.trend} (see \code{\link{msl.trend}}) and \dQuote{msl.forecast}
#' (see \code{\link{msl.forecast}}) objects. The same range of alternative pdf
#' plotting options are available via \code{\link{msl.pdf}}.
#'
#' @seealso \code{\link{msl.trend}}, \code{\link{msl.forecast}},
#' \code{\link{msl.pdf}}, \code{\link{Balt}}, \code{\link{s}}, \code{\link{t}}
#'
#' @examples
#' # -------------------------------------------------------------------------
#' # Isolate trend from Baltimore record, filling gaps with spline interpolation,
#' # 500 iterations and adding 1000 mm of slr to 2100. Use raw 'Balt.csv' data file.
#' # Note: ordinarily user would call 'File.csv' direct from working directory
#' # using the following sample code:
#' # s <- msl.trend('Balt.csv', fillgaps = 3, iter = 500, 'BALTIMORE, USA')
#' # t <- msl.forecast(s, slr = 1000)
#' # -------------------------------------------------------------------------
#'
#' data(s) # msl.trend object from above-mentioned example
#' data(t) # msl.forecast object from above-mentioned example
#' msl.plot(s) # default screen plot output, 3 panels, 95% confidence intervals
#' msl.plot(s, type = 2) # plot time series, 95% confidence intervals
#' msl.plot(s, type = 3) # plot instantaneous velocity, 95% confidence intervals
#' msl.plot(s, type = 4, ci = 2) # plot acceleration, 99% confidence intervals
#' msl.plot(t) # default screen plot output, 3 panels, 95% confidence intervals
#' msl.plot(t, type = 2) # plot time series, 95% confidence intervals
#' msl.plot(t, type = 3) # plot instantaneous velocity, 95% confidence intervals
#' msl.plot(t, type = 4, ci = 2) # plot acceleration, 99% confidence intervals
#'
#' @export
msl.plot <- function(x, type = 1, ci = 1) {
    # -----------------------------------------------------------------
    object <- x
    summ <- object$Summary
    grDevices::graphics.off()  # close active graphics windows
    # -----------------------------------------------------------------
    # object is a msl.trend or msl.forecast output dataframe
    # type alternatives = c(1,2,3,4) type = 1 (default)
    # default (3 plots: time series, velocity, acceleration)
    # type = 2 (single plot: time series) type = 3 (single plot: velocity)
    # type = 4 (single plot: acceleration)
    # ci alternatives = c('1','2') ci = 1 (default) (confidence interval = 95%)
    # ci = 2 (confidence interval = 99%)
    # -----------------------------------------------------------------
    # If not msl.trend or msl.forecast object
    if (class(object) == "msl.trend" | class(object) == "msl.forecast") {
        class(object) <- class(object)
    } else {
        stop("object is not an msl.trend or msl.forecast object: plotting terminated")
    }
    # -----------------------------------------------------------------
    # check type of object to direct to specific plots
    if (class(object) == "msl.trend") {
      tp <- 0
    } else {
      tp <- 1
      labc <- paste0("Projected Sea Level Rise = ",
                     object$Projected.SLR$Slr.mm, " mm")  # legend with slr
      labd <- paste0("Projection (SLR = ",
                     object$Projected.SLR$Slr.mm, " mm)")  # legend with slr
    }
    # -----------------------------------------------------------------
    # check if original time series contains missing values
    if (tp == 0) { # msl.trend object
      if (any(is.na(summ$MSL)) == TRUE) {
        p <- 0
      } else {
        p <- 1
      }
    }
    if (tp == 1) { # msl.forecast object
      if (any(is.na(summ$MSL[1:object$Historical.Record$Years])) == TRUE) {
        p <- 0
      } else {
        p <- 1
      }
    }
    # -----------------------------------------------------------------
    # specific settings, portions of object
    if (tp == 0) {
        n <- length(summ[, 1])
        n2 <- n - 3
        # dataframe without NA's for acceleration plot with msl.trend object
        summ2 <- summ[4:n2, ]
    } else {
        n <- length(summ[, 1])
        # dataframe without NA's for acceleration plot with msl.forecast object
        summ2 <- summ[4:n, ]
        # different plotting colours for historical
        summ3 <- summ[summ$Year <= object$Historical.Record$End, ]
        # different plotting colours for forecast
        summ4 <- summ[summ$Year >= object$Historical.Record$End, ]
    }
    # -----------------------------------------------------------------
    # Station name
    if (object$Station.Name == "Station Name not entered") {
        object$Station.Name <- NULL
    } else {
        # default setting is entered Station Name
        object$Station.Name <- object$Station.Name
    }
    # -----------------------------------------------------------------
    # type not entered or entered outside range
    if (type == 1 | type == 2 | type == 3 | type == 4) {
        type <- type
    } else {
        print("default TYPE setting applied")
        type <- 1
    }
    if (ci == 1 | ci == 2) {
        # If ci not entered or entered outside range
        ci <- ci
    } else {
        print("default CONFIDENCE INTERVAL setting applied")
        ci <- 1
    }
    # -----------------------------------------------------------------
    # Confidence interval input
    if (ci == 2) {
        ci = 2.575  # multiplication factor for 99% CI
        lab1 <- paste("99% Confidence Interval")
    } else {
        # default setting is 95% CI
        ci = 1.96
        lab1 <- paste("95% Confidence Interval")
    }
    # -----------------------------------------------------------------
    # check slope of graph for locating chart 1 legends
    if (summ$Trend[n] - summ$Trend[1] < 0) {
        laba = paste("topright")  # Chart description to topright
        labb = paste("bottomleft")  # Chart key to bottomleft
    } else {
        # default setting
        laba = paste("topleft")  # Chart description to topleft
        labb = paste("bottomright")  # Chart key to bottomright
    }
    # -----------------------------------------------------------------
    if (requireNamespace("plyr", quietly = TRUE)) {
      plyr::round_any
    }
    if (n < 100) {
        # setting up x-axis plotting parameters
        xtic = 10  # year ticks on x-axis
        xlo <- plyr::round_any(min(summ[, 1]), 10, floor)
        xhi <- plyr::round_any(max(summ[, 1]), 10, ceiling)
    } else {
        # default
        xtic = 20
        xlo <- plyr::round_any(min(summ[, 1]), 20, floor)
        xhi <- plyr::round_any(max(summ[, 1]), 20, ceiling)
    }
    xlim = c(xlo, xhi)
    # -----------------------------------------------------------------
    # plotting routines for msl.trend objects
    # -----------------------------------------------------------------
    # plot 3 charts: time series, velocity and acceleration
    if ((type == 1) && (tp == 0)) {
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mfrow = c(3, 1), las = 1)
        graphics::par(mar = c(0, 5.1, 2, 0.6), las = 1)  # chart 1 (time series)
        # -----------------------------------------------------------------
        # set y-axis scale
        if (p == 0) {
          # gap filling routine required
          M1 <- max(max(summ[, 2], na.rm = TRUE), max(summ[, 10]))
          M2 <- min(min(summ[, 2], na.rm = TRUE), min(summ[, 10]))
        }
        if (p == 1) {
          # no gaps in record
          M1 <- max(summ[, 2])
          M2 <- min(summ[, 2])
        }
        ylen <- M1 - M2
        if (ylen < 200) {
          # setting up y-axis plotting parameters
          ytic = 50  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 50, floor)
          yhi <- plyr::round_any(M1, 50, ceiling)
        } else {
          # default
          ytic = 100
          ylo <- plyr::round_any(M2, 100, floor)
          yhi <- plyr::round_any(M1, 100, ceiling)
        }
        ylim = c(ylo, yhi)
        # -----------------------------------------------------------------
        graphics::plot(summ[, 1], summ[, 2], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 2)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])),
                c((summ[, 3] + ci * summ[, 4]), rev((summ[, 3] - ci * summ[, 4]))),
                col = "azure3", border = NA)
        graphics::title(ylab = "Millimetres", font.lab = 2, cex.lab = 1.2, line = 3.5)
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        if (p == 0) {
            # gap filling routine required
            graphics::lines(summ[, 1], summ[, 10], col = "red")
            graphics::lines(summ[, 1], summ[, 2])
            graphics::lines(summ[, 1], summ[, 3], lwd = 2)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                   inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
            graphics::legend(labb, bg = "white",
                   legend = c("KEY", "Annual Average Data", "Gap Filling",
                              "MSL Trend", lab1, "Peak Rate"),
                   text.font = c(2, 1, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 1, 3),
                   lwd = c(1, 1, 1, 3, 8, 3),
                   col = c("black", "black", "red", "black", "azure3", "blue"),
                   cex = c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9))
        }
        if (p == 1) {
            # no gap filling routine required
            graphics::lines(summ[, 1], summ[, 2])
            graphics::lines(summ[, 1], summ[, 3], lwd = 2)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                             inset = c(-0.02, -0.01), bty = "n", text.font = 2,
                             cex = 1.2)
            graphics::legend(labb, bg = "white",
                   legend = c("KEY", "Annual Average Data", "MSL Trend", lab1,
                              "Peak Rate"), text.font = c(2, 1, 1, 1, 1),
                   lty = c(0, 1, 1, 1, 3), lwd = c(1, 1, 3, 8, 3),
                   col = c("black", "black", "black", "azure3", "blue"),
                   cex = c(0.9, 0.9, 0.9, 0.9, 0.9))
        }
        graphics::par(mar = c(0, 5.1, 0, 0.6), las = 1)  # chart 2 (velocity)
        ylen <- max(summ[, 5] + ci * summ[, 6]) - min(summ[, 5] - ci * summ[, 6])
        if (ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 1 & ylen <= 2) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        if (ylen > 2 & ylen <= 10) {
            # setting up y-axis plotting parameters
            ytic = 1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   1, ceiling)
        }
        if (ylen > 10 & ylen <= 20) {
            # setting up y-axis plotting parameters
            ytic = 2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   2, ceiling)
        }
        if (ylen > 20) {
            # setting up y-axis plotting parameters
            ytic = 5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   5, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   5, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ[, 1], summ[, 5], type = "line", lty = 0, xaxt = "n",
                       yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 5] + ci * summ[, 6]),
                                                rev((summ[, 5] - ci * summ[, 6]))),
                col = "azure3", border = NA)
        graphics::lines(summ[, 1], summ[, 5], lwd = 2)
        graphics::abline(h = 0, lty = 2)
        graphics::abline(h = max(summ[, 5]), lty = 3, col = "blue", lwd = 3)
        graphics::legend("topleft", legend = "Instantaneous Velocity",
                         inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
        graphics::title(ylab = "Millimetres/year", font.lab = 2, cex.lab = 1.2,
                        line = 3.5)
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(mar = c(3.6, 5.1, 0, 0.6), las = 1)  # chart 3 (acceleration)
        ylen <- max(summ2[, 7] + ci * summ2[, 8]) - min(summ2[, 7] - ci * summ2[, 8])
        if (ylen <= 0.1) {
            # setting up y-axis plotting parameters
            ytic = 0.02  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.02, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.02, ceiling)
        }
        if (ylen > 0.1 & ylen <= 0.2) {
            # setting up y-axis plotting parameters
            ytic = 0.05  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.05, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.05, ceiling)
        }
        if (ylen > 0.2 & ylen <= 0.5) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 0.5 & ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        if (ylen > 1) {
            # setting up y-axis plotting parameters
            ytic = 0.5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.5, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.5, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ2[, 1], summ2[, 7], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ2[, 1],
                  rev(summ2[, 1])), c((summ2[, 7] + ci * summ2[, 8]),
                                      rev((summ2[, 7] - ci * summ2[, 8]))),
                col = "azure3", border = NA)
        graphics::lines(summ2[, 1], summ2[, 7], lwd = 2)
        graphics::abline(h = 0, lty = 2)
        graphics::abline(h = max(summ2[, 7]), lty = 3, col = "blue", lwd = 3)
        graphics::legend("topleft", legend = "Instantaneous Acceleration",
               inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
        graphics::title(ylab = "Millimetres/year/year", font.lab = 2,
                        cex.lab = 1.2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
    }
    # -----------------------------------------------------------------
    # plot time series only
    if ((type == 2) && (tp == 0)) {
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2.1, 0.6), las = 1)
        # -----------------------------------------------------------------
        # set y-axis scale
        if (p == 0) {
          # gap filling routine required
          M1 <- max(max(summ[, 2], na.rm = TRUE), max(summ[, 10]))
          M2 <- min(min(summ[, 2], na.rm = TRUE), min(summ[, 10]))
        }
        if (p == 1) {
          # no gaps in record
          M1 <- max(summ[, 2])
          M2 <- min(summ[, 2])
        }
        ylen <- M1 - M2
        if (ylen < 200) {
          # setting up y-axis plotting parameters
          ytic = 20  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 20, floor)
          yhi <- plyr::round_any(M1, 20, ceiling)
        } else {
          # default
          ytic = 50
          ylo <- plyr::round_any(M2, 50, floor)
          yhi <- plyr::round_any(M1, 50, ceiling)
        }
        ylim = c(ylo, yhi)
        # -----------------------------------------------------------------
        graphics::plot(summ[, 1], summ[, 2], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 1.7)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 3] + ci * summ[, 4]),
                                                rev((summ[, 3] - ci * summ[, 4]))),
                col = "azure3", border = NA)
        graphics::title(ylab = "Millimetres", font.lab = 2, cex.lab = 1.2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        if (p == 0) {
            # gap filling routine required
            graphics::lines(summ[, 1], summ[, 10], col = "red")
            graphics::lines(summ[, 1], summ[, 2])
            graphics::lines(summ[, 1], summ[, 3], lwd = 2)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                             inset = c(-0.04, -0.01), bty = "n", text.font = 2,
                             cex = 1.3)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "Gap Filling",
                                        "MSL Trend", lab1),
                   text.font = c(2, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 1),
                   lwd = c(1, 1, 1, 3, 8),
                   col = c("black", "black", "red", "black", "azure3"),
                   cex = c(0.8, 0.8, 0.8, 0.8, 0.8))
        }
        if (p == 1) {
            # no gap filling routine required
            graphics::lines(summ[, 1], summ[, 2])
            graphics::lines(summ[, 1], summ[, 3], lwd = 2)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                   inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "MSL Trend",
                                        lab1),
                   text.font = c(2, 1, 1, 1), lty = c(0, 1, 1, 1),
                   lwd = c(1, 1, 3, 8), col = c("black", "black", "black", "azure3"),
                   cex = c(0.8, 0.8, 0.8, 0.8))
        }
        graphics::par(opar)  # restore original settings
    }
    # -----------------------------------------------------------------
    # plot velocity series only
    if ((type == 3) && (tp == 0)) {
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2.1, 0.6), las = 1)
        ylen <- max(summ[, 5] + ci * summ[, 6]) - min(summ[, 5] - ci * summ[, 6])
        if (ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 1 & ylen <= 2) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        if (ylen > 2 & ylen <= 5) {
            # setting up y-axis plotting parameters
            ytic = 0.5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.5, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.5, ceiling)
        }
        if (ylen > 5 & ylen <= 10) {
            # setting up y-axis plotting parameters
            ytic = 1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   1, ceiling)
        }
        if (ylen > 10) {
            # setting up y-axis plotting parameters
            ytic = 2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   2, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ[, 1], summ[, 5], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ",  ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 1.7)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 5] + ci * summ[, 6]),
                                                rev((summ[, 5] - ci * summ[, 6]))),
                col = "azure3", border = NA)
        graphics::lines(summ[, 1], summ[, 5], lwd = 2)
        graphics::abline(h = 0, lty = 2)
        graphics::abline(h = max(summ[, 5]), lty = 3, col = "blue", lwd = 3)
        graphics::legend("topleft", legend = "Instantaneous Velocity",
               inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
        graphics::legend("bottomright", bg = "white", legend = c("KEY", "MSL Velocity",
                                                       "Peak Velocity", lab1),
               text.font = c(2, 1, 1, 1), lty = c(0, 1, 3, 1), lwd = c(1, 3, 3, 8),
               col = c("black", "black", "blue", "azure3"),
               cex = c(0.8, 0.8, 0.8, 0.8))
        graphics::title(ylab = "Millimetres/year", font.lab = 2, cex.lab = 1.2,
                        line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
    }
    # -----------------------------------------------------------------
    # plot acceleration series only
    if ((type == 4) && (tp == 0)) {
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2.1, 0.6), las = 1)
        ylen <- max(summ2[, 7] + ci * summ2[, 8]) - min(summ2[, 7] - ci * summ2[, 8])
        if (ylen <= 0.1) {
            # setting up y-axis plotting parameters
            ytic = 0.01  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.01, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.01, ceiling)
        }
        if (ylen > 0.1 & ylen <= 0.2) {
            # setting up y-axis plotting parameters
            ytic = 0.02  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.02, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.02, ceiling)
        }
        if (ylen > 0.2 & ylen <= 0.5) {
            # setting up y-axis plotting parameters
            ytic = 0.05  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.05, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.05, ceiling)
        }
        if (ylen > 0.5 & ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 1) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ2[, 1], summ2[, 7], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 1.7)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ2[, 1], rev(summ2[, 1])),
              c((summ2[, 7] + ci * summ2[, 8]), rev((summ2[, 7] - ci * summ2[, 8]))),
                col = "azure3", border = NA)
        graphics::lines(summ2[, 1], summ2[, 7], lwd = 2)
        graphics::abline(h = 0, lty = 2)
        graphics::abline(h = max(summ2[, 7]), lty = 3, col = "blue", lwd = 3)
        graphics::legend("topleft", legend = "Instantaneous Acceleration",
               inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
        graphics::legend("bottomright", bg = "white",
                         legend = c("KEY", "MSL Acceleration", "Peak Acceleration",
                                    lab1),
               text.font = c(2, 1, 1, 1), lty = c(0, 1, 3, 1), lwd = c(1, 3, 3, 8),
               col = c("black", "black", "blue", "azure3"),
               cex = c(0.8, 0.8, 0.8, 0.8))
        graphics::title(ylab = "Millimetres/year/year", font.lab = 2,
                        cex.lab = 1.2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
    }
    # -----------------------------------------------------------------
    # plotting routines for msl.forecast objects
    # -----------------------------------------------------------------
    # plot 3 charts: time series, velocity and acceleration
    if ((type == 1) && (tp == 1)) {
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mfrow = c(3, 1), las = 1)
        graphics::par(mar = c(0, 5.1, 2, 0.6), las = 1)  # chart 1 (time series)
        # -----------------------------------------------------------------
        # set y-axis scale
        if (p == 0) {
          # gap filling routine required
          M1 <- max(max(summ$MSL, na.rm = TRUE),
                    max(summ$Trend + ci * summ$TrendSD),
                    max(summ$FilledTS, na.rm = TRUE))
          M2 <- min(min(summ$MSL, na.rm = TRUE),
                    min(summ$Trend - ci * summ$TrendSD),
                    min(summ$FilledTS, na.rm = TRUE))
        }
        if (p == 1) {
          # no gaps in record
          M1 <- max(max(summ$MSL, na.rm = TRUE), max(summ$Trend + ci * summ$TrendSD))
          M2 <- min(min(summ$MSL, na.rm = TRUE), min(summ$Trend - ci * summ$TrendSD))
        }
        ylen <- M1 - M2
        if (ylen <= 500) {
          # setting up y-axis plotting parameters
          ytic = 50  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 50, floor)
          yhi <- plyr::round_any(M1, 50, ceiling)
        }
        if (ylen > 500 & ylen <= 1000) {
          # setting up y-axis plotting parameters
          ytic = 100  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 100, floor)
          yhi <- plyr::round_any(M1, 100, ceiling)
        }
        if (ylen > 1000) {
          # setting up y-axis plotting parameters
          ytic = 200  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 200, floor)
          yhi <- plyr::round_any(M1, 200, ceiling)
        }
        ylim = c(ylo, yhi)
        # -----------------------------------------------------------------
        graphics::plot(summ$Year, summ$MSL, type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 2)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])),
                c((summ[, 3] + ci * summ[, 4]), rev((summ[, 3] - ci * summ[, 4]))),
                col = "azure3", border = NA)
        graphics::title(ylab = "Millimetres", font.lab = 2, cex.lab = 1.2, line = 3.5)
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        if (p == 0) {
            # gap filling routine required
            graphics::lines(summ[, 1], summ[, 9], col = "red")
            graphics::lines(summ[, 1], summ[, 2])
            graphics::lines(summ3[, 1], summ3[, 3], lwd = 2)
            graphics::lines(summ4[, 1], summ4[, 3], lwd = 2, col = "black", lty = 3)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                             inset = c(-0.02, -0.01), bty = "n", text.font = 2,
                             cex = 1.2)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "Gap Filling",
                                        "MSL Trend", lab1, "Projection"),
                   text.font = c(2, 1, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 1, 3),
                   lwd = c(1, 1, 1, 3, 8, 3), col = c("black", "black", "red",
                                                      "black", "azure3", "black"),
                   cex = c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9))
        }
        if (p == 1) {
            # no gap filling routine required
            graphics::lines(summ[, 1], summ[, 2])
            graphics::lines(summ3[, 1], summ3[, 3], lwd = 2)
            graphics::lines(summ4[, 1], summ4[, 3], lwd = 2, col = "black", lty = 3)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                             inset = c(-0.02, -0.01), bty = "n", text.font = 2,
                             cex = 1.2)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "MSL Trend",
                                        lab1, "Projection"),
                   text.font = c(2, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 3),
                   lwd = c(1, 1, 3, 8, 3), col = c("black", "black", "black",
                                                   "azure3", "blue"),
                   cex = c(0.9, 0.9, 0.9, 0.9, 0.9))
        }
        graphics::par(mar = c(0, 5.1, 0, 0.6), las = 1)  # chart 2 (velocity)
        ylen <- max(summ[, 5] + ci * summ[, 6]) - min(summ[, 5] - ci * summ[, 6])
        if (ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        if (ylen > 1 & ylen <= 2) {
            # setting up y-axis plotting parameters
            ytic = 0.5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.5, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.5, ceiling)
        }
        if (ylen > 2 & ylen <= 10) {
            # setting up y-axis plotting parameters
            ytic = 1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   1, ceiling)
        }
        if (ylen > 10 & ylen <= 20) {
            # setting up y-axis plotting parameters
            ytic = 2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   2, ceiling)
        }
        if (ylen > 20) {
            # setting up y-axis plotting parameters
            ytic = 5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   5, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   5, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ[, 1], summ[, 5], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 5] + ci * summ[, 6]),
                                                rev((summ[, 5] - ci * summ[, 6]))),
                col = "azure3", border = NA)
        graphics::lines(summ3[, 1], summ3[, 5], lwd = 2)
        graphics::lines(summ4[, 1], summ4[, 5], lwd = 2, col = "black", lty = 3)
        graphics::abline(h = 0, lty = 2)
        graphics::legend("topleft", legend = "Instantaneous Velocity",
                         inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
        graphics::title(ylab = "Millimetres/year", font.lab = 2, cex.lab = 1.2,
                        line = 3.5)
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(mar = c(3.6, 5.1, 0, 0.6), las = 1)  # chart 3 (acceleration)
        ylen <- max(summ2[, 7] + ci * summ2[, 8]) - min(summ2[, 7] - ci * summ2[, 8])
        if (ylen <= 0.1) {
            # setting up y-axis plotting parameters
            ytic = 0.02  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.02, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.02, ceiling)
        }
        if (ylen > 0.1 & ylen <= 0.2) {
            # setting up y-axis plotting parameters
            ytic = 0.05  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.05, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.05, ceiling)
        }
        if (ylen > 0.2 & ylen <= 0.5) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 0.5 & ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        if (ylen > 1) {
            # setting up y-axis plotting parameters
            ytic = 0.5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.5, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.5, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ2[, 1], summ2[, 7], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ2[, 1], rev(summ2[, 1])),
                          c((summ2[, 7] + ci * summ2[, 8]),
                            rev((summ2[, 7] - ci * summ2[, 8]))),
                col = "azure3", border = NA)
        graphics::lines(summ3[, 1], summ3[, 7], lwd = 2)
        graphics::lines(summ4[, 1], summ4[, 7], lwd = 2, col = "black", lty = 3)
        graphics::abline(h = 0, lty = 2)
        graphics::legend("topleft", legend = "Instantaneous Acceleration",
               inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
        graphics::legend("bottomleft", legend = labc, inset = c(-0.02, -0.01),
                         bty = "n", text.font = 2, cex = 1.2)
        graphics::title(ylab = "Millimetres/year/year", font.lab = 2, cex.lab = 1.2,
                        line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
    }
    # -----------------------------------------------------------------
    # plot time series only
    if ((type == 2) && (tp == 1)) {
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2.1, 0.6), las = 1)
        # -----------------------------------------------------------------
        # set y-axis scale
        if (p == 0) {
          # gap filling routine required
          M1 <- max(max(summ$MSL, na.rm = TRUE),
                    max(summ$Trend + ci * summ$TrendSD),
                    max(summ$FilledTS, na.rm = TRUE))
          M2 <- min(min(summ$MSL, na.rm = TRUE),
                    min(summ$Trend - ci * summ$TrendSD),
                    min(summ$FilledTS, na.rm = TRUE))
        }
        if (p == 1) {
          # no gaps in record
          M1 <- max(max(summ$MSL, na.rm = TRUE), max(summ$Trend + ci * summ$TrendSD))
          M2 <- min(min(summ$MSL, na.rm = TRUE), min(summ$Trend - ci * summ$TrendSD))
        }
        ylen <- M1 - M2
        if (ylen < 500) {
          # setting up y-axis plotting parameters
          ytic = 50  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 50, floor)
          yhi <- plyr::round_any(M1, 50, ceiling)
        }
        if (ylen > 500 & ylen <= 1000) {
          # setting up y-axis plotting parameters
          ytic = 100  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 100, floor)
          yhi <- plyr::round_any(M1, 100, ceiling)
        }
        if (ylen > 1000) {
          # setting up y-axis plotting parameters
          ytic = 200  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 200, floor)
          yhi <- plyr::round_any(M1, 200, ceiling)
        }
        ylim = c(ylo, yhi)
        # -----------------------------------------------------------------
        graphics::plot(summ[, 1], summ[, 2], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 1.7)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 3] + ci * summ[, 4]),
                                                rev((summ[, 3] - ci * summ[, 4]))),
                col = "azure3", border = NA)
        graphics::title(ylab = "Millimetres", font.lab = 2, cex.lab = 1.2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        if (p == 0) {
            # gap filling routine required
            graphics::lines(summ[, 1], summ[, 9], col = "red")
            graphics::lines(summ[, 1], summ[, 2])
            graphics::lines(summ3[, 1], summ3[, 3], lwd = 2)
            graphics::lines(summ4[, 1], summ4[, 3], lwd = 2, col = "black", lty = 3)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                   inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
            graphics::legend(labb, bg = "white",
                   legend = c("KEY", "Annual Average Data", "Gap Filling",
                              "MSL Trend", lab1, labd),
                   text.font = c(2, 1, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 1, 3),
                   lwd = c(1, 1, 1, 3, 8, 3), col = c("black", "black", "red",
                                                      "black", "azure3", "black"),
                   cex = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8))
        }
        if (p == 1) {
            # no gap filling routine required
            graphics::lines(summ[, 1], summ[, 2])
            graphics::lines(summ[, 1], summ[, 3], lwd = 2)
            graphics::lines(summ3[, 1], summ3[, 3], lwd = 2)
            graphics::lines(summ4[, 1], summ4[, 3], lwd = 2, col = "black", lty = 3)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                             inset = c(-0.04, -0.01), bty = "n", text.font = 2,
                             cex = 1.3)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "MSL Trend",
                                        lab1, "Projection"),
                   text.font = c(2, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 3),
                   lwd = c(1, 1, 3, 8, 3), col = c("black", "black", "black",
                                                   "azure3", "black"),
                   cex = c(0.8, 0.8, 0.8, 0.8, 0.8))
        }
        graphics::par(opar)  # restore original settings
    }
    # -----------------------------------------------------------------
    # plot velocity series only
    if ((type == 3) && (tp == 1)) {
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2.1, 0.6), las = 1)
        ylen <- max(summ[, 5] + ci * summ[, 6]) - min(summ[, 5] - ci * summ[, 6])
        if (ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 1 & ylen <= 2) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        if (ylen > 2 & ylen <= 5) {
            # setting up y-axis plotting parameters
            ytic = 0.5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.5, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.5, ceiling)
        }
        if (ylen > 5 & ylen <= 10) {
            # setting up y-axis plotting parameters
            ytic = 1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   1, ceiling)
        }
        if (ylen > 10) {
            # setting up y-axis plotting parameters
            ytic = 2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   2, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ[, 1], summ[, 5], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 1.7)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 5] + ci * summ[, 6]),
                                                rev((summ[, 5] - ci * summ[, 6]))),
                col = "azure3", border = NA)
        graphics::lines(summ3[, 1], summ3[, 5], lwd = 2)
        graphics::lines(summ4[, 1], summ4[, 5], lwd = 2, col = "black", lty = 3)
        graphics::abline(h = 0, lty = 2)
        graphics::legend("topleft", legend = "Instantaneous Velocity",
                         inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
        graphics::legend("bottomright", bg = "white",
               legend = c("KEY", "MSL Velocity", lab1, labd),
               text.font = c(2, 1, 1, 1), lty = c(0, 1, 1, 3), lwd = c(1, 3, 8, 3),
               col = c("black", "black", "azure3", "black"),
               cex = c(0.8, 0.8, 0.8, 0.8))
        graphics::title(ylab = "Millimetres/year", font.lab = 2, cex.lab = 1.2,
                        line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
    }
    # -----------------------------------------------------------------
    # plot acceleration series only
    if ((type == 4) && (tp == 1)) {
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2.1, 0.6), las = 1)
        ylen <- max(summ2[, 7] + ci * summ2[, 8]) - min(summ2[, 7] - ci * summ2[, 8])
        if (ylen <= 0.1) {
            # setting up y-axis plotting parameters
            ytic = 0.01  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.01, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.01, ceiling)
        }
        if (ylen > 0.1 & ylen <= 0.2) {
            # setting up y-axis plotting parameters
            ytic = 0.02  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.02, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.02, ceiling)
        }
        if (ylen > 0.2 & ylen <= 0.5) {
            # setting up y-axis plotting parameters
            ytic = 0.05  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.05, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.05, ceiling)
        }
        if (ylen > 0.5 & ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 1) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ2[, 1], summ2[, 7], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 1.7)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ2[, 1], rev(summ2[, 1])),
                          c((summ2[, 7] + ci * summ2[, 8]),
                            rev((summ2[, 7] - ci * summ2[, 8]))), col = "azure3",
                          border = NA)
        graphics::lines(summ3[, 1], summ3[, 7], lwd = 2)
        graphics::lines(summ4[, 1], summ4[, 7], lwd = 2, col = "black", lty = 3)
        graphics::abline(h = 0, lty = 2)
        graphics::legend("topleft", legend = "Instantaneous Acceleration",
               inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
        graphics::legend("bottomright", bg = "white",
               legend = c("KEY", "MSL Acceleration", lab1, labd),
               text.font = c(2, 1, 1, 1), lty = c(0, 1, 1, 3), lwd = c(1, 3, 8, 3),
               col = c("black", "black", "azure3", "black"),
               cex = c(0.8, 0.8, 0.8, 0.8))
        graphics::title(ylab = "Millimetres/year/year", font.lab = 2,
                        cex.lab = 1.2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
    }
}

#' Pdf plotting options.
#'
#' @param x object of class \dQuote{msl.trend} (see \code{\link{msl.trend}} and
#' \code{\link{s}}) or \dQuote{msl.forecast} (see \code{\link{msl.forecast}} and
#' \code{\link{t}}).
#'
#' @param file_name is a character string indicating the name of the pdf output
#' file. If this field is left blank the output file will be automatically saved
#' in the working directory under the default name "File1.pdf".
#'
#' @param type numeric, enables a user defined input to select the type of chart
#' to be plotted. The default setting (type = 1) provides 3 charts in the same
#' plot area with the time series in the top panel, instantaneous velocity in
#' the middle panel and instantaneous acceleration in the bottom panel. The
#' alternatives (2, 3 and 4) are single panel plots of time series, instantaneous
#' velocity and instantaneous acceleration, respectively.
#'
#' @param ci numeric, enables a user defined input to select the type of
#' confidence interval to be displayed on the plots. The default setting
#' (ci = 1) corresponds to a 95\% confidence interval whilst ci=2 provides a
#' 99\% confidence interval.
#'
#' @details This routine provides a range of pdf plotting options for both
#' \dQuote{msl.trend} (see \code{\link{msl.trend}}) and \dQuote{msl.forecast}
#' (see \code{\link{msl.forecast}}) objects. Three panel plots (type 1 or
#' default) are formatted with width = 16.54 inches and height = 20 inches.
#' Single panel plots (types 2, 3, 4) are formatted with width = 16.54 inches
#' and height = 15 inches. All plots are designed to be proportionally correct
#' when imported into documents and re-sized to the width of a standard A4 page.
#' The same range of alternative screen plotting options are available via
#' \code{\link{msl.plot}}.
#'
#' @seealso \code{\link{msl.trend}}, \code{\link{msl.forecast}},
#' \code{\link{msl.plot}}, \code{\link{Balt}}, \code{\link{s}}, \code{\link{t}}.
#'
#' @examples
#' # -------------------------------------------------------------------------
#' # Isolate trend from Baltimore record, filling gaps with spline interpolation,
#' # 500 iterations and adding 1000 mm of slr to 2100. Use raw 'Balt.csv' data file.
#' # Note: ordinarily user would call 'File.csv' direct from working directory
#' # using the following sample code:
#' # s <- msl.trend('Balt.csv', fillgaps = 3, iter = 500, 'BALTIMORE, USA')
#' # t <- msl.forecast(s, slr = 1000)
#' # -------------------------------------------------------------------------
#'
#' data(s) # msl.trend object from above-mentioned example
#' data(t) # msl.forecast object from above-mentioned example
#'
#' # default output, 3 panels, 95% confidence intervals.
#' msl.pdf(s)
#' # Check 'File1.pdf' in working directory
#'
#' # pdf plot time series, 95% confidence intervals.
#' msl.pdf(s, file_name = 'Series.pdf', type = 2)
#' # Check 'Series.pdf' file in working directory
#'
#' # pdf plot instantaneous velocity, 95% confidence intervals.
#' msl.pdf(s, file_name = 'Velocity.pdf', type = 3)
#' # Check 'Velocity.pdf' file in working directory
#'
#' # pdf plot instantaneous acceleration, 99% confidence intervals.
#' msl.pdf(s, file_name = 'Acceleration.pdf', type = 4, ci = 2)
#' # Check 'Acceleration.pdf' file in working directory
#'
#' # default output, 3 panels, 95% confidence intervals.
#' msl.pdf(t, file_name = 'Forecast.pdf')
#' # Check 'Forecast.pdf' file in working directory
#'
#' @export
msl.pdf <- function(x, file_name = " ", type = 1, ci = 1) {
    # -----------------------------------------------------------------
    object <- x
    summ <- object$Summary
    grDevices::graphics.off()  # close active graphics windows
    # -----------------------------------------------------------------
    # object is a msl.trend output dataframe
    # type alternatives = c(1,2,3,4) type = 1 (default)
    # default (3 plots: time series, velocity, acceleration)
    # type = 2 (single plot: time series) type = 3 (single plot: velocity)
    # type = 4 (single plot: acceleration)
    # ci alternatives = c('1','2') ci = 1 (default) (confidence interval = 95%)
    # ci = 2 (confidence interval = 99%)
    # -----------------------------------------------------------------
    # If not msl.trend or msl.forecast object
    if (class(object) == "msl.trend" | class(object) == "msl.forecast") {
      class(object) <- class(object)
    } else {
      stop("object is not an msl.trend or msl.forecast object: plotting terminated")
    }
    # -----------------------------------------------------------------
    # check type of object to direct to specific plots
    if (class(object) == "msl.trend") {
      tp <- 0
    } else {
      tp <- 1
      labc <- paste0("Projected Sea Level Rise = ",
                     object$Projected.SLR$Slr.mm, " mm")  # legend with slr
      labd <- paste0("Projection (SLR = ",
                     object$Projected.SLR$Slr.mm, " mm)")  # legend with slr
    }
    # -----------------------------------------------------------------
    # check if original time series contains missing values
    if (tp == 0) { # msl.trend object
      if (any(is.na(summ$MSL)) == TRUE) {
        p <- 0
      } else {
        p <- 1
      }
    }
    if (tp == 1) { # msl.forecast object
      if (any(is.na(summ$MSL[1:object$Historical.Record$Years])) == TRUE) {
        p <- 0
      } else {
        p <- 1
      }
    }
    # -----------------------------------------------------------------
    # specific settings, portions of object
    if (tp == 0) {
        n <- length(summ[, 1])
        n2 <- n - 3
        # dataframe without NA's for acceleration charts with msl.trend object
        summ2 <- summ[4:n2, ]
    } else {
        n <- length(summ[, 1])
        # dataframe without NA's for acceleration charts with msl.forecast object
        summ2 <- summ[4:n, ]
        # different plotting colours for historical
        summ3 <- summ[summ$Year <= object$Historical.Record$End, ]
        # different plotting colours for forecast
        summ4 <- summ[summ$Year >= object$Historical.Record$End, ]
    }
    # -----------------------------------------------------------------
    # Station name
    if (object$Station.Name == "Station Name not entered") {
        object$Station.Name <- NULL
    } else {
        # default setting is entered Station Name
        object$Station.Name <- object$Station.Name
    }
    # -----------------------------------------------------------------
    # If station name not entered
    if (file_name == " ") {
        file_name <- "File1.pdf"  # default
    } else {
        file_name <- file_name
    }
    # -----------------------------------------------------------------
    # If type not entered or entered outside range
    if (type == 1 | type == 2 | type == 3 | type == 4) {
        type <- type
    } else {
        print("default TYPE setting applied")
        type <- 1
    }
    if (ci == 1 | ci == 2) {
        # If ci not entered or entered outside range
        ci <- ci
    } else {
        print("default CONFIDENCE INTERVAL setting applied")
        ci <- 1
    }
    # -----------------------------------------------------------------
    # Confidence interval input
    if (ci == 2) {
        ci = 2.575  # multiplication factor for 99% CI
        lab1 <- paste("99% Confidence Interval")
    } else {
        # default setting is 95% CI
        ci = 1.96
        lab1 <- paste("95% Confidence Interval")
    }
    # -----------------------------------------------------------------
    # check slope of graph for locating chart 1 legends
    if (summ$Trend[n] - summ$Trend[1] < 0) {
        laba = paste("topright")  # Chart description to topright
        labb = paste("bottomleft")  # Chart key to bottomleft
    } else {
        # default setting
        laba = paste("topleft")  # Chart description to topleft
        labb = paste("bottomright")  # Chart key to bottomright
    }
    # -----------------------------------------------------------------
    if (requireNamespace("plyr", quietly = TRUE)) {
      plyr::round_any
    }
    if (n < 100) {
        # setting up x-axis plotting parameters
        xtic = 10  # year ticks on x-axis
        xlo <- plyr::round_any(min(summ[, 1]), 10, floor)
        xhi <- plyr::round_any(max(summ[, 1]), 10, ceiling)
    } else {
        # default
        xtic = 20
        xlo <- plyr::round_any(min(summ[, 1]), 20, floor)
        xhi <- plyr::round_any(max(summ[, 1]), 20, ceiling)
    }
    xlim = c(xlo, xhi)
    # -----------------------------------------------------------------
    # plotting routines for msl.trend objects
    # -----------------------------------------------------------------
    # plot 3 charts - time series, velocity and acceleration
    if ((type == 1) && (tp == 0)) {
        grDevices::pdf(file = file_name, width = 16.54, height = 20, pointsize = 32)
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mfrow = c(3, 1), las = 1)
        graphics::par(mar = c(0, 5.1, 2, 0.6), las = 1)  # chart 1 (time series)
        # -----------------------------------------------------------------
        # set y-axis scale
        if (p == 0) {
          # gap filling routine required
          M1 <- max(max(summ[, 2], na.rm = TRUE), max(summ[, 10]))
          M2 <- min(min(summ[, 2], na.rm = TRUE), min(summ[, 10]))
        }
        if (p == 1) {
          # no gaps in record
          M1 <- max(summ[, 2])
          M2 <- min(summ[, 2])
        }
        ylen <- M1 - M2
        if (ylen < 200) {
          # setting up y-axis plotting parameters
          ytic = 50  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 50, floor)
          yhi <- plyr::round_any(M1, 50, ceiling)
        } else {
          # default
          ytic = 100
          ylo <- plyr::round_any(M2, 100, floor)
          yhi <- plyr::round_any(M1, 100, ceiling)
        }
        ylim = c(ylo, yhi)
        # -----------------------------------------------------------------
        graphics::plot(summ[, 1], summ[, 2], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 2)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 3] + ci * summ[, 4]),
                                                rev((summ[, 3] - ci * summ[, 4]))),
                col = "azure3", border = NA)
        graphics::title(ylab = "Millimetres", font.lab = 2, cex.lab = 1.2, line = 3.5)
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        if (p == 0) {
            # gap filling routine required
            graphics::lines(summ[, 1], summ[, 10], col = "red")
            graphics::lines(summ[, 1], summ[, 2], lwd = 2)
            graphics::lines(summ[, 1], summ[, 3], lwd = 4)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                   inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "Gap Filling",
                                        "MSL Trend", lab1, "Peak Rate"),
                   text.font = c(2, 1, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 1, 3),
                   lwd = c(1, 2, 2, 4, 10, 4),
                   col = c("black", "black", "red", "black", "azure3", "blue"),
                   cex = c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9))
        }
        if (p == 1) {
            # no gap filling routine required
            graphics::lines(summ[, 1], summ[, 2], lwd = 2)
            graphics::lines(summ[, 1], summ[, 3], lwd = 4)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                   inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "MSL Trend",
                                        lab1, "Peak Rate"),
                   text.font = c(2, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 3),
                   lwd = c(1, 2, 4, 10, 4), col = c("black", "black", "black",
                                                    "azure3", "blue"),
                   cex = c(0.9, 0.9, 0.9, 0.9, 0.9))
        }
        graphics::par(mar = c(0, 5.1, 0, 0.6), las = 1)  # chart 2 (velocity)
        ylen <- max(summ[, 5] + ci * summ[, 6]) - min(summ[, 5] - ci * summ[, 6])
        if (ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 1 & ylen <= 2) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        if (ylen > 2 & ylen <= 10) {
            # setting up y-axis plotting parameters
            ytic = 1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   1, ceiling)
        }
        if (ylen > 10 & ylen <= 20) {
            # setting up y-axis plotting parameters
            ytic = 2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   2, ceiling)
        }
        if (ylen > 20) {
            # setting up y-axis plotting parameters
            ytic = 5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   5, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   5, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ[, 1], summ[, 5], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 5] + ci * summ[, 6]),
                                                rev((summ[, 5] - ci * summ[, 6]))),
                col = "azure3", border = NA)
        graphics::lines(summ[, 1], summ[, 5], lwd = 4)
        graphics::abline(h = 0, lty = 2)
        graphics::abline(h = max(summ[, 5]), lty = 3, col = "blue", lwd = 4)
        graphics::legend("topleft", legend = "Instantaneous Velocity",
               inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
        graphics::title(ylab = "Millimetres/year", font.lab = 2, cex.lab = 1.2,
                        line = 3.5)
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(mar = c(3.6, 5.1, 0, 0.6), las = 1)  # chart 3 (acceleration)
        ylen <- max(summ2[, 7] + ci * summ2[, 8]) - min(summ2[, 7] - ci * summ2[, 8])
        if (ylen <= 0.1) {
            # setting up y-axis plotting parameters
            ytic = 0.02  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.02, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.02, ceiling)
        }
        if (ylen > 0.1 & ylen <= 0.2) {
            # setting up y-axis plotting parameters
            ytic = 0.05  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.05, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.05, ceiling)
        }
        if (ylen > 0.2 & ylen <= 0.5) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 0.5 & ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        if (ylen > 1) {
            # setting up y-axis plotting parameters
            ytic = 0.5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.5, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.5, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ2[, 1], summ2[, 7], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ2[, 1], rev(summ2[, 1])),
                          c((summ2[, 7] + ci * summ2[, 8]),
                            rev((summ2[, 7] - ci * summ2[, 8]))), col = "azure3",
                          border = NA)
        graphics::lines(summ2[, 1], summ2[, 7], lwd = 4)
        graphics::abline(h = 0, lty = 2)
        graphics::abline(h = max(summ2[, 7]), lty = 3, col = "blue", lwd = 4)
        graphics::legend("topleft", legend = "Instantaneous Acceleration",
               inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
        graphics::title(ylab = "Millimetres/year/year", font.lab = 2,
                        cex.lab = 1.2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
        grDevices::dev.off()  # turn off pdf device
    }
    # -----------------------------------------------------------------
    # plot time series only
    if ((type == 2) && (tp == 0)) {
        grDevices::pdf(file = file_name, width = 16.54, height = 15, pointsize = 24)
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2.1, 0.6), las = 1)
        # -----------------------------------------------------------------
        # set y-axis scale
        if (p == 0) {
          # gap filling routine required
          M1 <- max(max(summ[, 2], na.rm = TRUE), max(summ[, 10]))
          M2 <- min(min(summ[, 2], na.rm = TRUE), min(summ[, 10]))
        }
        if (p == 1) {
          # no gaps in record
          M1 <- max(summ[, 2])
          M2 <- min(summ[, 2])
        }
        ylen <- M1 - M2
        if (ylen < 200) {
          # setting up y-axis plotting parameters
          ytic = 20  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 20, floor)
          yhi <- plyr::round_any(M1, 20, ceiling)
        } else {
          # default
          ytic = 50
          ylo <- plyr::round_any(M2, 50, floor)
          yhi <- plyr::round_any(M1, 50, ceiling)
        }
        ylim = c(ylo, yhi)
        # -----------------------------------------------------------------
        graphics::plot(summ[, 1], summ[, 2], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 1.7)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 3] + ci * summ[, 4]),
                                                rev((summ[, 3] - ci * summ[, 4]))),
                col = "azure3", border = NA)
        graphics::title(ylab = "Millimetres", font.lab = 2, cex.lab = 1.2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        if (p == 0) {
            # gap filling routine required
            graphics::lines(summ[, 1], summ[, 10], col = "red")
            graphics::lines(summ[, 1], summ[, 2], lwd = 2)
            graphics::lines(summ[, 1], summ[, 3], lwd = 4)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                   inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "Gap Filling",
                                        "MSL Trend", lab1),
                   text.font = c(2, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 1),
                   lwd = c(1, 2, 2, 4, 10),
                   col = c("black", "black", "red", "black", "azure3"),
                   cex = c(0.9, 0.9, 0.9, 0.9, 0.9))
        }
        if (p == 1) {
            # no gap filling routine required
            graphics::lines(summ[, 1], summ[, 2], lwd = 2)
            graphics::lines(summ[, 1], summ[, 3], lwd = 4)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                   inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "MSL Trend",
                                        lab1),
                   text.font = c(2, 1, 1, 1), lty = c(0, 1, 1, 1),
                   lwd = c(1, 2, 4, 10),
                   col = c("black", "black", "black", "azure3"),
                   cex = c(0.9, 0.9, 0.9, 0.9))
        }
        graphics::par(opar)  # restore original settings
        grDevices::dev.off()  # turn off pdf device
    }
    # -----------------------------------------------------------------
    # plot velocity series only
    if ((type == 3) && (tp == 0)) {
        grDevices::pdf(file = file_name, width = 16.54, height = 15, pointsize = 24)
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2.1, 0.6), las = 1)
        ylen <- max(summ[, 5] + ci * summ[, 6]) - min(summ[, 5] - ci * summ[, 6])
        if (ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 1 & ylen <= 2) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        if (ylen > 2 & ylen <= 5) {
            # setting up y-axis plotting parameters
            ytic = 0.5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.5, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.5, ceiling)
        }
        if (ylen > 5 & ylen <= 10) {
            # setting up y-axis plotting parameters
            ytic = 1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   1, ceiling)
        }
        if (ylen > 10) {
            # setting up y-axis plotting parameters
            ytic = 2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   2, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ[, 1], summ[, 5], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 1.7)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 5] + ci * summ[, 6]),
                                                rev((summ[, 5] - ci * summ[, 6]))),
                col = "azure3", border = NA)
        graphics::lines(summ[, 1], summ[, 5], lwd = 4)
        graphics::abline(h = 0, lty = 2)
        graphics::abline(h = max(summ[, 5]), lty = 3, col = "blue", lwd = 4)
        graphics::legend("topleft", legend = "Instantaneous Velocity",
               inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
        graphics::legend("bottomright", bg = "white", legend = c("KEY", "MSL Velocity",
                                                       "Peak Velocity", lab1),
               text.font = c(2, 1, 1, 1), lty = c(0, 1, 3, 1), lwd = c(1, 4, 4, 10),
               col = c("black", "black", "blue", "azure3"),
               cex = c(0.9, 0.9, 0.9, 0.9))
        graphics::title(ylab = "Millimetres/year", font.lab = 2, cex.lab = 1.2,
                        line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
        grDevices::dev.off()  # turn off pdf device
    }
    # -----------------------------------------------------------------
    # plot acceleration series only
    if ((type == 4) && (tp == 0)) {
        grDevices::pdf(file = file_name, width = 16.54, height = 15, pointsize = 24)
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2.1, 0.6), las = 1)
        ylen <- max(summ2[, 7] + ci * summ2[, 8]) - min(summ2[, 7] - ci * summ2[, 8])
        if (ylen <= 0.1) {
            # setting up y-axis plotting parameters
            ytic = 0.01  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.01, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.01, ceiling)
        }
        if (ylen > 0.1 & ylen <= 0.2) {
            # setting up y-axis plotting parameters
            ytic = 0.02  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.02, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.02, ceiling)
        }
        if (ylen > 0.2 & ylen <= 0.5) {
            # setting up y-axis plotting parameters
            ytic = 0.05  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.05, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.05, ceiling)
        }
        if (ylen > 0.5 & ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 1) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ2[, 1], summ2[, 7], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 1.7)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ2[, 1], rev(summ2[, 1])),
                          c((summ2[, 7] + ci * summ2[, 8]),
                            rev((summ2[, 7] - ci * summ2[, 8]))), col = "azure3",
                          border = NA)
        graphics::lines(summ2[, 1], summ2[, 7], lwd = 4)
        graphics::abline(h = 0, lty = 2)
        graphics::abline(h = max(summ2[, 7]), lty = 3, col = "blue", lwd = 4)
        graphics::legend("topleft", legend = "Instantaneous Acceleration",
               inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
        graphics::legend("bottomright", bg = "white",
                         legend = c("KEY", "MSL Acceleration", "Peak Acceleration",
                                    lab1),
               text.font = c(2, 1, 1, 1), lty = c(0, 1, 3, 1), lwd = c(1, 4, 4, 10),
               col = c("black", "black", "blue", "azure3"),
               cex = c(0.9, 0.9, 0.9, 0.9))
        graphics::title(ylab = "Millimetres/year/year", font.lab = 2,
                        cex.lab = 1.2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
        grDevices::dev.off()  # turn off pdf device
    }
    # -----------------------------------------------------------------
    # plotting routines for msl.forecast objects
    # -----------------------------------------------------------------
    # plot 3 charts - time series, velocity and acceleration
    if ((type == 1) && (tp == 1)) {
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        grDevices::pdf(file = file_name, width = 16.54, height = 20, pointsize = 32)
        graphics::par(mfrow = c(3, 1), las = 1)
        graphics::par(mar = c(0, 5.1, 2, 0.6), las = 1)  # chart 1 (time series)
        # -----------------------------------------------------------------
        # set y-axis scale
        if (p == 0) {
          # gap filling routine required
          M1 <- max(max(summ$MSL, na.rm = TRUE),
                    max(summ$Trend + ci * summ$TrendSD),
                    max(summ$FilledTS, na.rm = TRUE))
          M2 <- min(min(summ$MSL, na.rm = TRUE),
                    min(summ$Trend - ci * summ$TrendSD),
                    min(summ$FilledTS, na.rm = TRUE))
        }
        if (p == 1) {
          # no gaps in record
          M1 <- max(max(summ$MSL, na.rm = TRUE), max(summ$Trend + ci * summ$TrendSD))
          M2 <- min(min(summ$MSL, na.rm = TRUE), min(summ$Trend - ci * summ$TrendSD))
        }
        ylen <- M1 - M2
        if (ylen <= 500) {
          # setting up y-axis plotting parameters
          ytic = 50  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 50, floor)
          yhi <- plyr::round_any(M1, 50, ceiling)
        }
        if (ylen > 500 & ylen <= 1000) {
          # setting up y-axis plotting parameters
          ytic = 100  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 100, floor)
          yhi <- plyr::round_any(M1, 100, ceiling)
        }
        if (ylen > 1000) {
          # setting up y-axis plotting parameters
          ytic = 200  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 200, floor)
          yhi <- plyr::round_any(M1, 200, ceiling)
        }
        ylim = c(ylo, yhi)
        # -----------------------------------------------------------------
        graphics::plot(summ$Year, summ$MSL, type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 2)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 3] + ci * summ[, 4]),
                                                rev((summ[, 3] - ci * summ[, 4]))),
                col = "azure3", border = NA)
        graphics::title(ylab = "Millimetres", font.lab = 2, cex.lab = 1.2, line = 3.5)
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        if (p == 0) {
            # gap filling routine required
            graphics::lines(summ[, 1], summ[, 9], col = "red")
            graphics::lines(summ[, 1], summ[, 2], lwd = 2)
            graphics::lines(summ3[, 1], summ3[, 3], lwd = 4)
            graphics::lines(summ4[, 1], summ4[, 3], lwd = 4, col = "black", lty = 3)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                   inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "Gap Filling",
                                        "MSL Trend", lab1, "Projection"),
                   text.font = c(2, 1, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 1, 3),
                   lwd = c(1, 2, 2, 4, 10, 4),
                   col = c("black", "black", "red", "black", "azure3", "black"),
                   cex = c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9))
        }
        if (p == 1) {
            # no gap filling routine required
            graphics::lines(summ[, 1], summ[, 2], lwd = 2)
            graphics::lines(summ3[, 1], summ3[, 3], lwd = 4)
            graphics::lines(summ4[, 1], summ4[, 3], lwd = 4, col = "black", lty = 3)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                   inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "MSL Trend",
                                        lab1, "Projection"),
                   text.font = c(2, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 3),
                   lwd = c(1, 2, 4, 10, 4),
                   col = c("black", "black", "black", "azure3", "blue"),
                   cex = c(0.9, 0.9, 0.9, 0.9, 0.9))
        }
        graphics::par(mar = c(0, 5.1, 0, 0.6), las = 1)  # chart 2 (velocity)
        ylen <- max(summ[, 5] + ci * summ[, 6]) - min(summ[, 5] - ci * summ[, 6])
        if (ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        if (ylen > 1 & ylen <= 2) {
            # setting up y-axis plotting parameters
            ytic = 0.5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.5, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.5, ceiling)
        }
        if (ylen > 2 & ylen <= 10) {
            # setting up y-axis plotting parameters
            ytic = 1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   1, ceiling)
        }
        if (ylen > 10 & ylen <= 20) {
            # setting up y-axis plotting parameters
            ytic = 2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   2, ceiling)
        }
        if (ylen > 20) {
            # setting up y-axis plotting parameters
            ytic = 5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   5, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   5, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ[, 1], summ[, 5], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 5] + ci * summ[, 6]),
                                                rev((summ[, 5] - ci * summ[, 6]))),
                col = "azure3", border = NA)
        graphics::lines(summ3[, 1], summ3[, 5], lwd = 4)
        graphics::lines(summ4[, 1], summ4[, 5], lwd = 4, col = "black", lty = 3)
        graphics::abline(h = 0, lty = 2)
        graphics::legend("topleft", legend = "Instantaneous Velocity",
               inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
        graphics::title(ylab = "Millimetres/year", font.lab = 2, cex.lab = 1.2,
                        line = 3.5)
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(mar = c(3.6, 5.1, 0, 0.6), las = 1)  # chart 3 (acceleration)
        ylen <- max(summ2[, 7] + ci * summ2[, 8]) - min(summ2[, 7] - ci * summ2[, 8])
        if (ylen <= 0.1) {
            # setting up y-axis plotting parameters
            ytic = 0.02  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.02, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.02, ceiling)
        }
        if (ylen > 0.1 & ylen <= 0.2) {
            # setting up y-axis plotting parameters
            ytic = 0.05  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.05, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.05, ceiling)
        }
        if (ylen > 0.2 & ylen <= 0.5) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 0.5 & ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        if (ylen > 1) {
            # setting up y-axis plotting parameters
            ytic = 0.5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.5, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.5, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ2[, 1], summ2[, 7], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ2[, 1], rev(summ2[, 1])),
                          c((summ2[, 7] + ci * summ2[, 8]),
                            rev((summ2[, 7] - ci * summ2[, 8]))), col = "azure3",
                          border = NA)
        graphics::lines(summ3[, 1], summ3[, 7], lwd = 4)
        graphics::lines(summ4[, 1], summ4[, 7], lwd = 4, col = "black", lty = 3)
        graphics::abline(h = 0, lty = 2)
        graphics::legend("topleft", legend = "Instantaneous Acceleration",
               inset = c(-0.02, -0.01), bty = "n", text.font = 2, cex = 1.2)
        graphics::legend("bottomleft", legend = labc, inset = c(-0.02, -0.01),
                         bty = "n", text.font = 2, cex = 1.2)
        graphics::title(ylab = "Millimetres/year/year", font.lab = 2,
                        cex.lab = 1.2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
        grDevices::dev.off()  # turn off pdf device
    }
    # -----------------------------------------------------------------
    # plot time series only
    if ((type == 2) && (tp == 1)) {
        grDevices::pdf(file = file_name, width = 16.54, height = 15, pointsize = 24)
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2.1, 0.6), las = 1)
        # -----------------------------------------------------------------
        # set y-axis scale
        if (p == 0) {
          # gap filling routine required
          M1 <- max(max(summ$MSL, na.rm = TRUE),
                    max(summ$Trend + ci * summ$TrendSD),
                    max(summ$FilledTS, na.rm = TRUE))
          M2 <- min(min(summ$MSL, na.rm = TRUE),
                    min(summ$Trend - ci * summ$TrendSD),
                    min(summ$FilledTS, na.rm = TRUE))
        }
        if (p == 1) {
          # no gaps in record
          M1 <- max(max(summ$MSL, na.rm = TRUE), max(summ$Trend + ci * summ$TrendSD))
          M2 <- min(min(summ$MSL, na.rm = TRUE), min(summ$Trend - ci * summ$TrendSD))
        }
        ylen <- M1 - M2
        if (ylen < 500) {
          # setting up y-axis plotting parameters
          ytic = 50  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 50, floor)
          yhi <- plyr::round_any(M1, 50, ceiling)
        }
        if (ylen > 500 & ylen <= 1000) {
          # setting up y-axis plotting parameters
          ytic = 100  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 100, floor)
          yhi <- plyr::round_any(M1, 100, ceiling)
        }
        if (ylen > 1000) {
          # setting up y-axis plotting parameters
          ytic = 200  # year ticks on y-axis
          ylo <- plyr::round_any(M2, 200, floor)
          yhi <- plyr::round_any(M1, 200, ceiling)
        }
        ylim = c(ylo, yhi)
        # -----------------------------------------------------------------
        graphics::plot(summ[, 1], summ[, 2], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 1.7)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])), c((summ[, 3] + ci * summ[, 4]),
                                                rev((summ[, 3] - ci * summ[, 4]))),
                col = "azure3", border = NA)
        graphics::title(ylab = "Millimetres", font.lab = 2, cex.lab = 1.2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        if (p == 0) {
            # gap filling routine required
            graphics::lines(summ[, 1], summ[, 9], col = "red")
            graphics::lines(summ[, 1], summ[, 2], lwd = 2)
            graphics::lines(summ3[, 1], summ3[, 3], lwd = 4)
            graphics::lines(summ4[, 1], summ4[, 3], lwd = 4, col = "black", lty = 3)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                   inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "Gap Filling",
                                        "MSL Trend", lab1, labd),
                   text.font = c(2, 1, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 1, 3),
                   lwd = c(1, 2, 2, 4, 10, 4),
                   col = c("black", "black", "red", "black", "azure3", "black"),
                   cex = c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9))
        }
        if (p == 1) {
            # no gap filling routine required
            graphics::lines(summ[, 1], summ[, 2], lwd = 2)
            graphics::lines(summ3[, 1], summ3[, 3], lwd = 4)
            graphics::lines(summ4[, 1], summ4[, 3], lwd = 4, col = "black", lty = 3)
            graphics::legend(laba, legend = "Relative Mean Sea Level",
                   inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
            graphics::legend(labb, bg = "white",
                             legend = c("KEY", "Annual Average Data", "MSL Trend",
                                        lab1, "Projection"),
                   text.font = c(2, 1, 1, 1, 1), lty = c(0, 1, 1, 1, 3),
                   lwd = c(1, 2, 4, 10, 4),
                   col = c("black", "black", "black", "azure3", "black"),
                   cex = c(0.9, 0.9, 0.9, 0.9, 0.9))
        }
        graphics::par(opar)  # restore original settings
        grDevices::dev.off()  # turn off pdf device
    }
    # -----------------------------------------------------------------
    # plot velocity series only
    if ((type == 3) && (tp == 1)) {
        grDevices::pdf(file = file_name, width = 16.54, height = 15, pointsize = 24)
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2.1, 0.6), las = 1)
        ylen <- max(summ[, 5] + ci * summ[, 6]) - min(summ[, 5] - ci * summ[, 6])
        if (ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 1 & ylen <= 2) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        if (ylen > 2 & ylen <= 5) {
            # setting up y-axis plotting parameters
            ytic = 0.5  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   0.5, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   0.5, ceiling)
        }
        if (ylen > 5 & ylen <= 10) {
            # setting up y-axis plotting parameters
            ytic = 1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   1, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   1, ceiling)
        }
        if (ylen > 10) {
            # setting up y-axis plotting parameters
            ytic = 2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ[, 5] - ci * summ[, 6]) - ylen * 0.15,
                                   2, floor)
            yhi <- plyr::round_any(max(summ[, 5] + ci * summ[, 6]) + ylen * 0.15,
                                   2, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ[, 1], summ[, 5], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 1.7)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ[, 1], rev(summ[, 1])),
                          c((summ[, 5] + ci * summ[, 6]),
                            rev((summ[, 5] - ci * summ[, 6]))), col = "azure3",
                          border = NA)
        graphics::lines(summ3[, 1], summ3[, 5], lwd = 4)
        graphics::lines(summ4[, 1], summ4[, 5], lwd = 4, col = "black", lty = 3)
        graphics::abline(h = 0, lty = 2)
        graphics::legend("topleft", legend = "Instantaneous Velocity",
               inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
        graphics::legend("bottomright", bg = "white", legend = c("KEY", "MSL Velocity",
                                                       lab1, labd),
               text.font = c(2, 1, 1, 1), lty = c(0, 1, 1, 3), lwd = c(1, 4, 10, 4),
               col = c("black", "black", "azure3", "black"),
               cex = c(0.9, 0.9, 0.9, 0.9))
        graphics::title(ylab = "Millimetres/year", font.lab = 2, cex.lab = 1.2,
                        line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
        grDevices::dev.off()  # turn off pdf device
    }
    # -----------------------------------------------------------------
    # plot acceleration series only
    if ((type == 4) && (tp == 1)) {
        grDevices::pdf(file = file_name, width = 16.54, height = 15, pointsize = 24)
        opar <- graphics::par(no.readonly = TRUE)  # capture current settings
        graphics::par(mar = c(3.6, 5.1, 2.1, 0.6), las = 1)
        ylen <- max(summ2[, 7] + ci * summ2[, 8]) - min(summ2[, 7] - ci * summ2[, 8])
        if (ylen <= 0.1) {
            # setting up y-axis plotting parameters
            ytic = 0.01  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.01, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.01, ceiling)
        }
        if (ylen > 0.1 & ylen <= 0.2) {
            # setting up y-axis plotting parameters
            ytic = 0.02  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.02, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.02, ceiling)
        }
        if (ylen > 0.2 & ylen <= 0.5) {
            # setting up y-axis plotting parameters
            ytic = 0.05  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.05, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.05, ceiling)
        }
        if (ylen > 0.5 & ylen <= 1) {
            # setting up y-axis plotting parameters
            ytic = 0.1  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.1, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.1, ceiling)
        }
        if (ylen > 1) {
            # setting up y-axis plotting parameters
            ytic = 0.2  # year ticks on x-axis
            ylo <- plyr::round_any(min(summ2[, 7] - ci * summ2[, 8]) - ylen * 0.15,
                                   0.2, floor)
            yhi <- plyr::round_any(max(summ2[, 7] + ci * summ2[, 8]) + ylen * 0.15,
                                   0.2, ceiling)
        }
        ylim = c(ylo, yhi)
        graphics::plot(summ2[, 1], summ2[, 7], type = "line", lty = 0, xaxt = "n",
             yaxt = "n", xlab = " ", ylab = " ", xlim = xlim, ylim = ylim,
             main = object$Station.Name, cex.main = 1.7)
        graphics::rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "lightcyan1")
        graphics::polygon(c(summ2[, 1], rev(summ2[, 1])),
                          c((summ2[, 7] + ci * summ2[, 8]),
                            rev((summ2[, 7] - ci * summ2[, 8]))),
                          col = "azure3", border = NA)
        graphics::lines(summ3[, 1], summ3[, 7], lwd = 4)
        graphics::lines(summ4[, 1], summ4[, 7], lwd = 4, col = "black", lty = 3)
        graphics::abline(h = 0, lty = 2)
        graphics::legend("topleft", legend = "Instantaneous Acceleration",
               inset = c(-0.04, -0.01), bty = "n", text.font = 2, cex = 1.3)
        graphics::legend("bottomright", bg = "white",
                         legend = c("KEY", "MSL Acceleration", lab1, labd),
               text.font = c(2, 1, 1, 1), lty = c(0, 1, 1, 3),
               lwd = c(1, 4, 10, 4),
               col = c("black", "black", "azure3", "black"),
               cex = c(0.9, 0.9, 0.9, 0.9))
        graphics::title(ylab = "Millimetres/year/year", font.lab = 2,
                        cex.lab = 1.2, line = 3.5)
        graphics::title(xlab = "Year", font.lab = 2, cex.lab = 1.2, line = 2.5)
        graphics::axis(side = 1, at = seq(xlo, xhi, by = xtic))
        graphics::axis(side = 2, at = seq(ylo, yhi, by = ytic))
        graphics::par(opar)  # restore original settings
        grDevices::dev.off()  # turn off pdf device
    }
}

#' Summary outputs of decomposed time series.
#'
#' @param object of class \dQuote{msl.trend} (see \code{\link{msl.trend}}) or
#' \dQuote{msl.forecast} (see \code{\link{msl.forecast}}).
#'
#' @details This routine provides a screen summary of the respective outputs
#' from a \code{\link{msl.trend}} or \code{\link{msl.forecast}} object. The
#' summary produced is identical to str( ) for an object of class
#' \dQuote{msl.trend} (see \code{\link{msl.trend}}) or \dQuote{msl.forecast}
#' (see \code{\link{msl.forecast}}).
#'
#' @seealso \code{\link{msl.trend}}, \code{\link{msl.forecast}},
#' \code{\link{Balt}}, \code{\link{s}}, \code{\link{t}}.
#'
#' @examples
#' # -------------------------------------------------------------------------
#' # Isolate trend from Baltimore record, filling gaps with spline interpolation,
#' # 500 iterations and adding 1000 mm of slr to 2100. Use raw 'Balt.csv' data file.
#' # Note: ordinarily user would call 'File.csv' direct from working directory
#' # using the following sample code:
#' # s <- msl.trend('Balt.csv', fillgaps = 3, iter = 500, 'BALTIMORE, USA')
#' # t <- msl.forecast(s, slr = 1000)
#' # -------------------------------------------------------------------------
#'
#' data(s) # msl.trend object from above-mentioned example
#' data(t) # msl.forecast object from above-mentioned example
#' summary(s) # summary for object of class 'msl.trend' object
#' summary(t) # summary for object of class 'msl.forecast' object
#'
#' @export
summary <- function(object) {
    # -----------------------------------------------------------------
    # summary for msl.trend and msl.forecast objects
    # -----------------------------------------------------------------
    # If not msl.trend or msl.forecast object
    if (class(object) == "msl.trend" | class(object) == "msl.forecast") {
        class(object) <- class(object)
    } else {
        stop("object is not an msl.trend or msl.forecast object: no summary
             available")
    }
    # -----------------------------------------------------------------
    print(utils::str(object))
}

#' @importFrom changepoint cpt.var cpts
NULL

#' @importFrom forecast auto.arima
NULL

#' @importFrom plyr round_any
NULL

#' @importFrom Rssa grouping.auto igapfill reconstruct ssa
NULL

#' @importFrom tseries na.remove
NULL

#' @importFrom zoo na.approx na.spline
NULL

#' @importFrom grDevices dev.off graphics.off pdf
NULL

#' @importFrom graphics abline axis legend lines par plot polygon rect title
NULL

#' @importFrom stats predict sd smooth.spline spec.pgram ts
NULL

#' @importFrom utils read.csv str
NULL
