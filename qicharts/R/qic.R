#' Quality improvement charts
#'
#' Run and control charts for quality improvement and control
#'
#' @export
#' @import graphics
#' @import grDevices
#' @import stats
#' @importFrom utils tail
#'
#' @param y Numeric vector of counts or measures to plot. Mandatory.
#' @param n Numeric vector of sample sizes. Mandatory for P and U charts.
#' @param x Subgrouping vector used for aggregating data and making x-labels.
#'   Mandatory for Xbar and S charts.
#' @param data Data frame containing variables.
#' @param chart Type of control chart. Possible types are: \itemize{ \item
#'   "run": run chart (default). \item "i": individuals chart. \item "mr":
#'   moving range chart. \item "xbar": sample average chart. \item "s": sample
#'   standard deviation chart. \item "t": time between events chart. \item "p":
#'   proportions chart. \item "c": counts chart. \item "u": rates chart. \item
#'   "g": cases between events chart. }
#' @param notes Character vector of notes to be added to individual. data
#'   points.
#' @param cl Value specifying the center line (if known). Must be of length one
#'   or same as number of subgroups (for variable center line).
#' @param agg.fun String specifying the aggregate function if there is more than
#'   one value per subgroup. Possible values are  'mean' and 'sum'. Only
#'   relevant if you want to aggregate count data with run charts or I charts.
#'   If \code{agg.fun} = 'sum', the \code{n} argument (if provided) will be
#'   ignored.
#' @param ylim Range of y axis limits.
#' @param target Value specifying a target line to plot.
#' @param direction Value indication direction of improvement, 0 (down) or 1
#'   (up).
#' @param freeze Number identifying the last data point to include in
#'   calculations of center and limits (ignored if \code{breaks} argument is
#'   given).
#' @param breaks Numeric vector of break points. Useful for splitting graph in
#'   two or more sections with separate center line and control limits.
#' @param exclude Numeric vector of data points to exclude from calculations of
#'   center and control lines.
#' @param negy Logical value, if TRUE, the y axis is allowed to be negative
#'   (only relevant for I and Xbar charts).
#' @param dots.only Logical value. If TRUE, data points are not connected by
#'   lines and runs analysis is not performed. Useful for comparison and funnel
#'   plots.
#' @param multiply Integer indicating a number to multiply y axis by, e.g. 100
#'   for percents rather than proportions.
#' @param x.format Date format of x axis labels. See ?strftime for date formats.
#' @param nint Number indicating (approximately) the desired number of tick
#'   marks on the x axis.
#' @param cex Number indicating the amount by which text and symbols should be
#'   magnified.
#' @param main Character string specifying the title of the plot.
#' @param xlab Character string specifying the x axis label.
#' @param ylab Character string specifying the y axis label.
#' @param sub Character string specifying a subtitle to be printed in the lower
#'   left corner of the plot.
#' @param decimals Integer indicating the number of decimals shown for center
#'   and limits on the plot. Default behaviour is smart rounding to at least two
#'   significant digits.
#' @param pre.text Character string labelling pre-freeze period
#' @param post.text Character string labelling post-freeze period
#' @param llabs Character vector with four elements specifying labels for lower
#'   control limit, centre line, upper control limit and target line
#'   respectively
#' @param runvals Logical value, if TRUE, prints statistics from runs analysis
#'   on plot.
#' @param linevals Logical value, if TRUE, prints values for center and control
#'   lines on plot.
#' @param plot.chart Logical value, if TRUE, prints plot.
#' @param print.out Logical value, if TRUE, prints return value
#' @param prime Logical value, if TRUE, control limits incorporate
#'   between-subgroup variation as proposed by Laney (2002). This is recommended
#'   for data involving very large sample sizes. Only relevant for P and U
#'   charts.
#' @param standardised Logical value, if TRUE, creates a standardised control
#'   chart, where points are plotted in standard deviation units along with a
#'   center line at zero and control limits at 3 and -3. Only relevant for P, U
#'   and Xbar charts.
#' @param ... Further arguments to plot function.
#'
#' @details If \code{chart} is not specified, \code{qic()} plots a \strong{run
#'   chart}. Non-random variation will be marked by a dashed, yellow center line
#'   (the median) if either the longest run of data points above or below the
#'   median is longer than predicted or if the graph crosses the median fewer
#'   times than predicted (see references for details).
#'
#'   Only the \code{y} argument giving the count or measure of interest is
#'   mandatory for a run chart. If a denominator argument, \code{n}, is given,
#'   \eqn{y/n} will be plotted. If a subgrouping argument, \code{x}, is given,
#'   \eqn{sum(y)/sum(n)}, within each subgroup will be plotted. This behaviour
#'   can be modified using the \code{agg.fun} argument.
#'
#'   With \strong{controlcharts}, data aggregation by subgroups is handled
#'   according to chart type. For P, U, and I charts, data are aggregated as
#'   described for the run chart. For the C chart, the sum of counts,
#'   \code{sum(y)}, within each subgroups will be plotted.
#'
#'   For Xbar and S charts, the subgrouping argument, \code{x}, is mandatory.
#'   However, the sample size argument, \code{n}, is irrelevant and will be
#'   ignored.
#'
#'   The subgrouping argument, \code{x}, is irrelevant for T and G charts, and,
#'   if given, an error will occur if any subgroup has more than one element.
#'
#'   If more than one \code{note} is present within any subgroup, the first
#'   \code{note} (alphabetically) is chosen.
#'
#'   If both \code{prime} and \code{standardised} are \code{TRUE}, points are
#'   plotted in units corresponding to Laney's modified "standard deviation",
#'   which incorporates the variation between subgroups.
#'
#' @return A list of of class qic containing values and parameters of the qic
#'   plot.
#'
#' @references Runs analysis: \itemize{ \item Jacob Anhoej, Anne Vingaard Olesen
#'   (2014). Run Charts Revisited: A Simulation Study of Run Chart Rules for
#'   Detection of Non-Random Variation in Health Care Processes. PLoS ONE 9(11):
#'   e113825. doi: 10.1371/journal.pone.0113825 . \item Jacob Anhoej (2015).
#'   Diagnostic Value of Run Chart Analysis: Using Likelihood Ratios to Compare
#'   Run Chart Rules on Simulated Data Series. PLoS ONE 10(3): e0121349. doi:
#'   10.1371/journal.pone.0121349 \item Mark F. Schilling (2012). The Surprising
#'   Predictability of Long Runs. Math. Mag. 85, 141-149. \item Zhenmin Chen
#'   (2010). A note on the runs test. Model Assisted Statistics and Applications
#'   5, 73-77. } Calculation of control limits: \itemize{ \item   Douglas C.
#'   Montgomery (2009). Introduction to Statistical Process Control, Sixth
#'   Edition, John Wiley & Sons. \item James C. Benneyan (2001). Number-Between
#'   g-Type Statistical Quality Control Charts for Monitoring Adverse Events.
#'   Health Care Management Science 4, 305-318. \item Lloyd P. Provost, Sandra
#'   K. Murray (2011). The Health Care Data Guide: Learning from Data for
#'   Improvement. San Francisco: John Wiley & Sons Inc. \item David B. Laney
#'   (2002). Improved control charts for attributes. Quality Engineering, 14(4),
#'   531-537.}
#'
#' @examples
#' set.seed(1)
#' # Run chart of 24 samples of a random continuous variable
#' # with an approximate mean = 12 and standard deviation = 3.
#' y <- rnorm(24, 12, 3)
#' qic(y)
#'
#' # Add subgroup vector (dates) and a target
#' x <- seq.Date(as.Date('2013-08-04'), by = 'week', length = 24)
#' qic(y, x = x, target = 16)
#'
#' # Individuals control chart
#' qic(y, x = x, chart = 'i')
#'
#' # Xbar control chart, sample size = 5
#' y <- rnorm(5 * 24)
#' x <- rep(x, 5)
#' qic(y, x = x, chart = 'xbar')
#'
#' # Create data frame with counts and sample sizes by week
#' d <- data.frame(week = seq.Date(as.Date('2013-08-04'),
#'                                 by = 'week',
#'                                 length = 36),
#'                 y = c(rbinom(24, 20, 0.5), rbinom(12, 20, 0.8)),
#'                 n = round(rnorm(36, 20, 2)))
#'
#' # Proportions control chart
#' qic(y, n, x = week, data = d[1:24,], chart = 'p')
#'
#' # Introduce change in process performance
#' qic(y, n, x = week, data = d, chart = 'p')
#'
#' # Freeze baseline to first 24 samples
#' qic(y, n, x = week, data = d, chart = 'p', freeze = 24)
#'
#' # Break control chart before and after change
#' qic(y, n, x = week, data = d, chart = 'p', breaks = 24)
#'
#' # Introduce extreme sample value and notes
#' d$a <- ''
#' d$a[30] <- 'Extreme value'
#' d$y[30] <- 1
#' qic(y, n, x = week, data = d, chart = 'p',
#'     breaks = 24,
#'     notes = a)
#'
#' # Exclude value from calculations
#' d$a[30] <- 'Value excluded from calculations'
#' qic(y, n, x = week, data = d, chart = 'p',
#'     breaks = 24,
#'     notes = a,
#'     exclude = 30)

qic <- function(y,
                n,
                x,
                data,
                chart        = c('run',
                                 'i',
                                 'mr',
                                 'xbar',
                                 's',
                                 't',
                                 'p',
                                 'c',
                                 'u',
                                 'g'),
                notes        = NULL,
                cl           = NULL,
                agg.fun      = c('mean', 'sum'),
                ylim         = NULL,
                target       = NULL,
                direction    = NULL,
                freeze       = NULL,
                breaks       = NULL,
                exclude      = NULL,
                negy         = TRUE,
                dots.only    = FALSE,
                multiply     = 1,
                prime        = FALSE,
                standardised = FALSE,
                x.format     = '%Y-%m-%d',
                nint         = 5,
                cex          = 0.8,
                main,
                xlab         = 'Subgroup',
                ylab         = 'Indicator',
                sub          = NULL,
                decimals     = NULL,
                pre.text     = 'Before data',
                post.text    = 'After data',
                llabs        = c('LCL', 'CL', 'UCL', 'TRG'),
                runvals      = FALSE,
                linevals     = TRUE,
                plot.chart   = TRUE,
                print.out    = FALSE,
                ...) {

  # Select chart type
  type <- match.arg(chart)
  fn <- paste0('qic.', type)

  # Select aggregate function
  agg.fun <- match.arg(agg.fun)

  if(agg.fun != 'mean' & !missing(n)) {
    warning('\"n\" argument is irrelevant and will be ignored when \"agg.fun\" argument is provided.')
  }

  # Prepare chart title
  no_title <- missing(main)

  if(no_title) {
    main <- paste(toupper(type), "Chart of", deparse(substitute(y)))
    if(!missing(n)) {
      main <- paste(main, '/', deparse(substitute(n)))
    }
    if(multiply != 1) {
      main <- paste(main, 'x', multiply)
    }
    if(prime == TRUE)
      main <- paste(paste0(toupper(type), "'"), "Chart of", deparse(substitute(y)))

    if(standardised == TRUE)
      main <- paste("Standardised", main)
  }

  # if(no_title)
  #   main <- paste(toupper(type), "Chart of", deparse(substitute(y)))


  # Get data, sample sizes, subgroups, and notes
  if(!missing(data)){
    class(data) <- 'data.frame'
    y <- data[,deparse(substitute(y))]
    if(deparse(substitute(n)) %in% colnames(data))
      n <- data[,deparse(substitute(n))]
    if(deparse(substitute(x)) %in% colnames(data))
      x <- data[,deparse(substitute(x))]
    if(deparse(substitute(notes)) %in% colnames(data))
      notes <- data[,deparse(substitute(notes))]
    if(deparse(substitute(cl)) %in% colnames(data))
      cl <- data[,deparse(substitute(cl))]
  }

  # Check arguments
  if(missing(y))
    stop('\"y\" argument must be provided.')
  if(any(type == c('p', 'u')) & missing(n))
    stop('\"n\" argument must be provided for P and U charts.')
  if(any(type == c('xbar', 's')) & missing(x))
    stop('\"x\" argument must be provided for Xbar and S charts.')
  if(any(type == c('xbar', 's', 'c', 'g', 't')) & !missing(n)) {
    warning('\"n\" argument is only relevant for P, U, and I control charts and will be ignored.')
    n <- rep(1, length(y))
  }
  if(missing(n)) {
    n <- rep(1, length(y))
  } else {
    if(length(n) != length(y))
      stop('\"y\" and \"n\" arguments must have same length.')
  }
  if(missing(x)) {
    x <- 1:length(y)
  } else {
    if(is.factor(x))
      x <- droplevels(x)
    if(length(x) != length(y))
      stop('\"y\" and \"x\" arguments must have same length.')
  }
  if(length(target) > 1) {
    warning('\"target\" argument must be a single value. Argument ignored.')
    target <- NULL
  }
  if(!is.null(direction)) {
    if(!direction %in% c(0, 1)) {
      warning('\"direction\" must be 0 or 1. Argument ignored')
      direction <- NULL
    }
  }

  # Fix missing values
  cases <- complete.cases(y, n)
  y[!cases] <- NA
  n[!cases] <- NA

  # Clear freeze argument if breaks argument is given
  if(!is.null(breaks)){
    freeze <- NULL
    cl <- NULL
  }

  # Create data frame of values to analyse
  d <- data.frame(y.sum    = tapply(y, x, sum, na.rm = TRUE),
                  y.mean   = tapply(y, x, mean, na.rm = TRUE),
                  y.median = tapply(y, x, median, na.rm = TRUE),
                  y.sd     = tapply(y, x, sd, na.rm = TRUE),
                  y.n      = tapply(n, x, sum, na.rm = TRUE))

  # Check that subgroups are unique for T and G charts
  if(any(type == c('t', 'g')) & max(d$y.n, na.rm = TRUE) > 1)
    stop('The grouping argument, \"x\", must contain unique values for T and G charts.')

  # Replace NaN values with NA (if subgroup is empty)
  d[is.nan(d$y.mean),] <- NA

  # Get number of data points
  n.obs <- nrow(d)

  # Check for variable center line
  if(!length(cl) %in% c(0, 1, n.obs)){
    warning('\"cl\" argument should be either a single value or a vector of the same length as the number of subgroups. Argument ignored.')
    cl = NULL
  }

  if(length(cl) > 1 & type != 'run') {
    warning('Variable center line is only applicable to run charts.')
    cl = NULL
  }

  # Fix notes
  if(missing(notes)) {
    notes                  <- rep(NA, n.obs)
  } else {
    notes                  <- as.character(notes)
    notes[notes == '']     <- NA
    notes                  <- c(notes, rep(NA, length(y) - length(notes)))
    suppressWarnings(notes <- tapply(notes, x, min, na.rm = TRUE))
    dimnames(notes)        <- NULL
  }

  # Create indices of chart parts (if breaks argument is given)
  if(is.null(breaks))
    breaks <- n.obs
  breaks   <- breaks[order(breaks)]
  breaks   <- breaks[breaks > 1 & breaks < n.obs - 1]
  start    <- c(1, breaks + 1)
  end      <- c(breaks, n.obs)
  parts    <- list()

  for(i in 1:length(start)) {
    parts[[i]] <- seq(start[i], end[i])
  }

  # Perform data analysis according to chart type
  qic <- list()

  for(p in parts) {
    # Calculate indices of data points to exclude within each chart part
    ex <- exclude[exclude %in% p]
    if(length(ex)) {
      ex <- ex - min(p) + 1
    } else {
      ex <- NULL
    }
    # Build qic object
    y <- do.call(fn, list(d            = d[p,],
                          cl           = cl,
                          agg.fun      = agg.fun,
                          freeze       = freeze,
                          exclude      = ex,
                          prime        = prime,
                          standardised = standardised))
    qic$y   <- c(qic$y, y$y)
    qic$cl  <- c(qic$cl, y$cl)
    qic$lcl <- c(qic$lcl, y$lcl)
    qic$ucl <- c(qic$ucl, y$ucl)
  }

  # Suppress center and control lines in run charts if more than half
  #  the data points are of the same value
  if(max(table(qic$y)) >= (length(na.omit(qic$y)) / 2) & type == 'run') {
    qic$cl <- NA
    qic$ucl <- NA
    qic$lcl <- NA
  }

  # Prevent negative y axis if negy argument is FALSE
  if(!negy & min(qic$y, na.rm = TRUE) >= 0)
    qic$lcl[qic$lcl < 0] <- NA

  # Multiply y axis by multiply argument
  if(!standardised) {
    qic$y   <- as.vector(qic$y) * multiply
    qic$cl  <- as.vector(qic$cl) * multiply
    qic$lcl <- as.vector(qic$lcl) * multiply
    qic$ucl <- as.vector(qic$ucl) * multiply
    if(!is.null(target))
      target <- as.vector(target) * multiply
  }

  # Create x axis labels
  labels <- row.names(d)

  if(inherits(x, c('Date','POSIXct', 'POSIXt'))) {
    # labels <- as.Date(labels)
    labels <- as.POSIXct(labels)
    labels <- format(labels, format = x.format)
  }

  # Perform runs analysis
  if(is.null(ex)) {
    runs             <- sign(qic$y - qic$cl)
  } else {
    runs             <- sign(qic$y[-ex] - qic$cl[-ex])
  }

  runs               <- runs[runs != 0 & !is.na(runs)]
  n.useful           <- length(runs)

  if(n.useful) {
    run.lengths      <- rle(runs)$lengths
    n.runs           <- length(run.lengths)
    longest.run      <- max(run.lengths)
    longest.run.max  <- round(log2(n.useful)) + 3                # Schilling 2012
    n.crossings      <- max(n.runs - 1, 0)
    n.crossings.min  <- qbinom(0.05, max(n.useful - 1, 0), 0.5)  # Chen 2010 (7)
    runs.test        <- longest.run > longest.run.max |
      n.crossings < n.crossings.min
  } else {
    longest.run      <- NA
    longest.run.max  <- NA
    n.crossings      <- NA
    n.crossings.min  <- NA
    runs.test        <- FALSE
  }

  signals            <- which(qic$y > qic$ucl | qic$y < qic$lcl)

  # Complete qic object
  qic$agg.fun         <- agg.fun
  qic$target          <- target
  qic$direction       <- direction
  qic$n               <- as.vector(d$y.n)
  qic$labels          <- labels
  qic$notes           <- notes
  qic$parts           <- parts
  qic$freeze          <- freeze
  qic$exclude         <- exclude
  qic$n.obs           <- n.obs
  qic$n.useful        <- n.useful
  qic$longest.run     <- longest.run
  qic$longest.run.max <- longest.run.max
  qic$n.crossings     <- n.crossings
  qic$n.crossings.min <- n.crossings.min
  qic$runs.test       <- runs.test
  qic$signals         <- signals
  qic$dots.only       <- dots.only
  qic$decimals        <- decimals
  qic$runvals         <- runvals
  qic$linevals        <- linevals
  qic$main            <- main
  qic$xlab            <- xlab
  qic$ylab            <- ylab
  qic$sub             <- sub
  qic$pre.text        <- pre.text
  qic$post.text       <- post.text
  qic$llabs           <- llabs
  qic$ylim            <- ylim
  qic$nint            <- nint
  qic$cex             <- cex

  class(qic) <- 'qic'

  # Plot qic chart
  if(plot.chart)
    plot(qic, ...)

  # Return qic object
  if(print.out) {
    return(qic)
  } else {
    invisible(qic)
  }
}

qic.run <- function(d, freeze, cl, agg.fun, exclude, ...){
  # Calcutate indicator to plot

  switch(agg.fun,
         mean   = y <- d$y.sum / d$y.n,
         sum    = y <- d$y.sum)


  # Get number of subgroups
  y.length <- length(y)

  # Define subgroups to be used in calculations
  if(is.null(freeze))
    freeze <- y.length
  base <- 1:freeze
  if(!is.null(exclude))
    base <- base[-exclude]

  # Calculate center line
  if(is.null(cl))
    cl <- median(y[base], na.rm = TRUE)
  if(length(cl) == 1)
    cl <- rep(cl, y.length)

  # Calculate limits
  ucl <- rep(NA, y.length)
  lcl <- ucl

  # Return object to calling function
  return(list(y = y,
              cl = cl,
              lcl = lcl,
              ucl = ucl))
}

qic.i <- function(d, freeze, cl, agg.fun, exclude, ...) {

  # Get indicator to plot
  switch(agg.fun,
         mean   = y <- d$y.sum / d$y.n,
         sum    = y <- d$y.sum)

  # Get number of subgroups
  y.length <- length(y)

  # Define subgroups to be used in calculations
  if(is.null(freeze))
    freeze <- y.length
  base <- 1:freeze
  if(!is.null(exclude))
    base <- base[-exclude]

  # Calculate center line
  if(is.null(cl))
    cl <- mean(y[base], na.rm = TRUE)
  cl <- rep(cl, y.length)

  # Calculate moving ranges
  mr <- abs(diff(y[base]))
  amr <- mean(mr, na.rm = TRUE)

  # Calculate upper limit for moving ranges
  ulmr <- 3.267 * amr

  # Remove moving ranges bigger than ulmr and recalculate amr, Provost
  mr <- mr[mr < ulmr]
  amr <- mean(mr, na.rm = TRUE)

  # Calculate standard deviation, Montgomery, 6.33
  stdev <- amr/1.128
  stdev <- rep(stdev, y.length)

  # Calculate limits
  lcl <- cl - 3 * stdev
  ucl <- cl + 3 * stdev

  # Return object to calling function
  return(list(y = y,
              cl = cl,
              lcl = lcl,
              ucl = ucl))
}

qic.mr <- function(d, freeze, cl, agg.fun, exclude, ...) {

  # Calcutate indicator to plot
  switch(agg.fun,
         mean   = y <- d$y.sum / d$y.n,
         sum    = y <- d$y.sum)

  # y <- d$y.sum / d$y.n
  y <- c(NA, abs(diff(y)))

  # Get number of subgroups
  y.length <- length(y)

  # Define subgroups to be used in calculations
  if(is.null(freeze))
    freeze <- y.length
  base <- 1:freeze
  if(!is.null(exclude))
    base <- base[-exclude]

  # Calculate center line
  if(is.null(cl))
    cl <- mean(y[base], na.rm = TRUE)
  cl <- rep(cl, y.length)

  # Calculate upper limit for moving ranges
  ucl <- 3.27 * cl

  # Return object to calling function
  return(list(y = y,
              cl = cl,
              lcl = 0,
              ucl = ucl))
}

qic.t <- function(d, freeze, cl, exclude, ...) {

  # Get values to plot
  y <- d$y.mean

  if(min(y, na.rm = TRUE) < 0) {
    stop('Time between events cannot contain negative values')
  }

  if(min(y, na.rm = TRUE) == 0) {
    y[y == 0] <- 0.1
    warning('Time between events should not contain zero values. Zeros replaced by 0.1')
  }

  # Transform measures, Montgomery 7.28
  y <- y^(1 / 3.6)
  #   d$y.mean <- y
  d$y.sum <- y

  # Calculate center and limits for transformed values
  qic <- qic.i(d, freeze, cl, exclude, ...)

  # Back transform center line and limits
  y = qic$y^3.6
  cl = qic$cl^3.6
  ucl = qic$ucl^3.6
  lcl = qic$lcl^3.6
  lcl[lcl < 0 | is.nan(lcl)] <- NA

  return(list(y = y,
              cl = cl,
              lcl = lcl,
              ucl = ucl))
}

qic.xbar <- function(d, freeze, cl, exclude, standardised, ...){

  # Get values to plot
  y <- d$y.mean
  n <- d$y.n
  s <- d$y.sd

  # Get number of subgroups
  y.length <- length(y)

  # Define subgroups to be used in calculations
  if(is.null(freeze))
    freeze <- y.length
  base <- 1:freeze
  if(!is.null(exclude))
    base <- base[-exclude]

  # Calculate center line, Montgomery 6.30
  if(is.null(cl))
    cl <- sum(n[base] * y[base], na.rm = TRUE) / sum(n[base], na.rm = TRUE)
  cl <- rep(cl, y.length)

  # Calculate standard deviation and control limits, Montgomery 6.31
  stdev <- sqrt(sum(s[base]^2 * (n[base] - 1), na.rm = TRUE) /
                  sum(n[base] - 1, na.rm = TRUE))

  A3 <- a3(n)
  ucl <- cl + A3 * stdev
  lcl <- cl - A3 * stdev

  # Calculations for standardised control chart
  if(standardised) {
    y <- (y - cl) / (stdev / sqrt(n))
    cl <- rep(0, y.length)
    ucl <- rep(3, y.length)
    lcl <- rep(-3, y.length)
  }

  # Return object to calling function
  return(list(y = y,
              cl = cl,
              lcl = lcl,
              ucl = ucl))
}

qic.s <- function(d, freeze = NULL, cl, exclude, ...){

  # Get values to plot
  s <- d$y.sd
  n <- d$y.n
  excl <- which(is.na(s))

  # Get number of subgroups
  y.length <- length(s)

  # Define subgroups to be used in calculations
  if(is.null(freeze))
    freeze <- y.length
  base <- 1:freeze
  if(!is.null(exclude))
    base <- base[-exclude]

  # Calculate center line, Montgomery 6.31
  sbar <- sqrt(sum(s[base]^2 * (n[base] - 1), na.rm = TRUE) /
                 (sum(n[base], na.rm = TRUE) - y.length))
  cl <- rep(sbar, y.length)
  B3 <- b3(n)
  B4 <- b4(n)
  ucl <- B4 * sbar
  lcl <- B3 * sbar

  # Return object to calling function
  return(list(y = s,
              cl = cl,
              lcl = lcl,
              ucl = ucl))
}

qic.p <- function(d, freeze, cl, exclude, prime, standardised, ...){

  # Calcutate indicator to plot
  n <- d$y.sum
  N <- d$y.n
  y <- n / N

  # Get number of subgroups
  y.length <- length(y)

  # Define subgroups to be used in calculations
  if(is.null(freeze))
    freeze <- y.length
  base <- 1:freeze
  if(!is.null(exclude))
    base <- base[-exclude]

  # Calculate center line, Montgomery 7.7
  if(is.null(cl))
    cl <- sum(n[base], na.rm = TRUE) / sum(N[base], na.rm = TRUE)
  cl <- rep(cl, y.length)

  # Calculate standard deviation, Montgomery 7.8
  stdev <- sqrt(cl * (1 - cl) / N)

  # Calculate standard deviation for Laney's p-prime chart, incorporating
  # between-subgroup variation.
  if(prime) {
    z_i <- (y[base] - cl[base]) / stdev[base]
    sigma_z <- mean(abs(diff(z_i)), na.rm = TRUE) / 1.128
    stdev <- stdev * sigma_z
  }

  # Calculate limits
  ucl <- cl + 3 * stdev
  ucl[ucl > 1] <- NA
  lcl <- cl - 3 * stdev
  lcl[lcl < 0] <- NA

  # Calculations for standardised control chart, Montgomery 7.14
  if(standardised) {
    y <- (y - cl) / stdev # "z_i" in Montgomery
    cl <- rep(0, y.length)
    ucl <- rep(3, y.length)
    lcl <- rep(-3, y.length)
  }

  # Return object to calling function
  return(list(y = y,
              cl = cl,
              lcl = lcl,
              ucl = ucl))
}

qic.c <- function(d, freeze, cl, exclude, ...){

  # Calcutate indicator to plot
  y <- d$y.sum

  # Get number of subgroups
  y.length <- length(y)

  # Define subgroups to be used in calculations
  if(is.null(freeze))
    freeze <- y.length
  base <- 1:freeze
  if(!is.null(exclude))
    base <- base[-exclude]

  # Calculate center line
  if(is.null(cl))
    cl <- mean(y[base], na.rm = TRUE)
  cl <- rep(cl, y.length)

  # Calculate standard deviation, Montgomery 7.17
  stdev <- sqrt(cl)

  # Calculate limits
  ucl <- cl + 3 * stdev
  lcl <- cl - 3 * stdev
  lcl[lcl < 0] <- NA

  # Return object to calling function
  return(list(y = y,
              cl = cl,
              lcl = lcl,
              ucl = ucl))
}

qic.u <- function(d, freeze, cl, exclude, prime, standardised, ...){

  # Calcutate indicator to plot
  n <- d$y.sum
  N <- d$y.n

  y <- n / N

  # Get number of subgroups
  y.length <- length(y)

  # Define subgroups to be used in calculations
  if(is.null(freeze))
    freeze <- y.length
  base <- 1:freeze
  if(!is.null(exclude))
    base <- base[-exclude]


  # Calculate center line
  if(is.null(cl))
    cl <- sum(n[base], na.rm = TRUE) / sum(N[base], na.rm = TRUE)
  cl <- rep(cl, y.length)

  # Calculate standard deviation, Montgomery 7.19
  stdev <- sqrt(cl / N)

  # Calculate standard deviation for Laney's p-prime chart, incorporating
  # between-subgroup variation.
  if(prime) {
    z_i <- (y[base] - cl[base]) / stdev[base]
    sigma_z <- mean(abs(diff(z_i)), na.rm = TRUE) / 1.128
    stdev <- stdev * sigma_z
  }

  # Calculate limits
  ucl <- cl + 3 * stdev
  lcl <- cl - 3 * stdev
  lcl[lcl < 0] <- NA

  # Calculations for standardised control chart, Montgomery 7.20
  if(standardised) {
    y <- (y - cl) / stdev  # "u_i" in Montgomery
    cl <- rep(0, y.length)
    ucl <- rep(3, y.length)
    lcl <- rep(-3, y.length)
  }

  # Return object to calling function
  return(list(y = y,
              cl = cl,
              lcl = lcl,
              ucl = ucl))
}

qic.g <- function(d, freeze, cl, exclude, ...){

  # Calcutate indicator to plot
  y <- d$y.sum

  # Get number of subgroups
  y.length <- length(y)

  # Define subgroups to be used in calculations
  if(is.null(freeze))
    freeze <- y.length
  base <- 1:freeze
  if(!is.null(exclude))
    base <- base[-exclude]

  # Calculate center line
  if(is.null(cl))
    cl <- mean(y[base], na.rm = TRUE)
  cl <- rep(cl, y.length)

  # Calculate standard deviation:
  # Benneyan (2001)
  # Montgomery, p. 319
  stdev <- sqrt(cl * (cl + 1))

  # Calculate limits
  ucl <- cl + 3 * stdev
  lcl <- cl - 3 * stdev
  lcl[lcl < 0] <- NA

  #   Set center line to theoretical median, Provost (2011) p. 228
  cl <- 0.693 * cl

  # Return object to calling function
  return(list(y = y,
              cl = cl,
              lcl = lcl,
              ucl = ucl))
}

a3 <- function(n) {
  n[n == 0] <- NA
  tbl <- c(NA,
           2.659, 1.954, 1.628, 1.427, 1.287, 1.182,
           1.099, 1.032, 0.975, 0.927, 0.886, 0.850,
           0.817, 0.789, 0.763, 0.739, 0.718, 0.698,
           0.680, 0.663, 0.647, 0.633, 0.619, 0.606)
  x <- 3 / (4 * (n - 1)) * (4 * n - 3) / sqrt(n)
  w <- which(n <= 25)
  x[w] <- tbl[n[w]]
  x[is.nan(x)] <- NA
  return(x)
}

b3 <- function(n) {
  n[n == 0] <- NA
  tbl <- c(NA,
           0.000, 0.000, 0.000, 0.000, 0.030, 0.118,
           0.185, 0.239, 0.284, 0.321, 0.354, 0.382,
           0.406, 0.428, 0.448, 0.466, 0.482, 0.497,
           0.510, 0.523, 0.534, 0.545, 0.555, 0.565)
  x <- 1 - (3 / c4(n) / sqrt(2 * (n - 1)))
  w <- which(n <= 25)
  x[w] <- tbl[n[w]]
  x[is.nan(x)] <- NA
  return(x)
}

b4 <- function(n) {
  n[n == 0] <- NA
  tbl <- c(NA,
           3.267, 2.568, 2.266, 2.089, 1.970, 1.882,
           1.815, 1.761, 1.716, 1.679, 1.646, 1.618,
           1.594, 1.572, 1.552, 1.534, 1.518, 1.503,
           1.490, 1.477, 1.466, 1.455, 1.445, 1.435)
  x <- 1 + (3 / c4(n) / sqrt(2 * (n - 1)))
  w <- which(n <= 25)
  x[w] <- tbl[n[w]]
  x[is.nan(x)] <- NA
  return(x)
}

c4 <- function(n) {
  n[n == 0] <- NA
  tbl <- c(NA,
           0.7979, 0.8862, 0.9213, 0.9400, 0.9515, 0.9594,
           0.9650, 0.9693, 0.9727, 0.9754, 0.9776, 0.9794,
           0.9810, 0.9823, 0.9835, 0.9845, 0.9854, 0.9862,
           0.9869, 0.9876, 0.9882, 0.9887, 0.9892, 0.9896)

  x <- 4 * (n - 1) / (4 * n - 3)
  w <- which(n <= 25)
  x[w] <- tbl[n[w]]
  x[is.nan(x)] <- NA
  return(x)
}

#' Plot qic object
#'
#' @export
#' @importFrom utils tail
#'
#' @param x List object returned from the qic() function.
#' @param y Ignored. Included for compatibility with generic plot function.
#' @param ... Further arguments to plot function.
#'
#' @return Creates a qic plot.
#'
#' @examples
#' y <- rnorm(24)
#' p <- qic(y, plot.chart = FALSE)
#' plot(p)
#'
plot.qic <- function(x, y = NULL, ...) {
  col1            <- rgb(093, 165, 218, maxColorValue = 255) # blue
  col2            <- rgb(140, 140, 140, maxColorValue = 255) # grey
  col3            <- rgb(005, 151, 072, maxColorValue = 255) # green
  # col4            <- rgb(255, 165, 000, maxColorValue = 255) # yellow
  col4    <- rgb(241, 088, 084, maxColorValue = 255) # red
  n.obs           <- x$n.obs
  y               <- x$y
  cl              <- x$cl
  lcl             <- x$lcl
  ucl             <- x$ucl
  target          <- x$target
  direction       <- x$direction
  signals         <- x$signals
  runs.test       <- x$runs.test
  freeze          <- x$freeze
  parts           <- x$parts
  exclude         <- x$exclude
  labels          <- x$labels
  main            <- x$main
  xlab            <- x$xlab
  ylab            <- x$ylab
  sub             <- x$sub
  notes           <- x$notes
  n.useful        <- x$n.useful
  longest.run     <- x$longest.run
  longest.run.max <- x$longest.run.max
  n.crossings     <- x$n.crossings
  n.crossings.min <- x$n.crossings.min
  dots.only       <- x$dots.only
  decimals        <- x$decimals
  runvals         <- x$runvals
  linevals        <- x$linevals
  ylim            <- x$ylim
  pre.text        <- x$pre.text
  post.text       <- x$post.text
  llabs           <- x$llabs
  nint            <- x$nint
  cex             <- x$cex
  x               <- 1:n.obs
  ylim            <- range(ylim, y, ucl, lcl, cl, target, na.rm = TRUE)
  lwd             <- cex
  cex             <- par('cex') * cex
  cex2            <- cex * 0.9
  type            <- ifelse(dots.only, 'p', 'o')
  pch             <- ifelse(dots.only, 19, 20)

  # Setup plot margins
  mar <- par('mar')
  mar <- mar + c(-0.5, 1, -1, 3)

  if(runvals & !dots.only)
    mar <- mar + c(1.5, 0, 0, 0)

  if(!is.null(sub))
    mar <- mar + c(1, 0, 0, 0)

  op <- par(mar = mar, cex = cex, lwd = lwd)

  # setup empty plot area
  plot(x    = x,
       y    = y,
       type = 'n',
       axes = FALSE,
       ylim = ylim,
       xlab = '',
       ylab = '',
       ...)

  # add axes and title to plot
  if(dots.only) {
    at <- x
  } else {
    at <- axisTicks(range(x), log = FALSE, nint = nint)
  }

  axis(1,
       at = at,
       labels = labels[at],
       tcl = -0.2,
       lwd = 0,
       lwd.ticks = lwd,
       col = col2,
       ...)
  axis(2,
       tcl = -0.2,
       lwd = 0,
       lwd.ticks = lwd,
       col = col2,
       las = 2,
       ...)
  box(bty = 'l',
      lwd = lwd,
      col = col2)
  title(main = main,
        adj = 0,
        line = 1.7,
        # cex.main = cex * 1.25,
        font.main = 1)
  title(xlab = xlab, line = 2.8)
  title(ylab = ylab, line = 3.75)

  if(!is.null(sub))
    title(sub = sub,
          adj = 0,
          cex.sub = cex,
          line = ifelse(runvals & !dots.only, 5.8, 4.3))

  # colour center line according to runs analysis
  lty <- 1
  col <- col2

  if(runs.test & !dots.only) {
    col <- col4
    lty <- 5
  }

  # add lines to plot
  for(p in parts) {
    lines(p, cl[p], lty = lty, col = col, lwd = lwd * 1.5)
    lines(p, ucl[p], lty = 1, col = col2)
    lines(p, lcl[p], lty = 1, col = col2)
    lines(p, y[p], type = type, col = col1, lwd = lwd * 4, pch = pch)
  }

  # add target line
  if(!is.null(target))
    lines(x, rep(target, n.obs), lty = 3, col = col2, lwd = lwd * 1.5)

  # annotate before and after data if freeze argument is given
  if(!is.null(freeze)) {
    abline(v = freeze + 0.5,
           col = col2,
           lty = 3)
    mtext(pre.text,
          at = freeze / 2,
          cex = cex2,
          line = 0.7)
    mtext(post.text,
          at = freeze + (n.obs - freeze) / 2,
          cex = cex2,
          line = 0.7)
  }

  # colour data points outside sigma limits
  points(signals, y[signals], col = col4, pch = pch, cex = cex * 1.5)

  # mark excluded data points
  points(exclude, y[exclude], bg = 0, col = col2, pch = 21, cex = cex * 1.5)

  # add values for center and limits to the plot
  if(linevals) {
    val <- tail(na.omit(cl), 1)
    if(length(val) && !is.na(val)) {
      # mtext(paste('CL =', sround(val, decimals)),
      mtext(paste(llabs[2], '=', sround(val, decimals)),
            side = 4,
            at = val,
            las = 1,
            cex = cex2)
    }
    val <- tail((lcl), 1)
    if(length(val) && !is.na(val)) {
      mtext(paste(llabs[1], '=', sround(val, decimals)),
            side = 4,
            at = val,
            las = 1,
            cex = cex2)
    }
    val <- tail((ucl), 1)
    if(length(val) && !is.na(val)) {
      mtext(paste(llabs[3], '=', sround(val, decimals)),
            side = 4,
            at = val,
            las = 1,
            cex = cex2)
    }
    padj <- NA
    if(length(target) & !is.na(cl[1])) {
      ll <- c(lcl = tail(na.omit(lcl), 1),
              cl = tail(na.omit(cl), 1),
              ucl = tail(na.omit(ucl), 1))
      ll <- get(names(which.min(abs(tail(target, 1) - ll))))
      ll <- tail(ll, 1)
      padj <- ifelse(ll > tail(target, 1), 1.5, -0.5)
    }

    if(length(target)) {

      mtext(paste(llabs[4], '=', target), #sround(target, decimals)),
            side = 4,
            at = target,
            padj = padj,
            las = 1,
            cex = cex2)
    }
  }

  # Print statistics from runs analysis to plot
  if(runvals & !dots.only) {
    mtext(paste0('Obs. (useful) = ', sum(!is.na(y)),
                 ' (', n.useful, ')'),
          cex = cex2,
          side = 1,
          line = 4.5,
          adj = 0)
    mtext(paste0('Longest run (max) = ', longest.run,
                 ' (', longest.run.max, ')'),
          cex = cex2,
          side = 1,
          line = 4.5,
          adj = 0.5,
          col = ifelse(longest.run > longest.run.max, col4, 1))
    mtext(paste0('Crossings (min) = ', n.crossings,
                 ' (', n.crossings.min, ')'),
          cex = cex2,
          side = 1,
          adj = 1,
          line = 4.5,
          col = ifelse(n.crossings < n.crossings.min, col4, 1))
  }

  # Add notes to plot
  i.ann <- which(!is.na(notes))
  if(length(i.ann)) {
    x.ann <- x[i.ann]
    y.ann <- y[i.ann]
    mtext(notes,
          side = 3,
          at = x,
          adj = 0.5,
          cex = cex2)
    segments(x.ann,
             max(ylim, na.rm = TRUE) * 1.05,
             y1 = y.ann,
             lty = 3)#, lwd = lwd)
  }

  # Add direction arrow
  if(!is.null(direction) & length(na.omit(cl))) {
    a <- ifelse(direction, expression('' %->% ''), expression('' %<-% ''))
    mtext(a, side = 4, line = -1.2, at = tail(na.omit(cl), 1), col = col2)
  }

  par(op)
}

# Smart rounding for labels to at least 2 significant digits
sround <- function(x, dec) {
  if(is.null(dec)) {
    n <- nchar(as.character(floor(x)))
    return(signif(x, max(2, n)))
  }
  return(round(x, dec))
}
