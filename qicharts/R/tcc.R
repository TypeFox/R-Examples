#' Trellis Control Charts
#'
#' Run and control charts for multivariate data i trellis (grid) layout.
#'
#' @export
#' @import ggplot2
#' @import ggrepel
#' @importFrom utils tail
#'
#' @param n Numerator, numeric vector of counts or measures to plot. Mandatory.
#' @param d Denominator, numeric vector of subgroup sizes. Mandatory for P and U
#'   charts.
#' @param x Subgrouping vector used for aggregating data by subgroup and making
#'   x-labels. Mandatory for Xbar and S charts.
#' @param g1 Grouping vector 1 used for trellis layout (facets).
#' @param g2 Grouping vector 2 used for trellis layout (facets).
#' @param breaks Numeric vector of break points. Useful for splitting graph in
#'   two or more sections with separate center line and control limits.
#' @param notes Character vector of notes to be added to individual. data
#'   points.
#' @param data Data frame containing variables.
#' @param chart Type of control chart. Possible types are: \itemize{ \item
#'   "run": run chart (default). \item "i": individuals chart. \item "mr":
#'   moving range chart. \item "xbar": sample average chart. \item "s": sample
#'   standard deviation chart. \item "t": time between events chart. \item "p":
#'   proportions chart. \item "c": counts chart. \item "u": rates chart. \item
#'   "g": cases between events chart. }
#' @param multiply Integer indicating a number to multiply y axis by, e.g. 100
#'   for percents rather than proportions. See also \code{y.percent} argument.
#' @param freeze Number identifying the last data point to include in
#'   calculations of center and limits (ignored if \code{breaks} argument is
#'   given).
#' @param exclude Numeric vector of data points to exclude from runs analysis
#'   and calculations of center and control lines (same for each facet).
#' @param target Numeric value indicating a target value to be plotted as a
#'   horizontal line (same for each facet).
#' @param n.sum Logical value indicating whether the mean (default) or sum of
#'   numerator (n argument) per subgroup should be plotted. Only relevant for
#'   run, C, and I charts with multiple counts per subgroup.
#' @param y.neg Logical value. If TRUE (default), the y axis is allowed to be
#'   negative (only relevant for I and Xbar charts).
#' @param y.percent Logical. If TRUE, formats y axis labels as percent.
#' @param y.expand Numeric value to include in y axis. Useful e.g. for beginning
#'   y axis at zero.
#' @param x.pad Number indicating expansion of x axis to make room for center
#'   line labels.
#' @param x.date.format Date format of x axis labels. See \code{?strftime()} for
#'   possible date formats.
#' @param cl.lab Logical value. If TRUE (default), plots center line labels.
#' @param cl.decimals Number of decimals on center line labels.
#' @param main Character string specifying the title of the plot.
#' @param xlab Character string specifying the x axis label.
#' @param ylab Character string specifying the y axis label.
#' @param cex Number indicating the amount by which text should be magnified.
#' @param pex Number indicating the amount by which plotting symbols should be
#'   magnified.
#' @param prime Logical value, If TRUE (default unless \code{dots.only = TRUE}),
#'   control limits incorporate between-subgroup variation as proposed by Laney
#'   (2002). Only relevant for P and U charts.
#' @param flip Logical. If TRUE rotates the plot 90 degrees.
#' @param dots.only Logical value. If TRUE, data points are not connected by
#'   lines, prime is forced to be FALSE, and runs analysis is not performed.
#'   Useful for comparison and funnel plots.
#' @param print.summary Logical. If TRUE, prints summary of tcc object.
#' @param ... Further arguments to ggplot function.
#'
#' @details \code{tcc()} is a wrapper function for \code{ggplot2()} that makes
#'   multivariate run and control charts. It takes up to two grouping variables
#'   for multidimensional trellis plots.
#'
#'   Note that, in contrast to the qic() function, the prime argument defaults
#'   to TRUE, which means that control limits of P and U charts by default
#'   incorporate between-subgroup variation as proposed by Laney (2002).
#'
#' @return An object of class ggplot.
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
#' # Run chart of 24 random normal variables
#' tcc(rnorm(24))
#'
#' # Build data frame for examples
#' d <- data.frame(x = rep(1:24, 4),
#'                 mo = (rep(seq(as.Date('2014-1-1'),
#'                               length.out = 24,
#'                               by = 'month'),
#'                           4)),
#'                 n = rbinom(4 * 24, 100, 0.5),
#'                 d = round(runif(4 * 24, 90, 110)),
#'                 g1 = rep(c('a', 'b'), each = 48),
#'                 g2 = rep(c('A', 'B'), each = 24))
#'
#' # Run chart with two grouping variables
#' tcc(n, d, mo, g1 = g1, g2 = g2, data = d)
#'
#' # P chart
#' tcc(n, d, mo, g1 = g1, g2 = g2, data = d, chart = 'p')
#'
#' # P chart with baseline fixed to the first 12 data points
#' tcc(n, d, mo, g1 = g1, g2 = g2, data = d, chart = 'p', freeze = 12)
#'
#' # P chart with two breaks and summary output
#' tcc(n, d, mo, g1 = g1, g2 = g2, data = d, chart = 'p',
#'  breaks = c(12, 18), print.summary = TRUE)

tcc <- function(n, d, x, g1, g2, breaks, notes,
                data,
                chart         = c("run", "i", "mr", "xbar", "s",
                                  "t", "p", "c", "u", "g"),
                multiply      = 1,
                freeze        = NULL,
                exclude,
                target        = NA,
                n.sum         = FALSE,
                y.neg         = TRUE,
                y.percent     = FALSE,
                y.expand      = NULL,
                x.pad         = 1,
                x.date.format = NULL,
                cl.lab        = TRUE,
                cl.decimals   = 2,
                main,
                xlab          = 'Subgroup',
                ylab          = 'Value',
                cex           = 1,
                pex           = 1,
                prime         = TRUE,
                flip          = FALSE,
                dots.only     = FALSE,
                print.summary = FALSE,
                ...) {
  # Get chart type
  type <- match.arg(chart)
  fn <- paste0('tcc.', type)

  # Build chart title
  if(missing(main)) {
    main <- paste(toupper(type), "Chart of", deparse(substitute(n)))
    if(!missing(d)) {
      main <- paste(main, '/', deparse(substitute(d)))
    }
    if(multiply != 1) {
      main <- paste(main, 'x', multiply)
    }
  }

  # Get data
  if(!missing(data)){
    class(data) <- 'data.frame'
    n <- data[,deparse(substitute(n))]
    if(deparse(substitute(d)) %in% colnames(data))
      d <- data[,deparse(substitute(d))]
    if(deparse(substitute(x)) %in% colnames(data))
      x <- data[,deparse(substitute(x))]
    if(deparse(substitute(g1)) %in% colnames(data))
      g1 <- data[,deparse(substitute(g1))]
    if(deparse(substitute(g2)) %in% colnames(data))
      g2 <- data[,deparse(substitute(g2))]
    if(deparse(substitute(notes)) %in% colnames(data))
      notes <- data[,deparse(substitute(notes))]
  }

  # Set missing denominator
  no.d <- missing(d)
  if(no.d)
    d <- rep(1, length(n))

  # Set missing subgroups
  if(missing(x))
    x <- seq_along(n)
  if(missing(g1))
    g1 <- rep(1, length(n))
  if(missing(g2))
    g2 <- rep(1, length(n))
  if(missing(notes))
    notes <- rep(NA, length(n))

  # Fix missing values
  cases <- complete.cases(n, d)
  n[!cases] <- NA
  d[!cases] <- NA
  cases <- complete.cases(n, d)

  # Prevent prime if dots.only
  if(dots.only) {
    prime = FALSE
  }

  # Initialise data frame
  df <- data.frame(n, d, x, g1, g2, cases, notes)
  df <- droplevels(df)

  # Build breaks variable
  if(missing(breaks)) {
    df$breaks <- rep(1, nrow(df))
  } else {
    freeze <- NULL
    df <- split(df, list(df$g1, df$g2))
    df <- lapply(df, function(x) {
      if(!all(breaks %in% 2:(nrow(x) - 2))) {
        warning('Invalid \"breaks\" argument')
        x$breaks <- rep(1, nrow(x))
      } else {
        breaks <- c(0, breaks, nrow(x))
        breaks <- sort(breaks)
        breaks <- diff(breaks)
        breaks <- rep(seq_along(breaks), breaks)
        x$breaks <- breaks
      }
      return(x)
    })
    df <- do.call(rbind, df)
    df <- df[order(df$g1, df$g2, df$breaks, df$x), ]
  }

  # Calculate values to plot
  d1 <- aggregate(cbind(n, d) ~ x + g1 + g2 + breaks,
                  data = df,
                  FUN = sum,
                  na.rm = TRUE,
                  na.action = na.pass)
  d2 <- aggregate(cbind(s = n) ~ x + g1 + g2 + breaks,
                  data = df,
                  FUN = sd,
                  na.rm = TRUE,
                  na.action = na.pass)
  d3 <- aggregate(cbind(n.obs = cases) ~ x + g1 + g2 + breaks,
                  data = df,
                  FUN = sum,
                  na.rm = TRUE,
                  na.action = na.pass)
  try(d4 <- aggregate(notes ~ x + g1 + g2 + breaks,
                      data = df,
                      FUN = paste,
                      collapse=' | '),
      silent = TRUE)
  df <- merge(d1, d2)
  df <- merge(df, d3)

  if(exists('d4')) {
    df <- merge(df, d4, all = T)
  } else{
    df$notes <- NA
  }

  # Check that subgroups are unique for T and G charts
  if(any(type == c('t', 'g')) & max(df$n.obs, na.rm = TRUE) > 1)
    stop('The grouping argument, \"x\", must contain unique values for T and G charts.')

  df <- df[order(df$x), ]

  # Calculate y variable
  if(n.sum & no.d | type %in% c('c', 't', 'g')) {
    df$y <- df$n
  } else {
    df$y <- df$n / df$d
  }

  df$y[df$d == 0] <- NA

  # Build exclude variable
  df$exclude <- FALSE
  if(!missing(exclude)) {
    df <- split(df, list(df$g1, df$g2))
    df <- lapply(df, function(x) {x$exclude[exclude] <- TRUE; return(x)})
    df <- do.call(rbind, df)
  }

  # Complete data frame
  df                         <- split(df, list(df$g1, df$g2, df$breaks))
  df                         <- lapply(df, fn,
                                       freeze = freeze,
                                       prime = prime,
                                       n.sum)
  df                         <- lapply(df, runs.analysis)
  df                         <- do.call(rbind, df)
  row.names(df)              <- NULL
  num.cols                   <- c('n', 'd', 's', 'n.obs',
                                  'cl', 'lcl', 'ucl', 'y')
  mult.cols                  <- c('cl', 'lcl', 'ucl', 'y')
  df[num.cols]               <- sapply(df[num.cols], as.numeric)
  df[mult.cols]              <- sapply(df[mult.cols], function(x) x * multiply)
  df$limits.signal           <- df$y < df$lcl | df$y > df$ucl
  x                          <- is.na(df$limits.signal)
  df$limits.signal[x]        <- FALSE
  df$ucl[!is.finite(df$ucl)] <- NA
  df$lcl[!is.finite(df$lcl)] <- NA
  df$target                 <- as.numeric(target)

  # Prevent negative y axis if negy argument is FALSE
  if(!y.neg & min(df$y, na.rm = TRUE) >= 0)
    df$lcl[df$lcl < 0] <- NA

  # Build return value
  tcc <- list(df       = df,
              main     = main,
              xlab     = xlab,
              ylab     = ylab,
              n.sum    = n.sum,
              multiply = multiply,
              freeze   = freeze,
              y.neg    = y.neg,
              prime    = prime)

  # Plot and return
  p <- plot.tcc(tcc, cex = cex, pex = pex, cl.decimals = cl.decimals,
                x.pad = x.pad, cl.lab = cl.lab, y.expand = y.expand,
                y.percent = y.percent, x.date.format = x.date.format,
                flip = flip, dots.only = dots.only,
                ...)

  class(p) <- c('tcc', class(p))

  if(print.summary)
    print(summary(p))

  return(p)

  # if(print) {
  #   return(tcc)
  # } else {
  #   invisible(tcc)
  # }
}

tcc.run <- function(df, freeze, ...) {
  l <- nrow(df)
  if(is.null(freeze))
    freeze <- l
  base <- 1:freeze
  base <- base[!df$exclude]
  y    <- df$y
  cl   <- median(y[base], na.rm = TRUE)
  cl   <- rep(cl, l)
  lcl  <- NA
  ucl  <- NA
  df   <- cbind(df, cl, ucl, lcl)

  return(df)
}

tcc.i <- function(df, freeze, ...) {
  l <- nrow(df)
  if(is.null(freeze))
    freeze <- l
  base <- 1:freeze
  base <- base[!df$exclude]
  y    <- df$y

  # Calculate centre line
  cl   <- mean(y[base], na.rm = TRUE)
  cl   <- rep(cl, l)

  # Average moving range
  mr  <- abs(diff(y))
  amr <- mean(mr, na.rm = TRUE)

  # Upper limit for moving ranges
  ulmr <- 3.267 * amr

  # Remove moving ranges greater than ulmr and recalculate amr, Provost p.156
  mr  <- mr[mr < ulmr]
  amr <- mean(mr, na.rm = TRUE)

  # Calculate standard deviation, Montgomery, 6.33
  stdev <- amr / 1.128

  # Calculate limits
  lcl <- cl - 3 * stdev
  ucl <- cl + 3 * stdev

  # Build and return data frame
  df <- cbind(df, cl, ucl, lcl)
  return(df)
}

tcc.mr <- function(df, freeze, ...) {
  l <- nrow(df)
  if(is.null(freeze))
    freeze <- l
  base <- 1:freeze
  base <- base[!df$exclude]
  y    <- df$y
  y    <- c(NA, abs(diff(y)))
  df$y <- y

  # Calculate centre line
  cl <- mean(y[base], na.rm = TRUE)
  cl <- rep(cl, l)

  # Calculate upper limit for moving ranges
  lcl <- NA
  ucl <- 3.27 * cl

  # Build and return data frame
  df <- cbind(df, cl, ucl, lcl)
  return(df)
}

tcc.t <- function(df, freeze, ...) {
  if(min(df$y, na.rm = TRUE) < 0) {
    stop('Time between events cannot contain negative values')
  }

  if(min(df$y, na.rm = TRUE) == 0) {
    df$y[df$y == 0] <- 0.1
    warning('Time between events should not contain zero values. Zeros replaced by 0.1')
  }

  d   <- df
  d$y <- d$y^(1 / 3.6)

  d <- tcc.i(d, freeze, ...)

  # Back transform centre line and limits
  # y = d$y^3.6
  cl  <- d$cl^3.6
  ucl <- d$ucl^3.6
  lcl <- d$lcl^3.6
  lcl[lcl < 0 | is.nan(lcl)] <- NA

  # Build and return data frame
  df <- cbind(df, cl, ucl, lcl)
  return(df)
}

tcc.xbar <- function(df, freeze, ...){
  l <- nrow(df)
  if(is.null(freeze))
    freeze <- l
  base <- 1:freeze
  base <- base[!df$exclude]
  y    <- df$y
  n    <- df$n.obs
  s    <- df$s

  # Calculate centre line, Montgomery 6.30
  cl <- sum(n[base] * y[base], na.rm = TRUE) / sum(n[base], na.rm = TRUE)
  cl <- rep(cl, l)

  # Calculate standard deviation and control limits, Montgomery 6.31
  stdev <- sqrt(sum(s[base]^2 * (n[base] - 1), na.rm = TRUE) /
                  sum(n[base] - 1, na.rm = TRUE))
  A3    <- a3(n)
  ucl   <- cl + A3 * stdev
  lcl   <- cl - A3 * stdev

  # Build and return data frame
  df <- cbind(df, cl, ucl, lcl)
  return(df)
}

tcc.s <- function(df, freeze, ...){
  l <- nrow(df)
  if(is.null(freeze))
    freeze <- l
  base <- 1:freeze
  base <- base[!df$exclude]
  s    <- df$s
  df$y <- s
  n    <- df$n.obs

  # Calculate centre line, Montgomery 6.31
  sbar <- sqrt(sum(s[base]^2 * (n[base] - 1), na.rm = TRUE) /
                 (sum(n[base], na.rm = TRUE) - l))
  cl   <- rep(sbar, l)
  B3   <- b3(n)
  B4   <- b4(n)
  ucl  <- B4 * sbar
  lcl  <- B3 * sbar

  # Build and return data frame
  df <- cbind(df, cl, ucl, lcl)
  return(df)
}

tcc.p <- function(df, freeze, prime, ...) {
  l <- nrow(df)
  if(is.null(freeze))
    freeze <- l
  base  <- 1:freeze
  base  <- base[!df$exclude]
  n     <- df$n
  d     <- df$d
  y     <- df$y
  cl    <- sum(n[base], na.rm = TRUE) / sum(d[base], na.rm = TRUE)
  cl    <- rep(cl, l)
  stdev <- sqrt(cl * (1 - cl) / d)
  # Calculate standard deviation for Laney's p-prime chart, incorporating
  # between-subgroup variation.
  if(prime) {
    z_i     <- (y[base] - cl[base]) / stdev[base]
    sigma_z <- mean(abs(diff(z_i)), na.rm = TRUE) / 1.128
    stdev   <- stdev * sigma_z
  }

  ucl          <- cl + 3 * stdev
  lcl          <- cl - 3 * stdev
  ucl[ucl > 1] <- NA
  lcl[lcl < 0] <- NA

  # Build and return data frame
  df <- cbind(df, cl, ucl, lcl)
  return(df)
}

tcc.c <- function(df, freeze, ...){
  l        <- nrow(df)
  if(is.null(freeze))
    freeze <- l
  base <- 1:freeze
  base <- base[!df$exclude]
  y  <- df$y
  cl <- mean(y[base], na.rm = TRUE)
  cl <- rep(cl, l)

  # Calculate standard deviation, Montgomery 7.17
  stdev <- sqrt(cl)

  # Calculate limits
  ucl          <- cl + 3 * stdev
  lcl          <- cl - 3 * stdev
  lcl[lcl < 0] <- NA

  # Build and return data frame
  df <- cbind(df, cl, ucl, lcl)
  return(df)
}

tcc.u <- function(df, freeze, prime, ...){
  l        <- nrow(df)
  if(is.null(freeze))
    freeze <- l
  base <- 1:freeze
  base <- base[!df$exclude]
  n    <- df$n
  d    <- df$d
  y    <- n / d
  cl   <- sum(n[base], na.rm = TRUE) / sum(d[base], na.rm = TRUE)
  cl   <- rep(cl, l)

  # Calculate standard deviation, Montgomery 7.19
  stdev <- sqrt(cl / d)

  # Calculate standard deviation for Laney's u-prime chart, incorporating
  # between-subgroup variation.
  if(prime) {
    z_i     <- (y[base] - cl[base]) / stdev[base]
    sigma_z <- mean(abs(diff(z_i)), na.rm = TRUE) / 1.128
    stdev   <- stdev * sigma_z
  }

  # Calculate limits
  ucl          <- cl + 3 * stdev
  lcl          <- cl - 3 * stdev
  lcl[lcl < 0] <- NA

  # Build and return data frame
  df <- cbind(df, cl, ucl, lcl)
  return(df)
}

tcc.g <- function(df, freeze, ...){
  l <- nrow(df)
  if(is.null(freeze))
    freeze <- l
  base <- 1:freeze
  base <- base[!df$exclude]
  y    <- df$y

  # Calculate centre line
  cl <- mean(y[base], na.rm = TRUE)
  cl <- rep(cl, l)

  # Calculate standard deviation, Montgomery, p. 319
  stdev <- sqrt(cl * (cl + 1))

  # Calculate limits
  ucl          <- cl + 3 * stdev
  lcl          <- cl - 3 * stdev
  lcl[lcl < 0] <- NA

  # Set centre line to theoretical median, Provost (2011) p. 228
  cl <- 0.693 * cl

  # Build and return data frame
  df <- cbind(df, cl, ucl, lcl)
  return(df)
}

runs.analysis <- function(df) {
  y                  <- df$y[!df$exclude]
  cl                 <- df$cl[!df$exclude]
  runs               <- sign(y - cl)
  runs               <- runs[runs != 0 & !is.na(runs)]
  n.useful           <- length(runs)

  if(n.useful) {
    run.lengths      <- rle(runs)$lengths
    n.runs           <- length(run.lengths)
    longest.run      <- max(run.lengths)
    longest.run.max  <- round(log2(n.useful)) + 3                # Schilling 2012
    n.crossings      <- max(n.runs - 1, 0)
    n.crossings.min  <- qbinom(0.05, max(n.useful - 1, 0), 0.5)  # Chen 2010 (7)
    runs.signal      <- longest.run > longest.run.max ||
      n.crossings < n.crossings.min
  } else {
    longest.run      <- NA
    longest.run.max  <- NA
    n.crossings      <- NA
    n.crossings.min  <- NA
    runs.signal      <- FALSE
  }

  # run.lengths     <- rle(runs)$lengths
  # n.runs          <- length(run.lengths)
  # longest.run     <- max(run.lengths)
  # longest.run.max <- round(log2(n.useful)) + 3                # Schilling 2012
  # n.crossings     <- max(n.runs - 1, 0)
  # n.crossings.min <- qbinom(0.05, max(n.useful - 1, 0), 0.5)  # Chen 2010 (7)
  # runs.signal     <- longest.run > longest.run.max ||
  #   n.crossings < n.crossings.min
  df$runs.signal  <- runs.signal
  df$longest.run <- longest.run
  df$longest.run.max <- longest.run.max
  df$n.crossings <- n.crossings
  df$n.crossings.min <- n.crossings.min

  return(df)
}

a3 <- function(n) {
  n[n == 0]    <- NA
  tbl          <- c(NA,
                    2.659, 1.954, 1.628, 1.427, 1.287, 1.182,
                    1.099, 1.032, 0.975, 0.927, 0.886, 0.850,
                    0.817, 0.789, 0.763, 0.739, 0.718, 0.698,
                    0.680, 0.663, 0.647, 0.633, 0.619, 0.606)
  x            <- 3 / (4 * (n - 1)) * (4 * n - 3) / sqrt(n)
  w            <- which(n <= 25)
  x[w]         <- tbl[n[w]]
  x[is.nan(x)] <- NA
  return(x)
}

b3 <- function(n) {
  n[n == 0]    <- NA
  tbl          <- c(NA,
                    0.000, 0.000, 0.000, 0.000, 0.030, 0.118,
                    0.185, 0.239, 0.284, 0.321, 0.354, 0.382,
                    0.406, 0.428, 0.448, 0.466, 0.482, 0.497,
                    0.510, 0.523, 0.534, 0.545, 0.555, 0.565)
  x            <- 1 - (3 / c4(n) / sqrt(2 * (n - 1)))
  w            <- which(n <= 25)
  x[w]         <- tbl[n[w]]
  x[is.nan(x)] <- NA
  return(x)
}

b4 <- function(n) {
  n[n == 0]    <- NA
  tbl          <- c(NA,
                    3.267, 2.568, 2.266, 2.089, 1.970, 1.882,
                    1.815, 1.761, 1.716, 1.679, 1.646, 1.618,
                    1.594, 1.572, 1.552, 1.534, 1.518, 1.503,
                    1.490, 1.477, 1.466, 1.455, 1.445, 1.435)
  x            <- 1 + (3 / c4(n) / sqrt(2 * (n - 1)))
  w            <- which(n <= 25)
  x[w]         <- tbl[n[w]]
  x[is.nan(x)] <- NA
  return(x)
}

c4 <- function(n) {
  n[n == 0]   <- NA
  tbl         <- c(NA,
                   0.7979, 0.8862, 0.9213, 0.9400, 0.9515, 0.9594,
                   0.9650, 0.9693, 0.9727, 0.9754, 0.9776, 0.9794,
                   0.9810, 0.9823, 0.9835, 0.9845, 0.9854, 0.9862,
                   0.9869, 0.9876, 0.9882, 0.9887, 0.9892, 0.9896)

  x            <- 4 * (n - 1) / (4 * n - 3)
  w            <- which(n <= 25)
  x[w]         <- tbl[n[w]]
  x[is.nan(x)] <- NA
  return(x)
}

plot.tcc <- function(x,
                     y             = NULL,
                     cex           = 1,
                     pex           = 1,
                     y.expand      = NULL,
                     y.percent     = FALSE,
                     x.date.format = NULL,
                     x.pad         = 1,
                     cl.lab        = TRUE,
                     cl.decimals   = 2,
                     flip          = FALSE,
                     dots.only     = FALSE,
                     ...) {
  df      <- x$df
  main    <- x$main
  ylab    <- x$ylab
  xlab    <- x$xlab
  freeze  <- x$freeze
  col1    <- rgb(093, 165, 218, maxColorValue = 255) # blue
  # col2    <- rgb(255, 165, 000, maxColorValue = 255) # amber
  col2    <- rgb(241, 088, 084, maxColorValue = 255) # red
  col3    <- rgb(140, 140, 140, maxColorValue = 255) # grey
  col4    <- 'white'
  col5    <-  rgb(096, 189, 104, maxColorValue = 255) # green
  cols    <- c('col1' = col1,
               'col2' = col2,
               'col3' = col3,
               'col4' = col4)
  df$pcol <- ifelse(df$limits.signal, 'col2', 'col1')
  df$pcol <- ifelse(df$exclude, 'col4', df$pcol)
  df$lcol <- ifelse(df$runs.signal, 'col2', 'col3')
  p <- ggplot(df) +
    theme_bw(base_size = 11 * cex) +
    theme(panel.border     = element_rect(colour = 'grey70', size = 0.1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # axis.line        = element_blank(),
          axis.ticks       = element_line(colour = col3),
          axis.title.y     = element_text(margin = margin(r = 8)),
          axis.title.x     = element_text(margin = margin(t = 8)),
          plot.title       = element_text(hjust = 0, margin = margin(b = 8)),
          strip.background = element_rect(fill = 'grey90', colour = 'grey70'))
  p <- p +
    geom_line(aes_string(x = 'x', y = 'cl', colour = 'lcol', group = 'breaks'),
              lwd = 0.6 * cex,
              na.rm = TRUE) +
    geom_line(aes_string(x = 'x', y = 'lcl', group = 'breaks'),
              colour = col3,
              lwd = 0.3 * cex,
              na.rm = TRUE) +
    geom_line(aes_string(x = 'x', y = 'ucl', group = 'breaks'),
              colour = col3,
              lwd = 0.3 * cex,
              na.rm = TRUE) +
    geom_line(aes_string(x = 'x', y = 'target'),
              colour = col5,
              lty = 5,
              lwd = 0.3 * cex,
              na.rm = TRUE)

  if(!dots.only) {
    p <- p + geom_line(aes_string(x = 'x', y = 'y', group = 'breaks'),
                       colour = col1,
                       lwd = 1.1 * cex,
                       na.rm = TRUE)
  }

  p <- p + geom_point(aes_string(x = 'x', y = 'y',
                                 group = 'breaks',
                                 fill = 'pcol'),
                      colour = col1,
                      stroke = 0.5,
                      size = 2 * pex * cex,
                      shape = 21,
                      na.rm = TRUE) +
    # Colour centre line according to runs analysis
    scale_colour_manual(values = cols) +
    # Colour data points according to limits analysis
    scale_fill_manual(values = cols) +
    # Suppress legend
    guides(colour = FALSE, fill = FALSE)

  ng1 <- nlevels(droplevels(as.factor(df$g1))) > 1
  ng2 <- nlevels(droplevels(as.factor(df$g2))) > 1

  if(ng1 & ng2) {
    p <- p + facet_grid(g1 ~ g2, ...)
  } else if(ng1) {
    p <- p + facet_wrap(~ g1, ...)
  } else if(ng2) {
    p <- p + facet_wrap(~ g2, ...)
  } else {
    p <- p +
      theme(panel.border = element_blank(),
            axis.line.x = element_line(size = 0.1, colour = col3),
            axis.line.y = element_line(size = 0.1, colour = col3))
  }

  # Add vertical line to mark freeze point
  if(!is.null(freeze)) {
    f <- as.numeric(df$x[freeze])
    f <- f + as.numeric(df$x[freeze + 1] - df$x[freeze]) / 2
    p <- p + geom_vline(xintercept = f, size = 0.5, lty = 3, colour = col3)
  }

  # Format data axis
  if(!is.null(x.date.format)) {
    if(inherits(df$x, 'Date')) {
      p <- p + scale_x_date(date_labels = x.date.format)
    }
    if(inherits(df$x, 'POSIXt')) {
      p <- p + scale_x_datetime(date_labels = x.date.format)
    }
  }


  # Add title and axis labels
  p <- p +
    labs(title = main, x = xlab, y = ylab) +
    expand_limits(y = y.expand)
  # coord_cartesian(ylim = ylim)

  # Add center line label
  if(cl.lab) {
    m <- ifelse(y.percent, 100, 1)
    p <- p + geom_text(aes_string(x = 'tail(x, 1)', y = 'cl',
                                  label = 'paste0("  ", sround(cl * m, cl.decimals))',
                                  hjust = 0),
                       # hjust = -0.25,
                       # nudge_x = 1,
                       check_overlap = TRUE,
                       col = 'grey30',
                       size = 3 * cex)
    if(inherits(df$x, c('Date', 'POSIXt', 'numeric', 'integer'))) {
      p <- p + expand_limits(x = max(df$x) + diff(range(df$x)) * 0.08 * x.pad)
    }
  }

  # Add annotations
  if(sum(!is.na(df$notes))) {
    p <- p +
      geom_text_repel(aes_string(x = 'x', y = 'y', label = 'notes'),
                      size = 3 * cex,
                      # point.padding = unit(2, 'points'),
                      box.padding = unit(0.5, 'lines'),
                      na.rm = TRUE)
  }

  # Flip chart
  if(flip) {
    p <- p + coord_flip()
  }

  if(y.percent) {
    p <- p + scale_y_continuous(labels = scales::percent)
  }

  # Plot chart
  # plot(p)
  return(p)
}

#' Summarise Trellis Control Charts
#'
#' Summary function for tcc objects.
#'
#' @export
#'
#' @param object tcc object
#' @param ... Ignored. Included for compatibility with generic summary function.
#'
#' @return A data frame with summary statistics of the tcc object.
#'
#' @examples
#' # Build data frame for example
#' d <- data.frame(x = rep(1:24, 4),
#'                 mo = (rep(seq(as.Date('2014-1-1'),
#'                               length.out = 24,
#'                               by = 'month'),
#'                           4)),
#'                 n = rbinom(4 * 24, 100, 0.5),
#'                 d = round(runif(4 * 24, 90, 110)),
#'                 g1 = rep(c('a', 'b'), each = 48),
#'                 g2 = rep(c('A', 'B'), each = 24))
#'
#' # P chart
#' p <- tcc(n, d, mo, g1 = g1, g2 = g2, breaks = 12, data = d, chart = 'p')
#' plot(p)
#' summary(p)

summary.tcc <- function(object, ...) {
  d <- object$data
  x1 <- aggregate(n.obs ~ g1 + g2 + breaks,
                  data = d,
                  length)
  x2 <- aggregate(cbind(n.useful = y != cl) ~ g1 + g2 + breaks,
                  data = d[!d$exclude, ],
                  sum,
                  na.rm = TRUE,
                  na.action = na.pass)
  x3 <- aggregate(cbind(cl, lcl, ucl) ~ g1 + g2 + breaks,
                  data = d,
                  tail, 1,
                  na.action = na.pass)
  x4 <- aggregate(cbind(limits.signal, runs.signal) ~ g1 + g2 + breaks,
                  data = d,
                  max,
                  na.rm = TRUE)
  x4[4:5] <- apply(x4[4:5], 2, as.logical)
  x5 <- aggregate(cbind(longest.run, longest.run.max,
                        n.crossings, n.crossings.min) ~ g1 + g2 + breaks,
                  data = d,
                  tail, 1,
                  na.action = na.pass)

  x <- merge(x1, x2)
  x <- merge(x, x3)
  x <- merge(x, x4)
  x <- merge(x, x5)

  return(x)
}
