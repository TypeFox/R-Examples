#' Index zones.
#'
#' Generate a vector of zone "levels" from an input vector and defined
#' boundaries.
#'
#' @param x numeric; values to be "zoned".
#' @param zbounds numeric; values for zone boundaries.
#'
#' @examples
#' data(ridedata)
#'
#' ## Best used to append to existing data.
#' ridedata$zone <- zone_index(ridedata$power.W, c(100, 200, 300))
#'
#' ## How much distance was covered in each zone?
#' ridedata$delta.dist <- c(0, diff(ridedata$distance.km))
#' with(ridedata, tapply(delta.dist, zone, sum, na.rm = TRUE))  # Km.
#'
#' @return a numeric vector of zone values of the same length as \code{x}. The
#'   number of zone levels will be \code{length(zbounds) + 1}.
#'
#' @export
zone_index <- function(x, zbounds) {
  zbounds <- sort(unique(zbounds))    # Important!
  zone_index_(x, zbounds)             # Rcpp.
}
#' Calculate time in zones.
#'
#' Given a vector of zone boundaries, sums the time spent in each zone.
#'
#' @param data a "cycleRdata" object, produced from a \code{\link{read_ride}}
#'   function.
#' @param column the column name of the data to which the zone boundaries
#'   relate.
#' @param zbounds numeric; zone boundaries.
#' @param pct should percentage values be returned?
#' @param character.only are column name arguments given as character strings? A
#'   backdoor around non-standard evaluation. Mainly for internal use.
#'
#' @examples
#' data(ridedata)
#'
#' ## Time spent above and below critical power...
#' zone_time(ridedata, "power.W", zbounds = 300) / 60  # Minutes.
#'
#' ## Or with more zones...
#' zone_time(ridedata, "power.W", zbounds = c(100, 200, 300)) / 60
#'
#' ## Or given as a percentage...
#' zone_time(ridedata, "power.W", zbounds = c(100, 200, 300), pct = TRUE)
#'
#' @return a data frame of zone times.
#'
#' @export
zone_time            <- function(data, column = "power.W", zbounds, pct = FALSE,
                                 character.only = FALSE)
  UseMethod("zone_time", data)
#' @export
zone_time.default    <- function(data, column = "power.W", zbounds, pct = FALSE,
                                 character.only = FALSE)
  format_error()
#' @export
zone_time.cycleRdata <- function(data, column = "power.W", zbounds, pct = FALSE,
                                 character.only = FALSE) {
  if (!character.only)
    column <- as.character(substitute(column))

  data   <- data[data$delta.t <= 20, ]           # Remove sig. stops.
  data$z <- zone_index(data[, column], zbounds)  # Rcpp.
  out    <- tapply(data$delta.t, data$z, sum, na.rm = TRUE)
  # Convert to percentage?
  if (pct) {
    out <- (out / sum(out, na.rm = TRUE)) * 100
    out <- round(out)
  }
  names(out) <- if (length(out) == 1) "Zone 1" else
    paste("Zone", 1:(length(zbounds) + 1))

  out
}
#' Zone-time distribution plot.
#'
#' Display the time distribution of values within a dataset. The distribution
#' can also be partitioned into zones if the \code{zbounds} argument is not
#' \code{NULL}.
#'
#' @param data a "cycleRdata" object, produced from a \code{\link{read_ride}}
#'   function.
#' @param column column in \code{data} giving the values of interest. Needn't be
#'   quoted.
#' @param binwidth how should values in \code{column} be binned? E.g.
#'   \code{bindiwdth = 10} will create 10 watt bins if \code{column} is power
#'   data.
#' @param zbounds optional; a numeric vector of zone boundaries.
#' @param character.only are column name arguments given as character strings? A
#'   backdoor around non-standard evaluation.
#' @param ... arguments to be passed to \code{barplot()} and/or graphical
#'   parameters (\code{\link{par}}).
#'
#' @examples
#' data(ridedata)
#'
#' ## Using power.
#' zdist_plot(
#'  data = ridedata, column = power.W,
#'  binwidth = 10,  # 10 watt bins.
#'  zbounds = c(100, 200, 300),
#'  xlim = c(110, 500), xlab = "Power (Watts)",
#'  main = "Power distribution" # Argument passed to barplot.
#' )
#'
#' ## Using speed.
#' zdist_plot(
#'   data = ridedata, column = speed.kmh,
#'   binwidth = 2,  # 2 km/hr bins.
#'   zbounds = c(10, 20, 30),
#'   xlab = "Speed (km/hr)",
#'   main = "Speed distribution"
#' )
#'
#' ## Without zone colouring (produces a warning).
#' zdist_plot(
#'   data = ridedata, column = speed.kmh,
#'   binwidth = 5,  # 2 km/hr bins.
#'   xlab = "Speed (km/hr)", main = "Dull"
#' )
#'
#' @return nothing; a plot is sent to the current graphics device.
#'
#' @export
zdist_plot            <- function(data, column = "power.W", binwidth = 10, zbounds = NULL,
                                  character.only = FALSE, ...)
  UseMethod("zdist_plot", data)
#' @export
zdist_plot.default    <- function(data, column = "power.W", binwidth = 10, zbounds = NULL,
                                  character.only = FALSE, ...)
  format_error()
#' @export
zdist_plot.cycleRdata <- function(data, column = "power.W", binwidth = 10, zbounds = NULL,
                                  character.only = FALSE, ...) {
  if (!character.only)
    column    <- as.character(substitute(column))

  data_orig <- data                        # Needed later on.
  data      <- data[data$delta.t <= 20, ]  # Remove sig. stops.
  # Convert column of interest to bins.
  data[, column] <- floor(data[, column] / binwidth) * binwidth
  # Calculate percentage time in zbounds.
  t_total <- sum(data$delta.t, na.rm = TRUE)
  vals    <- as.data.frame(tapply(
    data$delta.t, data[, column], FUN = function(x)
    {sum(x, na.rm = TRUE) / t_total} * 100
  ))
  vals[, 2]      <- as.numeric(rownames(vals))
  colnames(vals) <- c("t_pct", column)
  # Plot -----------------------------------------------------------------------
  opar          <- par(no.readonly = TRUE); on.exit(par(opar))
  args          <- list(...)
  args_par      <- args[which(names(args) %in% names(par()))]
  par(args_par)
  # If a zone argument is given.
  if (!is.null(zbounds)) {
    vals$zone    <- factor(zone_index(vals[, column], zbounds))
    # Barplot arguments.
    args_default <- list(height    = vals$t_pct,
                         width     = binwidth,
                         names.arg = names(vals$t_pct),
                         xpd       = FALSE,
                         col       = base_pal(vals$zone, unique(vals$zone)),
                         ylab      = "Time (%)")
    args_usr      <- args[which(names(args) %in% names(formals(barplot.default)))]
    args_barplot  <- c(args_usr,
                       args_default[names(args_default) %notin% names(args_usr)])
    # Plot.
    do.call(barplot, args_barplot)
    # Annotate.
    z_pct     <- zone_time(data_orig, column, zbounds, TRUE, character.only = TRUE)
    z_labs    <- paste0(round(z_pct), "%")
    legend("topright",
           legend    = paste0("Z", unique(vals$zone), ": ", z_labs),
           bty       = "n",
           fill      = base_pal(z_labs, z_labs),
           y.intersp = 1.2)
  } else {  # If no zone argument is given.
    warning("No zone bounds given.", call. = FALSE)
    # Barplot arguments.
    args_default <- list(height    = vals$t_pct,
                         width     = binwidth,
                         names.arg = names(vals$t_pct),
                         xpd       = FALSE,
                         ylab      = "Time (%)")
    args_usr      <- args[which(names(args) %in% names(formals(barplot.default)))]
    args_barplot  <- c(args_usr,
                       args_default[names(args_default) %notin% names(args_usr)])
    # Plot.
    do.call(barplot, args_barplot)
  }
  invisible()
}
