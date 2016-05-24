#' @include cfDataList.R
NULL

#' Plot a windrose
#'
#' Plot a windrose showing the wind speed and direction for given facets using
#' \pkg{ggplot2}.
#'
#' This is intended to be used as a stand-alone function for any wind dataset. A
#' different windrose is plotted for each level of the faceting variable which
#' is coerced to a factor if necessary. The facets will generally be the station
#' where the data were collected, seasons or dates. Currently only one faceting
#' variable is allowed and is passed to \code{\link[ggplot2]{facet_wrap}} with
#' the formula \code{~facet}.
#'
#' @section Theme Selection:
#' For black and white windroses that may be preferred if plots are to be used
#' in journal articles for example, recommended \code{ggtheme}s are \code{'bw'},
#' \code{'linedraw'}, \code{'minimal'} or \code{'classic'} and
#' the \code{col_pal} should be \code{'Greys'}. Otherwise, any of the sequential
#' \code{\link[RColorBrewer]{RColorBrewer}} colour palettes are recommended for
#' colour plots.
#'
#' @return a \code{ggplot} object.
#'
#' @param speed numeric vector of wind speeds.
#' @param direction numeric vector of wind directions.
#' @param facet character or factor vector of the facets used to plot the various
#'              windroses.
#' @param n_directions the number of direction bins to plot (petals on the rose).
#'                     The number of directions defaults to 12.
#' @param n_speeds the number of equally spaced wind speed bins to plot. This is
#'                 used if \code{speed_cuts} is \code{NA} (default 5).
#' @param speed_cuts numeric vector containing the cut points for the wind speed
#'                 intervals, or \code{NA} (default).
#' @param calm_wind the upper limit for wind speed that is considered calm
#'                  (default 0).
#' @param variable_wind numeric code for variable winds (if applicable).
#' @param legend_title character string to be used for the legend title.
#' @param col_pal character string indicating the name of the
#'                \code{\link[RColorBrewer]{RColorBrewer}} colour palette to be
#'                used for plotting, see 'Theme Selection' below.
#' @param ggtheme character string (partially) matching the
#'                \code{\link[ggplot2]{ggtheme}} to be used for plotting, see
#'                'Theme Selection' below.
#' @param n_col The number of columns of plots (default 1).
#' @param ... further arguments passed to \code{\link[ggplot2]{theme}}.
#'
#' @seealso \code{\link[ggplot2]{theme}} for more possible arguments to pass to
#' \code{windrose}.
#'
#' @examples
#' # Create some dummy wind data with predominant south to westerly winds, and
#' # occasional yet higher wind speeds from the NE (not too dissimilar to
#' # Auckland).
#'
#' wind_df = data.frame(wind_speeds = c(rweibull(80, 2, 4), rweibull(20, 3, 9)),
#'                      wind_dirs = c(rnorm(80, 135, 55), rnorm(20, 315, 35)) %% 360,
#'                      station = rep(rep(c("Station A", "Station B"), 2),
#'                                    rep(c(40, 10), each = 2)))
#'
#' # Plot a simple windrose using all the defaults, ignoring any facet variable
#' with(wind_df, windrose(wind_speeds, wind_dirs))
#'
#' # Change the ggtheme and colour scheme for black and white figures
#' with(wind_df, windrose(wind_speeds, wind_dirs,
#'                        ggtheme = "bw",
#'                        col_pal = "Greys"))
#'
#' # Create custom speed bins and legend title
#' with(wind_df, windrose(wind_speeds, wind_dirs,
#'                        speed_cuts = c(3, 6, 9, 12),
#'                        legend_title = "Wind Speed\n(m/s)",
#'                        legend.title.align = .5))
#' # Note that underscore-separated arguments come from the windrose method, and
#' # period-separated arguments come from ggplot2::theme().
#'
#' # Include a facet variable with one level
#' with(wind_df, windrose(wind_speeds, wind_dirs, "Artificial Auckland Wind"))
#'
#' # Plot a windrose for each level of the facet variable (each station)
#' with(wind_df, windrose(wind_speeds, wind_dirs, station, n_col = 2))
#'
#' # Make all the text larger
#' library(ggplot2) # for element_text()
#' with(wind_df, windrose(wind_speeds, wind_dirs, station,
#'                        text = element_text(size = 16), n_col = 2))
#'
#' \dontrun{
#' # Save the plot as a png to the current working directory
#' library(ggplot2)
#' ggsave("my_windrose.png")
#' }
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales percent_format
#' @importFrom ggplot2 ggplot coord_polar geom_bar cut_interval aes_string
#' scale_x_discrete scale_fill_manual theme_grey theme_bw theme_classic
#' theme_gray theme_linedraw theme_light theme_minimal element_blank
#' element_text facet_wrap scale_y_continuous theme
#' @importFrom methods missingArg
#' @export
windrose = function(speed, direction, facet, n_directions = 12,
                    n_speeds = 5, speed_cuts = NA, col_pal = "GnBu",
                    ggtheme = c("grey", "gray", "bw", "linedraw",
                                "light", "minimal", "classic"),
                    legend_title = "Wind Speed", calm_wind = 0,
                    variable_wind = 990, n_col = 1, ...){

  if (missingArg(direction))
    stop("speed can't be missing")

  if (missingArg(direction))
    stop("direction can't be missing")

  include_facet = !missingArg(facet)
  if (include_facet){

    if (!is.character(facet) && !is.factor(facet))
      stop("the faceting variable needs to be character or factor")

    if (length(facet) == 1)
      facet = rep(facet, length(speed))

    if (length(facet) != length(speed))
      stop("the facet variable must be the same length as the wind
                   speeds")
  }

  if (!is.numeric(speed))
    stop("wind speed need to be numeric")

  if (!is.numeric(direction))
    stop("wind directions need to be numeric")

  if (length(speed) != length(direction))
    stop("wind speeds and directions must be the same length")

  if (any(
    (direction > 360 | direction < 0) & (direction != variable_wind),
    na.rm = TRUE)
  )
    stop("wind directions can't be outside the interval [0, 360]")

  if (!is.numeric(n_directions) || length(n_directions) != 1)
    stop("n_directions must be a numeric vector of length 1")

  if (!is.numeric(n_speeds) || length(n_speeds) != 1)
    stop("n_speeds must be a numeric vector of length 1")

  if (!is.numeric(variable_wind) || length(variable_wind) != 1)
    stop("variable_wind must be a numeric vector of length 1")

  if (!is.numeric(calm_wind) || length(calm_wind) != 1)
    stop("calm_wind must be a numeric vector of length 1")

  if (!is.character(legend_title) || length(legend_title) != 1)
    stop("legend title must be a single character string")

  if (n_directions > 180){
    n_directions = 180
    warning("using the maximum number of wind directions; 180")
  }

  if (n_directions < 4){
    n_directions = 4
    warning("using the minimum number of wind directions; 4")
  }

  if (!missing(speed_cuts) && length(speed_cuts) < 3){
    warning("using the minimum 3 speed cuts")
    speed_cuts = 3
  }

  ggtheme = match.arg(ggtheme)

  ## Optimising the input - select values for n_directions so that bins center
  ## on all N, E, S and W
  optimal_n_dir = seq(1, 45, 2) * 4
  if (is.na(match(n_directions, optimal_n_dir))){
    n_directions = optimal_n_dir[which.min(abs(n_directions - optimal_n_dir))]
    message("using the closest optimal number of wind directions (",
            n_directions, ")")
  }

  ## Remove the variable winds
  not_variable = (direction != variable_wind)
  speed = speed[not_variable]
  direction = direction[not_variable]

  if (include_facet)
    facet = facet[not_variable]

  ## Create factor variable for wind direction intervals
  dir_bin_width = 360 / n_directions
  dir_bin_cuts = seq(dir_bin_width / 2, 360 - dir_bin_width / 2, dir_bin_width)
  dir_intervals = findInterval(c(direction, dir_bin_cuts), dir_bin_cuts)
  dir_intervals[dir_intervals == n_directions] = 0
  factor_labs = paste(c(tail(dir_bin_cuts, 1), head(dir_bin_cuts, -1)),
                      dir_bin_cuts, sep = ", ")
  dir_bin = head(factor(dir_intervals, labels = paste0("(", factor_labs, "]")),
                 -n_directions)


  ## Create a factor variable for wind speed intervals
  if (!missing(speed_cuts)){
    if (speed_cuts[1] > min(speed, na.rm = TRUE))
      speed_cuts = c(0, speed_cuts)

    if (tail(speed_cuts, 1) < max(speed, na.rm = TRUE))
      speed_cuts = c(speed_cuts, max(speed, na.rm = TRUE))
    spd_bin = cut(speed, speed_cuts)
  } else
    spd_bin = cut_interval(speed, n_speeds)

  spd_cols = brewer.pal(length(levels(spd_bin)), col_pal)

  if (length(spd_cols) != length(levels(spd_bin)))
    spd_bin = cut_interval(speed, length(spd_cols))

  ## Create the dataframe suitable for plotting
  if (include_facet){
    ggplot_df = as.data.frame(table(dir_bin, spd_bin, facet))
    ggplot_df$proportion = unlist(by(ggplot_df$Freq, ggplot_df$facet,
                                     function(x) x / sum(x)),
                                  use.names = FALSE)
  } else {
    ggplot_df = data.frame(table(dir_bin, spd_bin))
    ggplot_df$proportion = ggplot_df$Freq / sum(ggplot_df$Freq)
  }

  ## (gg)Plot me
  ggwindrose = ggplot(data = ggplot_df,
                      aes_string(x = "dir_bin", fill = "spd_bin", y = "proportion")) +
    geom_bar(stat = "identity") +
    scale_x_discrete(breaks = levels(ggplot_df$dir_bin)[seq(1, n_directions,
                                                            n_directions / 4)],
                     labels = c("N", "E", "S", "W"), drop = FALSE) +
    scale_fill_manual(name = legend_title, values = spd_cols) +
    coord_polar(start = 2 * pi - pi / n_directions) +
    scale_y_continuous(labels = percent_format()) +
    eval(call(paste0("theme_", ggtheme))) +
    theme(axis.title = element_blank(), ...)

  if (include_facet)
    ggwindrose = ggwindrose + facet_wrap(~facet, ncol = n_col)

  return(ggwindrose)
}

#'@importFrom methods setGeneric isGeneric
if (!isGeneric("plot"))
  setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

#' Plot Clifro Wind Objects
#'
#' Various plot methods for exploring wind speed and direction patterns for
#' given CliFlo stations.
#'
#' If \code{x} is a \code{cfDataList}, by default the first datatype will be
#' plotted unless \code{y} is supplied.
#'
#' @note If \code{x} is a \code{cfDataList} object and \code{y} refers to a
#' \pkg{clifro} dataframe that is not a \code{cfWind} object then it will be
#' passed to another method, if available.
#'
#' The default \code{plot} method plots a different windrose for each CliFlo
#' station. The \code{direction_plot} method plots wind direction contours
#' through time to visualise temporal patterns in wind directions. The
#' \code{speed_plot} method plots the time series of wind speeds with a +/-
#' standard deviation region (if applicable).
#'
#' @section Theme Selection:
#' For black and white windroses that may be preferred if plots are to be used
#' in journal articles for example, recommended \code{ggtheme}s are \code{'bw'},
#' \code{'linedraw'}, \code{'minimal'} or \code{'classic'} and
#' the \code{col_pal} should be \code{'Greys'}. Otherwise, any of the sequential
#' \code{\link[RColorBrewer]{RColorBrewer}} colour palettes are recommended for
#' colour plots.
#'
#' @param x a \code{cfWind} or \code{cfDataList} object.
#' @param y missing if \code{x} is a .\code{cfWind} object, otherwise a number
#'          indicating the dataframe to plot in the \code{cfDataList} (defaults
#'          to 1).
#' @param n_directions the number of direction bins to plot (petals on the
#'                     rose). The number of directions defaults to 12.
#' @param n_speeds the number of equally spaced wind speed bins to plot. This is
#'                 used if \code{spd_cuts} is NA (default 5).
#' @param speed_cuts numeric vector containing the cut points for the wind speed
#'                 intervals, or NA (default).
#' @param col_pal character string indicating the name of the
#'                \code{\link[RColorBrewer]{RColorBrewer}} colour palette to be
#'                used for plotting, see 'Theme Selection' below.
#' @param ggtheme character string (partially) matching the
#'                \code{\link[ggplot2]{ggtheme}} to be used for plotting, see
#'                'Theme Selection' below.
#' @param n_col the number of columns of plots (default 1).
#' @param contours the number of contour lines to draw (default 10).
#' @param scales character string partially matching the \code{scales} argument
#'               in the \code{link[ggplot2]{facet_wrap}} function.
#' @param ... further arguments passed to \code{\link[ggplot2]{theme}}.
#'
#' @note Given a value on the x-axis, the ends of the density function along the
#'  y-axis are not constrained to be equal for any of the derivatives for the
#'  \code{direction_plot} method. That is, the contours at direction = 0, do not
#'  match the contours at direction = 360.
#'
#'  @seealso \code{\link{plot,cfDataList,missing-method}} for general
#'   information on default plotting of \code{cfData} and \code{cfDataList}
#'   objects, and the links within. See \code{\link{cf_query}} for creating
#'   \code{cfWind} objects or \code{\link{windrose}} for plotting any wind data.
#'   Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
#'   to these methods. \code{\link{summary,cfWind-method}} for summarising wind
#'   information at each CliFlo station.
#'
#' @aliases plot.cfWind
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales percent_format
#' @importFrom ggplot2 ggplot coord_polar geom_bar cut_interval aes_string
#' scale_x_discrete scale_fill_manual theme_grey theme_bw theme_classic
#' theme_gray theme_linedraw theme_light theme_minimal element_blank
#' element_text
#' @examples
#' \dontrun{
#' # Retrieve maximum wind gust data at the Reefton Ews station from CliFlo
#' # (public data)
#' reefton_wind = cf_query(cf_user(), cf_datatype(2, 2, 1, 1), cf_station(),
#'                         start_date = "2012-01-01-00")
#'
#' class(reefton_wind)
#'
#' # Examples of the default plots --------------------------------------------
#'
#' # Plot a windrose
#' plot(reefton_wind)
#'
#' # Plot the wind direction contours
#' direction_plot(reefton_wind)
#'
#' # Plot the wind speed time-series
#' speed_plot(reefton_wind)
#'
#' # Examples of changing the defaults ----------------------------------------
#'
#' # Plot black and white windroses
#' plot(reefton_wind, ggtheme = "bw", col_pal = "Greys")
#' plot(reefton_wind, ggtheme = "linedraw", col_pal = "Greys")
#' plot(reefton_wind, ggtheme = "classic", col_pal = "Greys")
#' plot(reefton_wind, ggtheme = "minimal", col_pal = "Greys")
#'
#' # Plot the wind directions using 20 contours and the ggtheme 'classic'
#' direction_plot(reefton_wind, ggtheme = "classic", contours = 20)
#'
#' # Enlarge all the text to 18pt
#' library(ggplot2) # for element_text() and geom_point()
#' direction_plot(reefton_wind, ggtheme = "classic", contours = 20,
#'                text = element_text(size = 18))
#'
#' # Include the actual observations in the plots
#' direction_plot(reefton_wind) + geom_point(alpha = .2, size = 3)
#'
#' speed_plot(reefton_wind, ggtheme = "classic", text = element_text(size = 16)) +
#'   geom_point(shape = 1, size = 3)
#' # or equivalently using base graphics:
#' plot(reefton_wind$Date, reefton_wind$Speed, type = 'o',
#'      xlab = NA, ylab = "Daily max gust (m/s)", las = 1, main = "Reefton Ews")
#'
#' # Example of plotting a cfDataList -----------------------------------------
#' # Collect both surface wind run and hourly surface wind observations from
#' # Reefton Ews
#' reefton_list = cf_query(cf_user(), cf_datatype(2, 1, 1:2, 1),
#'                         cf_station(), "2012-01-01 00", "2012-02-01 00")
#'
#' reefton_list
#'
#' class(reefton_list) #cfDataList
#'
#' # Plot the first (default) dataframe
#' plot(reefton_list) # Error - no wind directions for wind run datatypes
#' # Try speed_plot instead
#' speed_plot(reefton_list)
#'
#' # Plot the second dataframe in the cfDataList
#' plot(reefton_list, 2)           # identical to plot(reefton_list[2])
#' speed_plot(reefton_list, 2)     # identical to speed_plot(reefton_list[2])
#' direction_plot(reefton_list, 2) # identical to direction_plot(reefton_list[2])
#'
#' # Save the ggplot externally -----------------------------------------------
#'
#' # Save the plot as a png to the current working directory
#' library(ggplot2) # for ggsave()
#' ggsave("my_wind_plot.png")
#' }
#' @export
setMethod("plot",
          signature(x = "cfWind", y = "missing"),
          definition = function(x, y, n_directions = 12, n_speeds = 5,
                                speed_cuts = NULL, col_pal = "GnBu",
                                ggtheme = c("grey", "gray", "bw", "linedraw",
                                            "light", "minimal", "classic"),
                                n_col = 1, ...)
          {

            if (all(is.na(x[, 4])))
              stop("no wind speeds to plot - try direction_plot(x)")

            if (all(is.na(x[, 3])))
              stop("no wind directions to plot - try speed_plot(x)")

            wind_df = data.frame(facet = x[, 1], speed = x[, 4],
                                 dir = x[, 3])

            if (!is.character(col_pal) || length(col_pal) != 1)
              stop("col_pal must be a single character string")

            if (!is.numeric(n_directions) || length(n_directions) != 1)
              stop("n_directions must be a single number")

            if (!is.numeric(n_speeds) || length(n_speeds) != 1)
              stop("n_speeds must be a single number")

            if (n_directions > 180) {
              n_directions = 180
              message("using the maximum number of wind directions; 180")
            }

            if (n_directions < 4) {
              n_directions = 4
              message("using the minimum number of wind directions; 4")
            }

            if (!is.null(speed_cuts) && length(speed_cuts) < 3) {
              message("using the minimum 3 speed cuts")
              speed_cuts = 3
            }

            ggtheme = match.arg(ggtheme)

            ## Optimising the input - select values for n_directions so that
            ## bins center on all N, E, S and W
            optimal_n_dir = seq(1, 45, 2) * 4
            if (is.na(match(n_directions, optimal_n_dir))) {
              n_directions = optimal_n_dir[which.min(abs(n_directions -
                                                           optimal_n_dir))]
              message("using the closest optimal number of wind directions (",
                      n_directions, ")")
            }

            wind_df = subset(wind_df, dir != 990)
            dir_bin_width = 360/n_directions
            dir_bin_cuts = seq(dir_bin_width/2, 360 - dir_bin_width/2,
                               dir_bin_width)
            dir_intervals = findInterval(c(wind_df$dir, dir_bin_cuts),
                                         dir_bin_cuts)
            dir_intervals[dir_intervals == n_directions] = 0
            factor_labs = paste(c(tail(dir_bin_cuts, 1),
                                  head(dir_bin_cuts, -1)), dir_bin_cuts,
                                sep = ", ")
            dir_bin = head(factor(dir_intervals,
                                  labels = paste0("(", factor_labs, "]")),
                           -n_directions)

            if (!missing(speed_cuts)) {

              if (speed_cuts[1] > min(wind_df$speed, na.rm = TRUE))
                speed_cuts = c(0, speed_cuts)

              if (tail(speed_cuts, 1) < max(wind_df$speed, na.rm = TRUE))
                speed_cuts = c(speed_cuts, max(wind_df$speed, na.rm = TRUE))

              speed_bin = cut(wind_df$speed, speed_cuts)
            } else
              speed_bin = cut_interval(wind_df$speed, n_speeds)

            speed_cols = brewer.pal(length(levels(speed_bin)), col_pal)

            if (length(speed_cols) != length(levels(speed_bin)))
              speed_bin = cut_interval(wind_df$speed, length(speed_cols))

            ggplot_df = as.data.frame(table(dir_bin, speed_bin, wind_df$facet))
            names(ggplot_df)[3] = "facet"
            ggplot_df$proportion = unlist(by(ggplot_df$Freq, ggplot_df$facet,
                                             function(x) x / sum(x)),
                                          use.names = FALSE)
            ggplot(data = ggplot_df, aes_string(x = "dir_bin",
                                                fill = "speed_bin",
                                                y = "proportion")) +
              geom_bar(stat = "identity") +
              scale_x_discrete(breaks =
                                 levels(dir_bin)[seq(1, n_directions,
                                                     n_directions/4)],
                               labels = c("N", "E", "S", "W"), drop = FALSE) +
              scale_fill_manual(name = x@data_label, values = speed_cols) +
              coord_polar(start = 2 * pi - pi/n_directions) +
              facet_wrap(~facet, ncol = n_col) +
              scale_y_continuous(labels = percent_format()) +
              eval(call(paste0("theme_", ggtheme))) +
              theme(axis.title = element_blank(), ...)
          }
)

#'@importFrom methods setGeneric
setGeneric("direction_plot", function(x, y, ...) standardGeneric("direction_plot"))

#' @importFrom ggplot2 ggplot stat_density2d scale_y_continuous facet_wrap
#' scale_alpha_continuous theme aes_string theme_grey theme_bw
#' theme_classic theme_gray theme_linedraw theme_light theme_minimal
#' @importFrom methods setMethod
#' @rdname plot-cfWind-missing-method
#' @aliases direction_plot
#' @export
setMethod("direction_plot",
          signature(x = "cfWind", y = "missing"),
          definition = function (x, y, ggtheme = c("grey", "gray", "bw", "linedraw",
                                                "light", "minimal", "classic"),
                                 contours = 10, n_col = 1, ...)
          {
            n_col = n_col[[1]]
            contours = contours[[1]]

            if (!is.numeric(n_col))
              stop("number of columns must be a single number")

            if (!is.numeric(contours))
              stop("contours must be a single number")

            ggtheme = match.arg(ggtheme)

            if (all(is.na(x[, 3])))
              stop("no wind speeds to plot - try speed_plot(x)")

            df = data.frame(date = x[, 2], Direction = x[, 3], facet = x[, 1])

            ggplot(df, aes_string(x = "date", y = "Direction")) +
              stat_density2d(aes_string(alpha = "..level.."), bins = contours,
                             size = 1) +
              facet_wrap(~facet, ncol = n_col) +
              eval(call(paste0("theme_", ggtheme))) +
              scale_alpha_continuous(guide = FALSE) +
              scale_y_continuous(breaks = seq(0, 360, by = 90),
                                 labels = c("N", "E", "S", "W", "N")) +
              theme(axis.title.x = element_blank(), ...)
          }
)

#' @importFrom methods setMethod
#' @rdname plot-cfWind-missing-method
#' @export
setMethod("direction_plot",
          signature(x = "cfDataList", y = "missing"),
          definition = function (x, y, ...){
            if (is(x[1], "cfWind"))
              direction_plot(x[1], ...)
            else
              stop("no direction_plot method for class ", class(x[1])[[1]])
          }
)

#' @importFrom methods setMethod
#' @rdname plot-cfWind-missing-method
#' @export
setMethod("direction_plot",
          signature(x = "cfDataList", y = "numeric"),
          definition = function (x, y, ...){
            y = y[1]
            if (y > length(x))
              stop("y must be between 1 and ", length(x))

            if (is(x[y], "cfWind"))
              direction_plot(x[y], ...)
            else
              stop("no direction_plot method for class ", class(x[y])[[1]])
          }
)

#'@importFrom methods setGeneric
setGeneric("speed_plot", function(x, y, ...) standardGeneric("speed_plot"))

#' @importFrom ggplot2 ggplot geom_ribbon geom_line facet_wrap theme aes ylab
#' element_blank theme_grey theme_bw theme_classic theme_gray theme_linedraw
#' theme_light theme_minimal
#' @importFrom methods setMethod
#' @rdname plot-cfWind-missing-method
#' @aliases speed_plot
#' @export
setMethod("speed_plot",
          signature(x = "cfWind", y = "missing"),
          definition = function (x, y,
                                 ggtheme = c("grey", "gray", "bw", "linedraw",
                                             "light", "minimal", "classic"),
                                 scales = c("fixed", "free_x", "free_y", "free"),
                                 n_col = 1, ...)
          {
            include_sd = if (tolower(x@dt_name) == "max gust") FALSE else TRUE

            if (include_sd && any(is.na(x[, 6]))) {
              include_sd = FALSE
              message("unable to plot wind speed standard deviations due to NA's")
            }

            if (include_sd){
              df = data.frame(date = x[, 2], speed = x[, 4],
                              spd.sd = x[, 6], facet = x[, 1],
                              lower = x[, 4] - x[, 6],
                              upper = x[, 4] + x[, 6])
            } else {
              df = data.frame(date = x[, 2], speed = x[, 4],
                              spd.sd = x[, 6], facet = x[, 1])
            }

            if (!is.numeric(n_col) || length(n_col) != 1)
              stop("number of columns must be a single number")

            ggtheme = match.arg(ggtheme)
            scales = match.arg(scales)

            p = ggplot(data = df, aes_string(x = "date", y = "speed"))

            if (include_sd)
              p = p + geom_ribbon(aes_string(ymin = "lower",
                                      ymax = "upper"),
                                  alpha = 0.2)

            p = p + geom_line() +
              eval(call(paste0("theme_", ggtheme))) +
              ylab(x@data_label) +
              facet_wrap(~facet, scales = scales, ncol = n_col) +
              theme(axis.title.x = element_blank(), ...)

            return(p)
          }
)

#' @importFrom methods setMethod
#' @rdname plot-cfWind-missing-method
#' @export
setMethod("speed_plot",
          signature(x = "cfDataList", y = "missing"),
          definition = function (x, y, ...){
            if (is(x[1], "cfWind"))
              speed_plot(x[1], ...)
            else
              stop("no speed_plot method for class ", class(x[1])[[1]])
          }
)

#' @importFrom methods setMethod
#' @rdname plot-cfWind-missing-method
#' @export
setMethod("speed_plot",
          signature(x = "cfDataList", y = "numeric"),
          definition = function (x, y, ...){
            y = y[1]
            if (y > length(x))
              stop("y must be between 1 and ", length(x))

            if (is(x[y], "cfWind"))
              speed_plot(x[y], ...)
            else
              stop("no speed_plot method for class ", class(x[y])[[1]])
          }
)

#' Summarise Clifro Wind Data
#'
#' This is a summary method for \code{cfWind} objects.
#'
#' A dataframe is returned containing the percentage of calm days
#' (wind speed >= \code{calm_days}), percentage of variable days (wind speed =
#' 990), and quantiles from the empirical cumulative distribution functions for
#' each CliFlo station at which there is wind data.
#'
#' @param object a \code{cfWind} object.
#' @param calm_wind a single number containing the wind speed that is considered
#'                  calm.
#'
#' @importFrom stats quantile
#' @importFrom methods setMethod
#' @examples
#' \dontrun{
#' # Retrieve maximum wind gust data at the Reefton Ews station from CliFlo
#' # (public data)
#' reefton_wind = cf_query(cf_user(), cf_datatype(2, 2, 1, 1), cf_station(),
#'                         start_date = "2012-01-01-00")
#'
#' class(reefton_wind) # cfWind object
#'
#' # Summarise the information
#' summary(reefton_wind)
#' }
#' @seealso \code{\link{plot.cfWind}} for default plotting of
#'  clifro wind data, and \code{\link{cf_query}} for creating \code{cfWind}
#'  objects.
#' @export
setMethod("summary",
          "cfWind",
          function (object, calm_wind = 0)
          {
            if (!is.numeric(calm_wind) || length(calm_wind) != 1)
              stop("calm_wind must be a single number")
            summary_df = data.frame(Station = object[, 1], Wind = object[, 4],
                                    Direction = object[, 3])
            calm_days = with(summary_df, by(Wind, Station, FUN = function(x)
              sum(x <= calm_wind)/length(x)))
            variable_days = with(summary_df, by(Wind, Station, FUN = function(x)
              sum(x == 990)/length(x)))
            speed_ecdf = with(summary_df, do.call(rbind, by(Wind, Station,
                                                     FUN = function(x)
                                                       quantile(ecdf(x)))))
            colnames(speed_ecdf) = seq(0, 100, by = 25)
            round(cbind(speed_ecdf, calm = calm_days * 100,
                        variable = variable_days * 100), 1)
          }
)

# Precipitation -----------------------------------------------------------

#' Plot Rain Timeseries
#'
#' Plot the amount of rainfall (mm) through time, with optional available soil
#' water capacity and runoff amounts (if applicable).
#'
#' When there is a rain event, the amount of runoff, if any, is dependent on how
#' much capacity the soil has available for more water. If there is no available
#' water capacity left in the soil then more rain will lead to a runoff event.
#' If \code{include_runoff = TRUE}, the available water capacity is plotted as
#' negative values and the runoff as positive values to signify this negative
#' relationship.
#'
#' @param x a \code{cfRain} object.
#' @param y missing.
#' @param include_runoff a logical indicating whether to plot the soil moisture
#'                       deficit and runoff as well as the rainfall, if the data
#'                       is available (default \code{TRUE}).
#' @param ggtheme character string (partially) matching the
#'                \code{\link[ggplot2]{ggtheme}} to be used for plotting, see
#'                'Theme Selection' below.
#' @param scales character string partially matching the \code{scales} argument
#'               in the \code{link[ggplot2]{facet_wrap}} function.
#' @param n_col the number of columns of plots (default 1).
#' @param ... further arguments passed to \code{\link[ggplot2]{theme}}.
#'
#' @importFrom ggplot2 ggplot aes_string facet_wrap geom_ribbon
#' scale_fill_discrete ylab geom_line geom_point scale_colour_manual theme
#' theme_grey theme_bw theme_classic theme_gray theme_linedraw theme_light
#' theme_minimal
#' @importFrom methods setMethod
#' @aliases plot.cfRain
#' @seealso \code{\link{plot,cfDataList,missing-method}} for general
#'   information on default plotting of \code{cfData} and \code{cfDataList}
#'   objects, and the links within. See \code{\link{cf_query}} for creating
#'   \code{cfRain} objects.
#'
#'   Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
#'   to these methods.
#' @examples
#' \dontrun{
#' # Retrieve public rain data for a month from CliFlo (at Reefton Ews station)
#' reefton_rain = cf_query(cf_user(), cf_datatype(3, 1, 1), cf_station(),
#'                         start_date = "2012-08-01-00",
#'                         end_date = "2012-09-01-00")
#'
#' class(reefton_rain) # cfRain object
#'
#' # Plot the rain data using the defaults
#' plot(reefton_rain)
#'
#' # Change the ggtheme and enlarge the text
#' library(ggplot2) # for element_text()
#' plot(reefton_rain, ggtheme = "bw", text = element_text(size = 16))
#'
#' # Save the plot as a png to the current working directory
#' library(ggplot2) # for ggsave()
#' ggsave("my_rain_plot.png")
#' }
#' @export
setMethod("plot",
          signature(x = "cfRain", y = "missing"),

          definition = function (x, y, include_runoff = TRUE,
                                 ggtheme = c("grey", "gray", "bw", "linedraw",
                                             "light", "minimal", "classic"),
          scales = c("fixed", "free_x", "free_y", "free"), n_col = 1,
          ...)
{
  include_runoff = include_runoff[[1]]
  n_col = n_col[[1]]

  if (!is.logical(include_runoff))
    stop("include_runoff must be either TRUE or FALSE")

  if (!is.numeric(n_col))
    stop("n_col must be a single number")

  ggtheme = match.arg(ggtheme)
  scales = match.arg(scales)

  if(!pmatch("Deficit", names(x), 0) && include_runoff){
    message("no runoff data to plot")
    include_runoff = FALSE
  }

  df = as(x, "data.frame")
  names(df)[1:3] = c("station", "date", "amount")

  if (include_runoff) {
    names(df)[5:6] = c("deficit", "runoff")
    df$deficit = -df$deficit
    df = melt(df[, c("station", "date", "amount", "runoff",
                     "deficit")], id = c("date", "station"))
    df$variable = factor(df$variable, labels = c("Rain",
                                                 "Soil runoff", "Soil deficit (AWC)"))
  }

  df$zero = 0

  p = ggplot(df, aes_string(x = "date"))

  if (include_runoff)
    p = p + geom_ribbon(aes_string(ymin = "zero", ymax = "value",
                                   fill = "variable"), alpha = 0.5) +
    ylab("Amount (mm)")
  else
    p = p + geom_ribbon(aes_string(ymin = "zero", ymax = "amount")) +
    geom_line(aes_string(y = "amount")) +
    ylab("Rain (mm)")
  p = p + facet_wrap(~station, scales = scales, ncol = n_col) +
    eval(call(paste0("theme_", ggtheme))) +
    theme(axis.title.x = element_blank(), legend.title = element_blank(), ...)

  return(p)
}
)

# Temperature -------------------------------------------------------------

#' Plot Screen Observations
#'
#' Plot temperature data from screen observations (degrees celsuis) through time.
#'
#' Temperature data from screen observations include the air, and wet bulb,
#' temperature at the time the measurement was taken (dry bulb and wet bulb
#' respectively), and the dew point. The dew point is the air temperature at
#' which dew starts to form. That is the temperature to which a given air parcel
#' must be cooled at constant pressure and constant water vapour content in
#' order for saturation to occur.
#'
#' The resulting figure plots the dry bulb, wet bulb and dew point temperatures
#' on the same  scale, for each station.
#'
#' @references \href{http://cliflo.niwa.co.nz/pls/niwp/wh.do_help?id=ls_scr1}{Screen Observation details}.
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_string geom_line ylab theme facet_wrap
#' theme_grey theme_bw theme_classic theme_gray theme_linedraw theme_light
#' theme_minimal
#'
#' @param x a cfScreenObs object.
#' @param y missing.
#' @param ggtheme character string (partially) matching the
#'                \code{\link[ggplot2]{ggtheme}} to be used for plotting, see
#'                'Theme Selection' below.
#' @param scales character string partially matching the \code{scales} argument
#'               in the \code{link[ggplot2]{facet_wrap}} function.
#' @param n_col the number of columns of plots (default 1).
#' @param ... further arguments passed to \code{\link[ggplot2]{theme}}.
#'
#' @aliases plot.cfScreenObs
#' @seealso \code{\link{plot,cfDataList,missing-method}} for general
#'   information on default plotting of \code{cfData} and \code{cfDataList}
#'   objects, and the links within. See \code{\link{cf_query}} for creating
#'   \code{cfScreenObs} objects.
#'
#'   Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
#'   to these methods.
#' @examples
#' \dontrun{
#' # Retrieve public temperature data from screen observations for the last week
#' # at Reefton Ews station
#'
#' # Subtract 7 days from today's date to get the start date
#' last_week = paste(as.character(Sys.Date() - 7), 0)
#'
#' reefton_screenobs = cf_query(cf_user(), cf_datatype(4, 1, 1), cf_station(),
#'                              start_date = last_week)
#'
#' class(reefton_screenobs) # cfScreenObs object
#'
#' # Plot the temperature data using the defaults
#' plot(reefton_screenobs)
#'
#' # Enlarge the text and add the observations as points
#' library(ggplot2) # for element_text() and geom_point()
#' plot(reefton_screenobs, ggtheme = "bw", text = element_text(size = 16)) +
#'   geom_point(size = 3, shape = 1)
#'
#' # Save the plot as a png to the current working directory
#' library(ggplot2) # for ggsave()
#' ggsave("my_screenobs_plot.png")
#' }
#' @export
setMethod("plot",
          signature(x = "cfScreenObs", y = "missing"),

          definition = function (x, y,
                                 ggtheme = c("grey", "gray", "bw", "linedraw",
                                             "light", "minimal", "classic"),
                                 scales = c("fixed", "free_x", "free_y", "free"),
                                 n_col = 1, ...)
          {
            if (!is.numeric(n_col))
              stop("n_col must be a single number")

            ggtheme = match.arg(ggtheme)
            scales = match.arg(scales)

            x.df = as.data.frame(x)
            names(x.df) = c("station", "date", "Dry bulb",
                            "Wet bulb", "rh", "Dew point")

            gg_df = melt(x.df, id.vars = 1:2,
                         measure.vars = c(3, 4, 6))

            p = ggplot(gg_df, aes_string(x = "date", y = "value",
                                         colour = "variable")) +
              facet_wrap(~station, scales = scales, ncol = n_col) +
              geom_line(size = 1, alpha = .7) +
              ylab(x@plot_label) +
              eval(call(paste0("theme_", ggtheme))) +
              theme(axis.title.x = element_blank(),
                    legend.title = element_blank(),
                    ...)
            return(p)
          }
)

#' Plot Temperature Range
#'
#' Plot minimum and maximum temperature data for a given period (degrees
#' celsuis) through time, for each chosen CliFlo station.
#'
#' This plotting method shows the temperature extremes as a grey region on the
#' plot, with a black line indicating the average temperature (if available).
#'
#' @importFrom ggplot2 ggplot geom_ribbon aes_string facet_wrap geom_line ylab
#' theme theme_grey theme_bw theme_classic theme_gray theme_linedraw theme_light
#' theme_minimal aes
#'
#' @param x a cfTemp object.
#' @param y missing.
#' @param ggtheme character string (partially) matching the
#'                \code{\link[ggplot2]{ggtheme}} to be used for plotting, see
#'                'Theme Selection' below.
#' @param scales character string partially matching the \code{scales} argument
#'               in the \code{link[ggplot2]{facet_wrap}} function.
#' @param n_col the number of columns of plots (default 1).
#' @param ... further arguments passed to \code{\link[ggplot2]{theme}}.
#'
#' @aliases plot.cfTemp
#' @seealso \code{\link{plot,cfDataList,missing-method}} for general
#'   information on default plotting of \code{cfData} and \code{cfDataList}
#'   objects, and the links within. See \code{\link{cf_query}} for creating
#'   \code{cfTemp} objects.
#'
#'   Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
#'   to these methods.
#' @examples
#' \dontrun{
#' # Retrieve public hourly minimum and maximum temperature data for the last
#' week at Reefton Ews station
#'
#' # Subtract 7 days from today's date to get the start date
#' last_week = paste(as.character(Sys.Date() - 7), 0)
#'
#' reefton_temp = cf_query(cf_user(), cf_datatype(4, 2, 2), cf_station(),
#'                         start_date = last_week)
#'
#' class(reefton_temp) # cfTemp object
#'
#' # Plot the temperature data using the defaults
#' plot(reefton_temp)
#'
#' # Enlarge the text and add the observations as points
#' library(ggplot2) # for element_text() and geom_point()
#' plot(reefton_temp, ggtheme = "bw", text = element_text(size = 16)) +
#'   geom_point(size = 3, shape = 1)
#'
#' # Save the plot as a png to the current working directory
#' library(ggplot2) # for ggsave()
#' ggsave("my_temperature_plot.png")
#' }
#' @export
setMethod("plot",
          signature(x = "cfTemp", y = "missing"),

          definition = function (x, y,
                                 ggtheme = c("grey", "gray", "bw", "linedraw",
                                             "light", "minimal", "classic"),
                                 scales = c("fixed", "free_x", "free_y", "free"),
                                 n_col = 1, ...){
            if (!is.numeric(n_col))
              stop("n_col must be a single number")

            ggtheme = match.arg(ggtheme)
            scales = match.arg(scales)

            x_df = as(x, "data.frame")
            names(x_df)[c(2, 3, 5, 9)] = c("date", "max", "min", "mean")

            if (!all(is.na(x_df[, 9]))){
              p = ggplot(x_df, aes_string(x = "date", y = "mean")) +
                geom_ribbon(aes_string(ymax = "max", ymin = "min"), alpha = .3) +
                geom_line()
            } else {
              p = ggplot(x_df, aes_string(x = "date")) +
                geom_ribbon(aes_string(ymax = "max", ymin = "min"), alpha = .3)
            }

            p = p +
              ylab(x@plot_label) +
              facet_wrap(~Station, scales = scales, ncol = n_col) +
              eval(call(paste0("theme_", ggtheme))) +
              theme(axis.title.x = element_blank(),
                    legend.title = element_blank(),
                    ...)
            return(p)
          }
)

#' Plot Earth Temperatures
#'
#' Plot the earth temperature for a given depth (degrees celsuis) through time,
#' for each chosen CliFlo station.
#'
#' @importFrom ggplot2 ggplot aes_string geom_line ylab theme element_blank
#' theme_grey theme_bw theme_classic theme_gray theme_linedraw theme_light
#' theme_minimal facet_wrap ylab
#'
#' @param x a cfEarthTemp object.
#' @param y missing.
#' @param ggtheme character string (partially) matching the
#'                \code{\link[ggplot2]{ggtheme}} to be used for plotting, see
#'                'Theme Selection' below.
#' @param scales character string partially matching the \code{scales} argument
#'               in the \code{link[ggplot2]{facet_wrap}} function.
#' @param n_col the number of columns of plots (default 1).
#' @param ... further arguments passed to \code{\link[ggplot2]{theme}}.
#'
#' @aliases plot.cfEarthTemp
#' @seealso \code{\link{plot,cfDataList,missing-method}} for general
#'   information on default plotting of \code{cfData} and \code{cfDataList}
#'   objects, and the links within. See \code{\link{cf_query}} for creating
#'   \code{cfEarthTemp} objects.
#'
#'   Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
#'   to these methods.
#' @examples
#' \dontrun{
#' # Retrieve public earth temperature data for the last 30 days at Reefton Ews
#' # station, at a depth of 10cm
#'
#' # Subtract 30 days from today's date to get the start date
#' last_month = paste(as.character(Sys.Date() - 30), 0)
#'
#' reefton_earth = cf_query(cf_user(), cf_datatype(4, 3, 2), cf_station(),
#'                          start_date = last_month)
#'
#' class(reefton_earth) # cfTemp object
#'
#' # Plot the temperature data using the defaults
#' plot(reefton_earth)
#'
#' # Enlarge the text and add the observations as points
#' library(ggplot2) # for element_text() and geom_point()
#' plot(reefton_earth, ggtheme = "bw", text = element_text(size = 16)) +
#'   geom_point(size = 3, shape = 1)
#'
#' # Save the plot as a png to the current working directory
#' library(ggplot2) # for ggsave()
#' ggsave("my_earthTemp_plot.png")
#' }
#' @export
setMethod("plot",
          signature(x = "cfEarthTemp", y = "missing"),

          definition = function (x, y,
                                 ggtheme = c("grey", "gray", "bw", "linedraw",
                                             "light", "minimal", "classic"),
                                 scales = c("fixed", "free_x", "free_y", "free"),
                                 n_col = 1, ...){
            if (!is.numeric(n_col))
              stop("n_col must be a single number")

            ggtheme = match.arg(ggtheme)
            scales = match.arg(scales)

            x_df = as(x, "data.frame")
            names(x_df)[2:3] = c("date", "temp")

            p = ggplot(x_df, aes_string(x = "date", y = "temp")) +
              geom_line() +
              eval(call(paste0("theme_", ggtheme))) +
              theme(axis.title.x = element_blank(), ...) +
              ylab(x@plot_label) +
              facet_wrap(~Station, scales = scales, ncol = n_col)

            return(p)
          }
)

# Sunshine ----------------------------------------------------------------


#' Plot Sunshine Hours
#'
#' Plot the duration of accumulated bright sunshine hours through time.
#'
#' @importFrom ggplot2 ggplot aes_string geom_line ylab theme element_blank
#' theme_grey theme_bw theme_classic theme_gray theme_linedraw theme_light
#' theme_minimal facet_wrap
#'
#' @param x a cfSunshine object.
#' @param y missing.
#' @param ggtheme character string (partially) matching the
#'                \code{\link[ggplot2]{ggtheme}} to be used for plotting, see
#'                'Theme Selection' below.
#' @param scales character string partially matching the \code{scales} argument
#'               in the \code{link[ggplot2]{facet_wrap}} function.
#' @param n_col the number of columns of plots (default 1).
#' @param ... further arguments passed to \code{\link[ggplot2]{theme}}.
#'
#' @aliases plot.cfSunshine
#' @seealso \code{\link{plot,cfDataList,missing-method}} for general
#'   information on default plotting of \code{cfData} and \code{cfDataList}
#'   objects, and the links within. See \code{\link{cf_query}} for creating
#'   \code{cfSunshine} objects.
#'
#'   Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
#'   to these methods.
#' @examples
#' \dontrun{
#' # Retrieve public hourly sunshine data for the last 7 days at Reefton Ews
#' # station
#'
#' # Subtract 7 days from today's date to get the start date
#' last_week = paste(as.character(Sys.Date() - 7), 0)
#'
#' reefton_sun = cf_query(cf_user(), cf_datatype(5, 1, 2), cf_station(),
#'                        start_date = last_week)
#'
#' class(reefton_sun) # cfSunshine object
#'
#' # Plot the temperature data using the defaults
#' plot(reefton_sun)
#'
#' # Enlarge the text and add the observations as points
#' library(ggplot2) # for element_text() and geom_point()
#' plot(reefton_sun, ggtheme = "bw", text = element_text(size = 16)) +
#'   geom_point(size = 3, shape = 1)
#'
#' # Save the plot as a png to the current working directory
#' library(ggplot2) # for ggsave()
#' ggsave("my_sunshine_plot.png")
#' }
#' @export
setMethod("plot",
          signature(x = "cfSunshine", y = "missing"),

          definition = function (x, y,
                                 ggtheme = c("grey", "gray", "bw", "linedraw",
                                             "light", "minimal", "classic"),
                                 scales = c("fixed", "free_x", "free_y", "free"),
                                 n_col = 1, ...){
            if (!is.numeric(n_col))
              stop("n_col must be a single number")

            ggtheme = match.arg(ggtheme)
            scales = match.arg(scales)

            x_df = as(x, "data.frame")
            names(x_df)[2:3] = c("date", "amount")

            p = ggplot(x_df, aes_string(x = "date", y = "amount")) +
              geom_line() +
              eval(call(paste0("theme_", ggtheme))) +
              ylab(x@data_label) +
              facet_wrap(~Station, scales = scales, ncol = n_col) +
              theme(axis.title.x = element_blank(), ...)

            return(p)
          }
)

# Pressure ----------------------------------------------------------------

#' Plot Mean Sea Level Atmospheric Pressure
#'
#' Plot the MSL atmospheric pressure through time.
#'
#' @importFrom ggplot2 ggplot aes_string geom_line ylab theme element_blank
#' theme_grey theme_bw theme_classic theme_gray theme_linedraw theme_light
#' theme_minimal facet_wrap
#'
#' @param x a cfPressure object.
#' @param y missing.
#' @param ggtheme character string (partially) matching the
#'                \code{\link[ggplot2]{ggtheme}} to be used for plotting, see
#'                'Theme Selection' below.
#' @param scales character string partially matching the \code{scales} argument
#'               in the \code{link[ggplot2]{facet_wrap}} function.
#' @param n_col the number of columns of plots (default 1).
#' @param ... further arguments passed to \code{\link[ggplot2]{theme}}.
#'
#' @aliases plot.cfPressure
#' @seealso \code{\link{plot,cfDataList,missing-method}} for general
#'   information on default plotting of \code{cfData} and \code{cfDataList}
#'   objects, and the links within. See \code{\link{cf_query}} for creating
#'   \code{cfPressure} objects.
#'
#'   Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
#'   to these methods.
#' @examples
#' \dontrun{
#' # Retrieve public hourly atmospheric pressure data for the last 30 days at
#' # Reefton Ews station
#'
#' # Subtract 30 days from today's date to get the start date
#' last_month = paste(as.character(Sys.Date() - 30), 0)
#'
#' reefton_pressure = cf_query(cf_user(), cf_datatype(7, 1, 1), cf_station(),
#'                             start_date = last_month)
#'
#' class(reefton_pressure) # cfPressure object
#'
#' # Plot the atmospheric pressure data using the defaults
#' plot(reefton_pressure)
#'
#' # Enlarge the text and add the observations as points
#' library(ggplot2) # for element_text() and geom_point()
#' plot(reefton_pressure, ggtheme = "bw", text = element_text(size = 16)) +
#'   geom_point(size = 3, shape = 1)
#'
#' # Save the plot as a png to the current working directory
#' library(ggplot2) # for ggsave()
#' ggsave("my_pressure_plot.png")
#' }
#' @export
setMethod("plot",
          signature(x = "cfPressure", y = "missing"),

          definition = function (x, y,
                                 ggtheme = c("grey", "gray", "bw", "linedraw",
                                             "light", "minimal", "classic"),
                                 scales = c("fixed", "free_x", "free_y", "free"),
                                 n_col = 1, ...){
            if (!is.numeric(n_col))
              stop("n_col must be a single number")

            ggtheme = match.arg(ggtheme)
            scales = match.arg(scales)

            x_df = as(x, "data.frame")
            names(x_df)[2:3] = c("date", "pressure")

            p = ggplot(x_df, aes_string(x = "date", y = "pressure")) +
              geom_line() +
              eval(call(paste0("theme_", ggtheme))) +
              ylab(x@data_label) +
              facet_wrap(~Station, scales = scales, ncol = n_col) +
              theme(axis.title.x = element_blank(),
                    ...)

            return(p)
          }
)

#' Default Clifro Plotting
#'
#' Plot \pkg{clifro} data based on the datatype.
#'
#' @param x a \code{cfData} or \code{cfDataList} object.
#' @param y missing for \code{cfData} objects, or a number representing the
#'          dataframe to plot if \code{x} is a \code{cfDataList} object.
#' @param ... arguments passed onto the different plotting methods.
#'
#' These methods are intended to simplify the data visualisation and exploration
#' of CliFlo data. The type of plot is determined by the type of the data output
#' from a \pkg{clifro} query. All of these methods plot individual plots for
#' each CliFlo station (if there is more than one in the query). If \code{x} is
#' a \code{cfDataList}, by default the first datatype will be plotted unless
#' \code{y} is supplied.
#'
#' The following table links the datatypes to the corresponding plot methods:
#'
#' \tabular{ll}{
#' \strong{Datatype} \tab \strong{Method}\cr
#' Wind \tab \code{\link{plot.cfWind}} for windrose, wind speed and direction
#' contour plots\cr
#' Rain \tab \code{\link{plot.cfRain}} for plotting rainfall (mm) through time\cr
#' Screen Obs \tab \code{\link{plot.cfScreenObs}} for timeseries plots of air,
#' wet bulb, and dew-point temperature plots\cr
#' Max/Min Temp \tab \code{\link{plot.cfTemp}} for maximum, minimum and
#' average temperature timeseries plots\cr
#' Earth Temp \tab \code{\link{plot.cfEarthTemp}} for earth temperature
#' timeseries plots\cr
#' Sunshine \tab \code{\link{plot.cfSunshine}} for accumulated, hourly or daily
#' sunshine, timeseries plots\cr
#' Pressure \tab \code{\link{plot.cfPressure}} for mean sea level atmospheric
#' pressure timeseries plots\cr
#' Other data \tab No default plot methods\cr
#' }
#'
#'
#' @seealso \code{\link{cf_query}} to retrieve the CliFlo data and create
#'  \code{cfData} objects.
#'
#'   Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
#'   to these methods.
#' @name plot.cfDataList
#' @rdname plot.cfDataList
#' @aliases plot,cfDataList,numeric-method
#' @importFrom methods setMethod
#' @importFrom graphics plot
setMethod("plot",
          signature(x = "cfDataList", y = "numeric"),

          definition = function (x, y, ...){

            y = y[1]

            if (y > length(x))
              stop("y needs to be a number between 1 and ", length(x))

            plot(x[y], ...)
          }
)

#' @rdname plot.cfDataList
#' @aliases plot,cfDataList,missing-method
#' @importFrom methods setMethod
#' @importFrom graphics plot
setMethod("plot",
          signature(x = "cfDataList", y = "missing"),

          definition = function (x, y, ...){
            plot(x[1], ...)
          }
)

#' @importFrom methods setMethod
#' @importFrom graphics plot
#' @rdname plot.cfDataList
#' @aliases plot,cfOther,missing-method
setMethod("plot",
          signature(x = "cfOther", y = "missing"),

          definition = function (x, y){
            stop("no default methods for this clifro datatype.")
          }
)
