#' Function to create an Abanico Plot.
#'
#' A plot is produced which allows comprehensive presentation of data precision
#' and its dispersion around a central value as well as illustration of a
#' kernel density estimate, histogram and/or dot plot of the dose values.
#'
#' The Abanico Plot is a combination of the classic Radial Plot
#' (\code{plot_RadialPlot}) and a kernel density estimate plot (e.g
#' \code{plot_KDE}). It allows straightforward visualisation of data precision,
#' error scatter around a user-defined central value and the combined
#' distribution of the values, on the actual scale of the measured data (e.g.
#' seconds, equivalent dose, years). The principle of the plot is shown in
#' Galbraith & Green (1990). The function authors are thankful for the
#' thoughtprovocing figure in this article. \cr The semi circle (z-axis) of the
#' classic Radial Plot is bent to a straight line here, which actually is the
#' basis for combining this polar (radial) part of the plot with any other
#' cartesian visualisation method (KDE, histogram, PDF and so on). Note that
#' the plot allows dispaying two measures of distribution. One is the 2-sigma
#' bar, which illustrates the spread in value errors, and the other is the
#' polygon, which stretches over both parts of the Abanico Plot (polar and
#' cartesian) and illustrates the actual spread in the values themselfes. \cr
#' Since the 2-sigma-bar is a polygon, it can be (and is) filled with shaded
#' lines. To change density (lines per inch, default is 15) and angle (default
#' is 45 degrees) of the shading lines, specify these parameters. See
#' \code{?polygon()} for further help. \cr The Abanico Plot supports other than
#' the weighted mean as measure of centrality. When it is obvious that the data
#' is not (log-)normally distributed, the mean (weighted or not) cannot be a
#' valid measure of centrality and hence central dose. Accordingly, the median
#' and the weighted median can be chosen as well to represent a proper measure
#' of centrality (e.g. \code{centrality = "median.weighted"}). Also
#' user-defined numeric values (e.g. from the central age model) can be used if
#' this appears appropriate. \cr The proportion of the polar part and the
#' cartesian part of the Abanico Plot can be modfied for display reasons
#' (\code{plot.ratio = 0.75}). By default, the polar part spreads over 75 \%
#' and leaves 25 \% for the part that shows the KDE graph.\cr\cr
#' A statistic summary, i.e. a collection of statistic measures of
#' centrality and dispersion (and further measures) can be added by specifying
#' one or more of the following keywords: \code{"n"} (number of samples),
#' \code{"mean"} (mean De value), \code{"mean.weighted"} (error-weighted mean),
#' \code{"median"} (median of the De values), \code{"sdrel"} (relative standard
#' deviation in percent), \code{"sdrel.weighted"} (error-weighted relative
#' standard deviation in percent), \code{"sdabs"} (absolute standard deviation),
#' \code{"sdabs.weighted"} (error-weighted absolute standard deviation),
#' \code{"serel"} (relative standard error), \code{"serel.weighted"} (
#' error-weighted relative standard error), \code{"seabs"} (absolute standard
#' error), \code{"seabs.weighted"} (error-weighted absolute standard error),
#' \code{"in.2s"} (percent of samples in 2-sigma range),
#' \code{"kurtosis"} (kurtosis) and \code{"skewness"} (skewness). \cr\cr
#'
#' The optional parameter \code{layout} allows to modify the entire plot more
#' sophisticated. Each element of the plot can be addressed and its properties
#' can be defined. This includes font type, size and decoration, colours and
#' sizes of all plot items. To infer the definition of a specific layout style
#' cf. \code{get_Layout()} or type eg. for the layout type \code{"journal"}
#' \code{get_Layout("journal")}. A layout type can be modified by the user by
#' assigning new values to the list object.\cr\cr It is possible for the
#' z-scale to specify where ticks are to be drawn by using the parameter
#' \code{at}, e.g. \code{at = seq(80, 200, 20)}, cf. function documentation of
#' \code{axis}. Specifying tick positions manually overrides a
#' \code{zlim}-definition.
#'
#' @param data \code{\link{data.frame}} or \code{\linkS4class{RLum.Results}}
#' object (required): for \code{data.frame} two columns: De (\code{data[,1]})
#' and De error (\code{data[,2]}). To plot several data sets in one plot the
#' data sets must be provided as \code{list}, e.g. \code{list(data.1, data.2)}.
#'
#' @param na.rm \code{\link{logical}} (with default): exclude NA values
#' from the data set prior to any further operations.
#'
#' @param log.z \code{\link{logical}} (with default): Option to display the
#' z-axis in logarithmic scale. Default is \code{TRUE}.
#'
#' @param z.0 \code{\link{character}} or \code{\link{numeric}}: User-defined
#' central value, used for centering of data. One out of \code{"mean"},
#' \code{"mean.weighted"} and \code{"median"} or a numeric value (not its
#' logarithm). Default is \code{"mean.weighted"}.
#'
#' @param dispersion \code{\link{character}} (with default): measure of
#' dispersion, used for drawing the scatter polygon. One out of \code{"qr"}
#' (quartile range), \code{"pnn"} (symmetric percentile range with nn the lower
#' percentile, e.g. \code{"p05"} depicting the range between 5 and 95 %),
#' \code{"sd"} (standard deviation) and \code{"2sd"} (2 standard deviations),
#' default is \code{"qr"}. Note that \code{"sd"} and \code{"2sd"} are only
#' meaningful in combination with \code{"z.0 = 'mean'"} because the unweighted
#' mean is used to center the polygon.
#'
#' @param plot.ratio \code{\link{numeric}}: Relative space, given to the radial
#' versus the cartesian plot part, deault is \code{0.75}.
#'
#' @param rotate \code{\link{logical}}: Option to turn the plot by 90 degrees.
#'
#' @param mtext \code{\link{character}}: additional text below the plot title.
#'
#' @param summary \code{\link{character}} (optional): add statistic measures of
#' centrality and dispersion to the plot. Can be one or more of several
#' keywords. See details for available keywords.
#'
#' @param summary.pos \code{\link{numeric}} or \code{\link{character}} (with
#' default): optional position coordinates or keyword (e.g. \code{"topright"})
#' for the statistical summary. Alternatively, the keyword \code{"sub"} may be
#' specified to place the summary below the plot header. However, this latter
#' option in only possible if \code{mtext} is not used.
#'
#' @param legend \code{\link{character}} vector (optional): legend content to
#' be added to the plot.
#'
#' @param legend.pos \code{\link{numeric}} or \code{\link{character}} (with
#' default): optional position coordinates or keyword (e.g. \code{"topright"})
#' for the legend to be plotted.
#'
#' @param stats \code{\link{character}}: additional labels of statistically
#' important values in the plot. One or more out of the following:
#' \code{"min"}, \code{"max"}, \code{"median"}.
#'
#' @param rug \code{\link{logical}}: Option to add a rug to the KDE part, to
#' indicate the location of individual values.
#'
#' @param kde \code{\link{logical}}: Option to add a KDE plot to the dispersion
#' part, default is \code{TRUE}.
#'
#' @param hist \code{\link{logical}}: Option to add a histogram to the
#' dispersion part. Only meaningful when not more than one data set is plotted.
#'
#' @param dots \code{\link{logical}}: Option to add a dot plot to the
#' dispersion part. If number of dots exceeds space in the dispersion part, a
#' square indicates this.
#'
#' @param y.axis \code{\link{logical}}: Option to hide y-axis labels. Useful
#' for data with small scatter.
#'
#' @param error.bars \code{\link{logical}}: Option to show De-errors as error
#' bars on De-points. Useful in combination with \code{y.axis = FALSE, bar.col
#' = "none"}.
#'
#' @param bar \code{\link{numeric}} (with default): option to add one or more
#' dispersion bars (i.e., bar showing the 2-sigma range) centered at the
#' defined values. By default a bar is drawn according to \code{"z.0"}. To omit
#' the bar set \code{"bar = FALSE"}.
#'
#' @param bar.col \code{\link{character}} or \code{\link{numeric}} (with
#' default): colour of the dispersion bar. Default is \code{"grey60"}.
#'
#' @param polygon.col \code{\link{character}} or \code{\link{numeric}} (with
#' default): colour of the polygon showing the data scatter. Sometimes this
#' polygon may be omitted for clarity. To disable it use \code{FALSE} or
#' \code{polygon = FALSE}. Default is \code{"grey80"}.
#'
#' @param line \code{\link{numeric}}: numeric values of the additional lines to
#' be added.
#'
#' @param line.col \code{\link{character}} or \code{\link{numeric}}: colour of
#' the additional lines.
#'
#' @param line.label \code{\link{character}}: labels for the additional lines.
#'
#' @param grid.col \code{\link{character}} or \code{\link{numeric}} (with
#' default): colour of the grid lines (originating at [0,0] and strechting to
#' the z-scale). To disable grid lines use \code{FALSE}. Default is
#' \code{"grey"}.
#'
#' @param frame \code{\link{numeric}} (with default): option to modify the
#' plot frame type. Can be one out of \code{0} (no frame), \code{1} (frame
#' originates at 0,0 and runs along min/max isochrons), \code{2} (frame
#' embraces the 2-sigma bar), \code{3} (frame embraces the entire plot as a
#' rectangle).Default is \code{1}.
#'
#' @param bw \code{\link{character}} (with default): bin-width for KDE, choose
#' a numeric value for manual setting.
#'
#' @param output \code{\link{logical}}: Optional output of numerical plot
#' parameters. These can be useful to reproduce similar plots. Default is
#' \code{FALSE}.
#'
#' @param \dots Further plot arguments to pass. \code{xlab} must be a vector of
#' length 2, specifying the upper and lower x-axes labels.
#'
#' @return returns a plot object and, optionally, a list with plot calculus
#' data.
#'
#' @section Function version: 0.1.7
#'
#' @author Michael Dietze, GFZ Potsdam (Germany),\cr Sebastian Kreutzer,
#' IRAMAT-CRP2A, Universite Bordeaux Montaigne (France)\cr Inspired by a plot
#' introduced by Galbraith & Green (1990)
#'
#' @seealso \code{\link{plot_RadialPlot}}, \code{\link{plot_KDE}},
#' \code{\link{plot_Histogram}}
#'
#' @references Galbraith, R. & Green, P., 1990. Estimating the component ages
#' in a finite mixture. International Journal of Radiation Applications and
#' Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 17 (3),
#' 197-206.
#'
#' Dietze, M., Kreutzer, S., Burow, C., Fuchs, M.C., Fischer, M., Schmidt, C., 2015.
#' The abanico plot: visualising chronometric data with individual standard errors.
#' Quaternary Geochronology. doi:10.1016/j.quageo.2015.09.003
#'
#' @examples
#'
#' ## load example data and recalculate to Gray
#' data(ExampleData.DeValues, envir = environment())
#' ExampleData.DeValues <- ExampleData.DeValues$CA1
#'
#' ## plot the example data straightforward
#' plot_AbanicoPlot(data = ExampleData.DeValues)
#'
#' ## now with linear z-scale
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  log.z = FALSE)
#'
#' ## now with output of the plot parameters
#' plot1 <- plot_AbanicoPlot(data = ExampleData.DeValues,
#'                           output = TRUE)
#' str(plot1)
#' plot1$zlim
#'
#' ## now with adjusted z-scale limits
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  zlim = c(10, 200))
#'
#' ## now with adjusted x-scale limits
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  xlim = c(0, 20))
#'
#' ## now with rug to indicate individual values in KDE part
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  rug = TRUE)
#'
#' ## now with a smaller bandwidth for the KDE plot
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  bw = 0.04)
#'
#' ## now with a histogram instead of the KDE plot
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  hist = TRUE,
#'                  kde = FALSE)
#'
#' ## now with a KDE plot and histogram with manual number of bins
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  hist = TRUE,
#'                  breaks = 20)
#'
#' ## now with a KDE plot and a dot plot
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  dots = TRUE)
#'
#' ## now with user-defined plot ratio
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  plot.ratio = 0.5)

#' ## now with user-defined central value
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  z.0 = 70)
#'
#' ## now with median as central value
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  z.0 = "median")
#'
#' ## now with the 17-83 percentile range as definition of scatter
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  z.0 = "median",
#'                  dispersion = "p17")
#'
#' ## now with user-defined green line for minimum age model
#' CAM <- calc_CentralDose(ExampleData.DeValues,
#'                         plot = FALSE)
#'
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  line = CAM,
#'                  line.col = "darkgreen",
#'                  line.label = "CAM")
#'
#' ## now create plot with legend, colour, different points and smaller scale
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  legend = "Sample 1",
#'                  col = "tomato4",
#'                  bar.col = "peachpuff",
#'                  pch = "R",
#'                  cex = 0.8)
#'
#' ## now without 2-sigma bar, polygon, grid lines and central value line
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  bar.col = FALSE,
#'                  polygon.col = FALSE,
#'                  grid.col = FALSE,
#'                  y.axis = FALSE,
#'                  lwd = 0)
#'
#' ## now with direct display of De errors, without 2-sigma bar
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  bar.col = FALSE,
#'                  ylab = "",
#'                  y.axis = FALSE,
#'                  error.bars = TRUE)
#'
#' ## now with user-defined axes labels
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  xlab = c("Data error (%)",
#'                           "Data precision"),
#'                  ylab = "Scatter",
#'                  zlab = "Equivalent dose (Gy)")
#'
#' ## now with minimum, maximum and median value indicated
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  stats = c("min", "max", "median"))
#'
#' ## now with a brief statistical summary as subheader
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  summary = c("n", "in.2s"))
#'
#' ## now with another statistical summary
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  summary = c("mean.weighted", "median"),
#'                  summary.pos = "topleft")
#'
#' ## now a plot with two 2-sigma bars for one data set
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  bar = c(30, 100))
#'
#' ## now the data set is split into sub-groups, one is manipulated
#' data.1 <- ExampleData.DeValues[1:30,]
#' data.2 <- ExampleData.DeValues[31:62,] * 1.3
#'
#' ## now a common dataset is created from the two subgroups
#' data.3 <- list(data.1, data.2)
#'
#' ## now the two data sets are plotted in one plot
#' plot_AbanicoPlot(data = data.3)
#'
#' ## now with some graphical modification
#' plot_AbanicoPlot(data = data.3,
#'                  z.0 = "median",
#'                  col = c("steelblue4", "orange4"),
#'                  bar.col = c("steelblue3", "orange3"),
#'                  polygon.col = c("steelblue1", "orange1"),
#'                  pch = c(2, 6),
#'                  angle = c(30, 50),
#'                  summary = c("n", "in.2s", "median"))
#'
#' ## create Abanico plot with predefined layout definition
#' plot_AbanicoPlot(data = ExampleData.DeValues,
#'                  layout = "journal")
#'
#' ## now with predefined layout definition and further modifications
#' plot_AbanicoPlot(data = data.3,
#'                  z.0 = "median",
#'                  layout = "journal",
#'                  col = c("steelblue4", "orange4"),
#'                  bar.col = adjustcolor(c("steelblue3", "orange3"),
#'                                        alpha.f = 0.5),
#'                  polygon.col = c("steelblue3", "orange3"))
#'
#' ## for further information on layout definitions see documentation
#' ## of function get_Layout()
#'
#' @export
plot_AbanicoPlot <- function(
  data,
  na.rm = TRUE,
  log.z = TRUE,
  z.0 = "mean.weighted",
  dispersion = "qr",
  plot.ratio = 0.75,
  rotate = FALSE,
  mtext,
  summary,
  summary.pos,
  legend,
  legend.pos,
  stats,
  rug = FALSE,
  kde = TRUE,
  hist = FALSE,
  dots = FALSE,
  y.axis = TRUE,
  error.bars = FALSE,
  bar,
  bar.col,
  polygon.col,
  line,
  line.col,
  line.label,
  grid.col,
  frame = 1,
  bw = "SJ",
  output = FALSE,
  ...
) {
  ## check data and parameter consistency--------------------------------------

  ## Homogenise input data format
  if(is(data, "list") == FALSE) {
    data <- list(data)
    }

  ## Check input data
  for(i in 1:length(data)) {
    if(is(data[[i]], "RLum.Results") == FALSE &
         is(data[[i]], "data.frame") == FALSE) {
      stop(paste("[plot_AbanicoPlot()] Input data format is neither",
                 "'data.frame' nor 'RLum.Results'"))
    } else {
      if(is(data[[i]], "RLum.Results") == TRUE) {
        data[[i]] <- get_RLum(data[[i]])[,c(1:2)]
      }
    }
  }

  ## Check input data
  for(i in 1:length(data)) {
    if(is(data[[i]], "RLum.Results") == FALSE &
         is(data[[i]], "data.frame") == FALSE) {
      stop(paste("[plot_AbanicoPlot()] Input data format is neither",
                 "'data.frame' nor 'RLum.Results'"))
    } else {
      if(is(data[[i]], "RLum.Results") == TRUE) {
        data[[i]] <- get_RLum(data[[i]])[,c(1:2)]
      }
    }
  }

  ## optionally, remove NA-values
  if(na.rm == TRUE) {
    for(i in 1:length(data)) {

      n.NA <- sum(!complete.cases(data[[i]]))

      if(n.NA == 1) {message(paste0("[plot_AbanicoPlot()] data set (", i, "): 1 NA value excluded."))
      } else if(n.NA > 1) {
        message(paste0("[plot_AbanicoPlot()] data set (", i,"): ",n.NA, " NA values excluded."))
      }

      data[[i]] <- na.exclude(data[[i]])
    }
  }

  ##AFTER NA removal, we should check the data set carefully again ...
  ##(1)
  ##check if there is still data left in the entire set
  if(all(sapply(data, nrow) == 0)){

    warning("[plot_AbanicoPlot()] Nothing plotted, your data set is empty!")
    return(NULL)

  }
  ##(2)
  ##check for sets with only 1 row or 0 rows at all
  else if(any(sapply(data, nrow) <= 1)){

    ##select problematic sets and remove the entries from the list
    NArm.id <- which(sapply(data, nrow) <= 1)
    data[NArm.id] <- NULL

    warning(paste0("[plot_AbanicoPlot()] Data sets ",paste(NArm.id, collapse = ", ")," are found to be empty or consisting of only 1 row. Sets removed!"))

    rm(NArm.id)

    ##unfortunately, the data set might become now empty at all
    if(length(data) == 0){
      warning("[plot_AbanicoPlot()] After removing invalid entries, nothing is plotted!")
      return(NULL)

    }

  }

  ## check for zero-error values
  for(i in 1:length(data)) {

    if(length(data[[i]]) < 2) {
      stop("Data without errors cannot be displayed!")
    }

    if(sum(data[[i]][,2] == 0) > 0) {
      data[[i]] <- data[[i]][data[[i]][,2] > 0,]

      if(nrow(data[[i]]) < 1) {
        stop("Data set contains only values with zero errors.")
      }

      warning("Values with zero errors cannot be displayed and were removed!")
    }
  }

  ## save original plot parameters and restore them when the function ends or stops
  par.old.full <- par(no.readonly = TRUE)

  ##this ensures, that par() for several plots on one page is  respected ...
  if(sum(par()$mfrow) == 2 & sum(par()$mfcol) == 2){
    on.exit(par(par.old.full))
  }

  ## check/set layout definitions
  if("layout" %in% names(list(...))) {
    layout = get_Layout(layout = list(...)$layout)
  } else {
    layout <- get_Layout(layout = "default")
  }

  if(missing(stats) == TRUE) {
    stats <- numeric(0)
  }

  if(missing(bar) == TRUE) {
    bar <- rep(TRUE, length(data))
  }

  if(missing(bar.col) == TRUE) {
    bar.fill <- rep(x = rep(x = layout$abanico$colour$bar.fill,
                    length.out = length(data)), length(bar))
    bar.line <- rep(rep(layout$abanico$colour$bar.line,
                    length.out = length(data)), length(bar))
  } else {
    bar.fill <- bar.col
    bar.line <- NA
  }

    if(missing(polygon.col) == TRUE) {
      polygon.fill <- rep(layout$abanico$colour$poly.fill,
                          length.out = length(data))
      polygon.line <- rep(layout$abanico$colour$poly.line,
                          length.out = length(data))
    } else {
    polygon.fill <- polygon.col
    polygon.line <- NA
  }

  if(missing(grid.col) == TRUE) {
    grid.major <- layout$abanico$colour$grid.major
    grid.minor <- layout$abanico$colour$grid.minor
  } else {
    if(length(grid.col) == 1) {
      grid.major <- grid.col[1]
      grid.minor <- grid.col[1]
    } else {
      grid.major <- grid.col[1]
      grid.minor <- grid.col[2]
    }
  }

  if(missing(summary) == TRUE) {
    summary <- c("n", "in.ci")
  }

  if(missing(summary.pos) == TRUE) {
    summary.pos <- "sub"
  }

  if(missing(mtext) == TRUE) {
    mtext <- ""
  }

  ## create preliminary global data set
  De.global <- data[[1]][,1]
  if(length(data) > 1) {
    for(i in 2:length(data)) {
      De.global <- c(De.global, data[[i]][,1])
    }
  }

  ## calculate major preliminary tick values and tick difference
  extraArgs <- list(...)
  if("zlim" %in% names(extraArgs)) {
    limits.z <- extraArgs$zlim
  } else {
    z.span <- (mean(De.global) * 0.5) / (sd(De.global) * 100)
    z.span <- ifelse(z.span > 1, 0.9, z.span)
    limits.z <- c((ifelse(min(De.global) <= 0, 1.1, 0.9) - z.span) *
                    min(De.global),
                  (1.1 + z.span) * max(De.global))
  }

  if("at" %in% names(extraArgs)) {
    ticks <- extraArgs$at
  } else {
    ticks <- round(pretty(limits.z, n = 5), 3)
  }

  if("breaks" %in% names(extraArgs)) {
    breaks <- extraArgs$breaks
  } else {
    breaks <- "Sturges"
  }

  ## check/set bw-parameter
  for(i in 1:length(data)) {
    bw.test <- try(density(x = data[[i]][,1],
                           bw = bw),
                   silent = TRUE)
    if(grepl(pattern = "Error", x = bw.test[1]) == TRUE) {
      bw <- "nrd0"
      warning("Option for bw not possible. Set to nrd0!")
    }
  }

  if ("fun" %in% names(extraArgs)) {
    fun <- list(...)$fun

  } else {
    fun <- FALSE
  }

  ## check for negative values, stoppp function, but do not stop
  if(min(De.global) <= 0) {
    message("\n [plot_AbanicoPlot()] Data contains negative or zero values. Nothing plotted!")
    return(NULL)
  }

  ## calculate and append statistical measures --------------------------------

  ## z-values based on log-option
  z <- lapply(1:length(data), function(x){
    if(log.z == TRUE) {
      log(data[[x]][,1])
      } else {
        data[[x]][,1]
        }
    })
  if(is(z, "list") == FALSE) {
    z <- list(z)
    }
  data <- lapply(1:length(data), function(x) {
    cbind(data[[x]], z[[x]])
    })
  rm(z)

  ## calculate dispersion based on log-option
  se <- lapply(1:length(data), function(x){
    if(log.z == TRUE) {
      data[[x]][,2] / data[[x]][,1]
      } else {
        data[[x]][,2]
        }
    })
  if(is(se, "list") == FALSE) {
    se <- list(se)
    }
  data <- lapply(1:length(data), function(x) {
    cbind(data[[x]], se[[x]])
    })
  rm(se)

  ## calculate initial data statistics
  stats.init <- list(NA)
  for(i in 1:length(data)) {
    stats.init[[length(stats.init) + 1]] <- calc_Statistics(data = data[[i]][,3:4])
  }
  stats.init[[1]] <- NULL

  ## calculate central values
  if(z.0 == "mean") {
    z.central <- lapply(1:length(data), function(x){
      rep(stats.init[[x]]$unweighted$mean,
          length(data[[x]][,3]))})
  } else if(z.0 == "median") {
    z.central <- lapply(1:length(data), function(x){
      rep(stats.init[[x]]$unweighted$median,
          length(data[[x]][,3]))})
  } else  if(z.0 == "mean.weighted") {
    z.central <- lapply(1:length(data), function(x){
      rep(stats.init[[x]]$weighted$mean,
          length(data[[x]][,3]))})
  } else if(is.numeric(z.0) == TRUE) {
    z.central <- lapply(1:length(data), function(x){
      rep(ifelse(log.z == TRUE,
                 log(z.0),
                 z.0),
          length(data[[x]][,3]))})
  } else {
    stop("Value for z.0 not supported!")
  }

  data <- lapply(1:length(data), function(x) {
    cbind(data[[x]], z.central[[x]])})
  rm(z.central)

  ## calculate precision
  precision <- lapply(1:length(data), function(x){
    1 / data[[x]][,4]})
  if(is(precision, "list") == FALSE) {precision <- list(precision)}
  data <- lapply(1:length(data), function(x) {
    cbind(data[[x]], precision[[x]])})
  rm(precision)

  ## calculate standardised estimate
  std.estimate <- lapply(1:length(data), function(x){
    (data[[x]][,3] - data[[x]][,5]) / data[[x]][,4]})
  if(is(std.estimate, "list") == FALSE) {std.estimate <- list(std.estimate)}
  data <- lapply(1:length(data), function(x) {
    cbind(data[[x]], std.estimate[[x]])})

  ## append empty standard estimate for plotting
  data <- lapply(1:length(data), function(x) {
    cbind(data[[x]], std.estimate[[x]])})
  rm(std.estimate)

  ## append optional weights for KDE curve
  if("weights" %in% names(extraArgs)) {
    if(extraArgs$weights == TRUE) {
      wgt <- lapply(1:length(data), function(x){
        (1 / data[[x]][,2]) / sum(1 / data[[x]][,2]^2)
      })

      if(is(wgt, "list") == FALSE) {
        wgt <- list(wgt)
      }

      data <- lapply(1:length(data), function(x) {
        cbind(data[[x]], wgt[[x]])})

      rm(wgt)
    } else {
      wgt <- lapply(1:length(data), function(x){
        rep(x = 1, times = nrow(data[[x]])) /
          sum(rep(x = 1, times = nrow(data[[x]])))
      })

      if(is(wgt, "list") == FALSE) {
        wgt <- list(wgt)
      }

      data <- lapply(1:length(data), function(x) {
        cbind(data[[x]], wgt[[x]])})

      rm(wgt)
    }
  } else {
    wgt <- lapply(1:length(data), function(x){
      rep(x = 1, times = nrow(data[[x]])) /
        sum(rep(x = 1, times = nrow(data[[x]])))
    })

    if(is(wgt, "list") == FALSE) {
      wgt <- list(wgt)
    }

    data <- lapply(1:length(data), function(x) {
      cbind(data[[x]], wgt[[x]])})

    rm(wgt)
  }

  ## generate global data set
  data.global <- cbind(data[[1]],
                       rep(x = 1, times = nrow(data[[1]])))
  colnames(data.global) <- rep("", 10)

  if(length(data) > 1) {
    for(i in 2:length(data)) {
      data.add <- cbind(data[[i]],
                        rep(x = i, times = nrow(data[[i]])))
      colnames(data.add) <- rep("", 10)
      data.global <- rbind(data.global,
                           data.add)
    }
  }

  ## create column names
  colnames(data.global) <- c("De",
                             "error",
                             "z",
                             "se",
                             "z.central",
                             "precision",
                             "std.estimate",
                             "std.estimate.plot",
                             "weights",
                             "data set")

  ## calculate global data statistics
  stats.global <- calc_Statistics(data = data.global[,3:4])

  ## calculate global central value
  if(z.0 == "mean") {
    z.central.global <- stats.global$unweighted$mean
  } else if(z.0 == "median") {
    z.central.global <- stats.global$unweighted$median
  } else  if(z.0 == "mean.weighted") {
    z.central.global <- stats.global$weighted$mean
  } else if(is.numeric(z.0) == TRUE) {
    z.central.global <- ifelse(log.z == TRUE,
                               log(z.0),
                               z.0)
  } else {
    stop("Value for z.0 not supported!")
  }

  ## create column names
  for(i in 1:length(data)) {
    colnames(data[[i]]) <- c("De",
                             "error",
                             "z",
                             "se",
                             "z.central",
                             "precision",
                             "std.estimate",
                             "std.estimate.plot",
                             "weights")
  }

  ## re-calculate standardised estimate for plotting
  for(i in 1:length(data)) {
    data[[i]][,8] <- (data[[i]][,3] - z.central.global) / data[[i]][,4]
  }

  data.global.plot <- data[[1]][,8]
  if(length(data) > 1) {
    for(i in 2:length(data)) {
      data.global.plot <- c(data.global.plot, data[[i]][,8])
    }
  }
  data.global[,8] <- data.global.plot

  ## print message for too small scatter
  if(max(abs(1 / data.global[6])) < 0.02) {
    small.sigma <- TRUE
    message("[plot_AbanicoPlot()] Attention, small standardised estimate scatter. Toggle off y.axis?")
  }

  ## read out additional arguments---------------------------------------------
  extraArgs <- list(...)

  main <- if("main" %in% names(extraArgs)) {
    extraArgs$main
  } else {
      expression(paste(D[e], " distribution"))
    }

  sub <- if("sub" %in% names(extraArgs)) {
    extraArgs$sub
  } else {
      ""
    }

  if("xlab" %in% names(extraArgs)) {
    if(length(extraArgs$xlab) != 2) {
      if (length(extraArgs$xlab) == 3) {
        xlab <- c(extraArgs$xlab[1:2], "Density")
      } else {
        stop("Argmuent xlab is not of length 2!")
      }
    } else {xlab <- c(extraArgs$xlab, "Density")}
  } else {
    xlab <- c(if(log.z == TRUE) {
      "Relative standard error (%)"
    } else {
      "Standard error"
    },
    "Precision",
    "Density")
  }

  ylab <- if("ylab" %in% names(extraArgs)) {
    extraArgs$ylab
  } else {
    "Standardised estimate"
  }

  zlab <- if("zlab" %in% names(extraArgs)) {
    extraArgs$zlab
  } else {
    expression(paste(D[e], " (Gy)"))
  }

  if("zlim" %in% names(extraArgs)) {
    limits.z <- extraArgs$zlim
  } else {
    z.span <- (mean(data.global[,1]) * 0.5) / (sd(data.global[,1]) * 100)
    z.span <- ifelse(z.span > 1, 0.9, z.span)
    limits.z <- c((0.9 - z.span) * min(data.global[[1]]),
                  (1.1 + z.span) * max(data.global[[1]]))
  }

  if("xlim" %in% names(extraArgs)) {
    limits.x <- extraArgs$xlim
  } else {
    limits.x <- c(0, max(data.global[,6]) * 1.05)
  }

  if(limits.x[1] != 0) {
    limits.x[1] <- 0
    warning("Lower x-axis limit not set to zero, issue corrected!")
  }

  if("ylim" %in% names(extraArgs)) {
    limits.y <- extraArgs$ylim
  } else {
    y.span <- (mean(data.global[,1]) * 10) / (sd(data.global[,1]) * 100)
    y.span <- ifelse(y.span > 1, 0.98, y.span)
    limits.y <- c(-(1 + y.span) * max(abs(data.global[,7])),
                  (1 + y.span) * max(abs(data.global[,7])))
  }

  cex <- if("cex" %in% names(extraArgs)) {
    extraArgs$cex
  } else {
    1
  }

  lty <- if("lty" %in% names(extraArgs)) {
    extraArgs$lty
  } else {
    rep(rep(2, length(data)), length(bar))
  }

  lwd <- if("lwd" %in% names(extraArgs)) {
    extraArgs$lwd
  } else {
    rep(rep(1, length(data)), length(bar))
  }

  pch <- if("pch" %in% names(extraArgs)) {
    extraArgs$pch
  } else {
    rep(20, length(data))
  }

  if("col" %in% names(extraArgs)) {
    bar.col <- extraArgs$col
    kde.line <- extraArgs$col
    kde.fill <- NA
    value.dot <- extraArgs$col
    value.bar <- extraArgs$col
    value.rug <- extraArgs$col
    summary.col <- extraArgs$col
    centrality.col <- extraArgs$col
  } else {
    if(length(layout$abanico$colour$bar) == 1) {
      bar.col <- 1:length(data)
    } else {
      bar.col <- layout$abanico$colour$bar.col
    }

    if(length(layout$abanico$colour$kde.line) == 1) {
      kde.line <- 1:length(data)
    } else {
      kde.line <- layout$abanico$colour$kde.line
    }

    if(length(layout$abanico$colour$kde.fill) == 1) {
      kde.fill <- rep(layout$abanico$colour$kde.fill, length(data))
    } else {
      kde.fill <- layout$abanico$colour$kde.fill
    }

    if(length(layout$abanico$colour$value.dot) == 1) {
      value.dot <- 1:length(data)
    } else {
      value.dot <- layout$abanico$colour$value.dot
    }

    if(length(layout$abanico$colour$value.bar) == 1) {
      value.bar <- 1:length(data)
    } else {
      value.bar <- layout$abanico$colour$value.bar
    }

    if(length(layout$abanico$colour$value.rug) == 1) {
      value.rug <- 1:length(data)
    } else {
      value.rug <- layout$abanico$colour$value.rug
    }

    if(length(layout$abanico$colour$summary) == 1) {
      summary.col <- 1:length(data)
    } else {
      summary.col <- layout$abanico$colour$summary
    }

    if(length(layout$abanico$colour$centrality) == 1) {
      centrality.col <- rep(x = 1:length(data), times = length(bar))
    } else {
      centrality.col <- rep(x = layout$abanico$colour$centrality,
                            times = length(bar))
    }
  }

  ## update central line colour
  centrality.col <- rep(centrality.col, length(bar))

  tck <- if("tck" %in% names(extraArgs)) {
    extraArgs$tck
  } else {
    NA
  }

  tcl <- if("tcl" %in% names(extraArgs)) {
    extraArgs$tcl
  } else {
    -0.5
  }

  ## define auxiliary plot parameters -----------------------------------------

  ## create empty plot to update plot parameters
  if(rotate == FALSE) {
    plot(NA,
         xlim = c(limits.x[1], limits.x[2] * (1 / plot.ratio)),
         ylim = limits.y,
         main = "",
         sub = "",
         xlab = "",
         ylab = "",
         xaxs = "i",
         yaxs = "i",
         frame.plot = FALSE,
         axes = FALSE)
  } else {
    plot(NA,
         xlim = limits.y,
         ylim = c(limits.x[1], limits.x[2] * (1 / plot.ratio)),
         main = "",
         sub = "",
         xlab = "",
         ylab = "",
         xaxs = "i",
         yaxs = "i",
         frame.plot = FALSE,
         axes = FALSE)
  }

  ## calculate conversion factor for plot coordinates
  f <- 0

  ## calculate major and minor z-tick values
  if("at" %in% names(extraArgs)) {
    tick.values.major <- extraArgs$at
    tick.values.minor <- extraArgs$at
  } else {
    tick.values.major <- signif(pretty(limits.z, n = 5), 3)
    tick.values.minor <- signif(pretty(limits.z, n = 25), 3)
  }

  tick.values.major <- tick.values.major[tick.values.major >=
                                           min(tick.values.minor)]
  tick.values.major <- tick.values.major[tick.values.major <=
                                           max(tick.values.minor)]
  tick.values.major <- tick.values.major[tick.values.major >=
                                           limits.z[1]]
  tick.values.major <- tick.values.major[tick.values.major <=
                                           limits.z[2]]
  tick.values.minor <- tick.values.minor[tick.values.minor >=
                                           limits.z[1]]
  tick.values.minor <- tick.values.minor[tick.values.minor <=
                                           limits.z[2]]


  if(log.z == TRUE) {

    tick.values.major[which(tick.values.major==0)] <- 1
    tick.values.minor[which(tick.values.minor==0)] <- 1

    tick.values.major <- log(tick.values.major)
    tick.values.minor <- log(tick.values.minor)
  }

  ## calculate z-axis radius
  r <- max(sqrt((limits.x[2])^2 + (data.global[,7] * f)^2))

  ## create z-axes labels
  if(log.z == TRUE) {
    label.z.text <- signif(exp(tick.values.major), 3)
  } else {
    label.z.text <- signif(tick.values.major, 3)
  }

  ## calculate node coordinates for semi-circle
  ellipse.values <- c(min(ifelse(log.z == TRUE,
                                 log(limits.z[1]),
                                 limits.z[1]),
                          tick.values.major,
                          tick.values.minor),
                      max(ifelse(log.z == TRUE,
                                 log(limits.z[2]),
                                 limits.z[2]),
                          tick.values.major,
                          tick.values.minor))

  ## correct for unpleasant value
  ellipse.values[ellipse.values == -Inf] <- 0

  if(rotate == FALSE) {
    ellipse.x <- r / sqrt(1 + f^2 * (ellipse.values - z.central.global)^2)
    ellipse.y <- (ellipse.values - z.central.global) * ellipse.x
  } else {
    ellipse.y <- r / sqrt(1 + f^2 * (ellipse.values - z.central.global)^2)
    ellipse.x <- (ellipse.values - z.central.global) * ellipse.y
  }

  ellipse <- cbind(ellipse.x, ellipse.y)


  ## calculate statistical labels
  if(length(stats == 1)) {stats <- rep(stats, 2)}
  stats.data <- matrix(nrow = 3, ncol = 3)
  data.stats <- as.numeric(data.global[,1])

  if("min" %in% stats == TRUE) {
    stats.data[1, 3] <- data.stats[data.stats == min(data.stats)][1]
    stats.data[1, 1] <- data.global[data.stats == stats.data[1, 3], 6][1]
    stats.data[1, 2] <- data.global[data.stats == stats.data[1, 3], 8][1]
  }

  if("max" %in% stats == TRUE) {
    stats.data[2, 3] <- data.stats[data.stats == max(data.stats)][1]
    stats.data[2, 1] <- data.global[data.stats == stats.data[2, 3], 6][1]
    stats.data[2, 2] <- data.global[data.stats == stats.data[2, 3], 8][1]
  }

  if("median" %in% stats == TRUE) {
    stats.data[3, 3] <- data.stats[data.stats == quantile(data.stats, 0.5, type = 3)]
    stats.data[3, 1] <- data.global[data.stats == stats.data[3, 3], 6][1]
    stats.data[3, 2] <- data.global[data.stats == stats.data[3, 3], 8][1]
  }

  ## re-calculate axes limits if necessary
  if(rotate == FALSE) {
    limits.z.x <- range(ellipse[,1])
    limits.z.y <- range(ellipse[,2])
  } else {
    limits.z.x <- range(ellipse[,2])
    limits.z.y <- range(ellipse[,1])
  }

  if(!("ylim" %in% names(extraArgs))) {
    if(limits.z.y[1] < 0.66 * limits.y[1]) {
      limits.y[1] <- 1.8 * limits.z.y[1]
    }
    if(limits.z.y[2] > 0.77 * limits.y[2]) {
      limits.y[2] <- 1.3 * limits.z.y[2]
    }

    if(rotate == TRUE) {
      limits.y <- c(-max(abs(limits.y)), max(abs(limits.y)))
    }

  }
  if(!("xlim" %in% names(extraArgs))) {
    if(limits.z.x[2] > 1.1 * limits.x[2]) {
      limits.x[2] <- limits.z.x[2]
    }
  }

  ## calculate and paste statistical summary
  De.stats <- matrix(nrow = length(data), ncol = 18)
  colnames(De.stats) <- c("n",
                          "mean",
                          "mean.weighted",
                          "median",
                          "median.weighted",
                          "kde.max",
                          "sd.abs",
                          "sd.rel",
                          "se.abs",
                          "se.rel",
                          "q25",
                          "q75",
                          "skewness",
                          "kurtosis",
                          "sd.abs.weighted",
                          "sd.rel.weighted",
                          "se.abs.weighted",
                          "se.rel.weighted")

  for(i in 1:length(data)) {
    statistics <- calc_Statistics(data[[i]])
    statistics.2 <- calc_Statistics(data[[i]][,3:4])
    
    De.stats[i,1] <- statistics$weighted$n
    De.stats[i,2] <- statistics.2$unweighted$mean
    De.stats[i,3] <- statistics.2$weighted$mean
    De.stats[i,4] <- statistics.2$unweighted$median
    De.stats[i,7] <- statistics$unweighted$sd.abs
    De.stats[i,8] <- statistics$unweighted$sd.rel
    De.stats[i,9] <- statistics$unweighted$se.abs
    De.stats[i,10] <- statistics$weighted$se.rel
    De.stats[i,11] <- quantile(data[[i]][,1], 0.25)
    De.stats[i,12] <- quantile(data[[i]][,1], 0.75)
    De.stats[i,13] <- statistics$unweighted$skewness
    De.stats[i,14] <- statistics$unweighted$kurtosis
    De.stats[i,15] <- statistics$weighted$sd.abs
    De.stats[i,16] <- statistics$weighted$sd.rel
    De.stats[i,17] <- statistics$weighted$se.abs
    De.stats[i,18] <- statistics$weighted$se.rel
    
    ## account for log.z-option
    if(log.z == TRUE) {
      De.stats[i,2:4] <- exp(De.stats[i,2:4])
    }

    ##kdemax - here a little doubled as it appears below again
    De.density <-density(x = data[[i]][,1],
                         kernel = "gaussian",
                         bw = bw,
                         from = limits.z[1],
                         to = limits.z[2])

    De.stats[i,6] <- De.density$x[which.max(De.density$y)]
  }

  label.text = list(NA)

  if(summary.pos[1] != "sub") {
    n.rows <- length(summary)

    for(i in 1:length(data)) {
      stops <- paste(rep("\n", (i - 1) * n.rows), collapse = "")

      summary.text <- character(0)

      for(j in 1:length(summary)) {
        summary.text <- c(summary.text,
                          paste(
                            "",
                            ifelse("n" %in% summary[j] == TRUE,
                                   paste("n = ",
                                         De.stats[i,1],
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("mean" %in% summary[j] == TRUE,
                                   paste("mean = ",
                                         round(De.stats[i,2], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("mean.weighted" %in% summary[j] == TRUE,
                                   paste("weighted mean = ",
                                         round(De.stats[i,3], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("median" %in% summary[j] == TRUE,
                                   paste("median = ",
                                         round(De.stats[i,4], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("kdemax" %in% summary[j] == TRUE,
                                   paste("kdemax = ",
                                         round(De.stats[i,6], 2),
                                         " \n ",
                                         sep = ""),
                                   ""),
                            ifelse("sdabs" %in% summary[j] == TRUE,
                                   paste("sd = ",
                                         round(De.stats[i,7], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("sdrel" %in% summary[j] == TRUE,
                                   paste("rel. sd = ",
                                         round(De.stats[i,8], 2), " %",
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("seabs" %in% summary[j] == TRUE,
                                   paste("se = ",
                                         round(De.stats[i,9], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("serel" %in% summary[j] == TRUE,
                                   paste("rel. se = ",
                                         round(De.stats[i,10], 2), " %",
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("skewness" %in% summary[j] == TRUE,
                                   paste("skewness = ",
                                         round(De.stats[i,13], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("kurtosis" %in% summary[j] == TRUE,
                                   paste("kurtosis = ",
                                         round(De.stats[i,14], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("in.2s" %in% summary[j] == TRUE,
                                   paste("in 2 sigma = ",
                                         round(sum(data[[i]][,7] > -2 &
                                                     data[[i]][,7] < 2) /
                                                 nrow(data[[i]]) * 100 , 1),
                                         " %",
                                         sep = ""),
                                   ""),
                            ifelse("sdabs.weighted" %in% summary[j] == TRUE,
                                   paste("abs. weighted sd = ",
                                         round(De.stats[i,15], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("sdrel.weighted" %in% summary[j] == TRUE,
                                   paste("rel. weighted sd = ",
                                         round(De.stats[i,16], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("seabs.weighted" %in% summary[j] == TRUE,
                                   paste("abs. weighted se = ",
                                         round(De.stats[i,17], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("serel.weighted" %in% summary[j] == TRUE,
                                   paste("rel. weighted se = ",
                                         round(De.stats[i,18], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            sep = ""))
      }

      summary.text <- paste(summary.text, collapse = "")

      label.text[[length(label.text) + 1]] <- paste(stops,
                                                    summary.text,
                                                    stops,
                                                    sep = "")
    }
  } else {
    for(i in 1:length(data)) {

      summary.text <- character(0)

      for(j in 1:length(summary)) {
        summary.text <- c(summary.text,
                          ifelse("n" %in% summary[j] == TRUE,
                                 paste("n = ",
                                       De.stats[i,1],
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("mean" %in% summary[j] == TRUE,
                                 paste("mean = ",
                                       round(De.stats[i,2], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("mean.weighted" %in% summary[j] == TRUE,
                                 paste("weighted mean = ",
                                       round(De.stats[i,3], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("median" %in% summary[j] == TRUE,
                                 paste("median = ",
                                       round(De.stats[i,4], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("kdemax" %in% summary[j] == TRUE,
                                 paste("kdemax = ",
                                       round(De.stats[i,6], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("sdrel" %in% summary[j] == TRUE,
                                 paste("rel. sd = ",
                                       round(De.stats[i,8], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("sdabs" %in% summary[j] == TRUE,
                                 paste("abs. sd = ",
                                       round(De.stats[i,7], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("serel" %in% summary[j] == TRUE,
                                 paste("rel. se = ",
                                       round(De.stats[i,10], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("seabs" %in% summary[j] == TRUE,
                                 paste("abs. se = ",
                                       round(De.stats[i,9], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("skewness" %in% summary[j] == TRUE,
                                 paste("skewness = ",
                                       round(De.stats[i,13], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("kurtosis" %in% summary[j] == TRUE,
                                 paste("kurtosis = ",
                                       round(De.stats[i,14], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("in.2s" %in% summary[j] == TRUE,
                                 paste("in 2 sigma = ",
                                       round(sum(data[[i]][,7] > -2 &
                                                   data[[i]][,7] < 2) /
                                               nrow(data[[i]]) * 100 , 1),
                                       " %   ",
                                       sep = ""),
                                 ""),
                          ifelse("sdabs.weighted" %in% summary[j] == TRUE,
                                 paste("abs. weighted sd = ",
                                       round(De.stats[i,15], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("sdrel.weighted" %in% summary[j] == TRUE,
                                 paste("rel. weighted sd = ",
                                       round(De.stats[i,16], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("seabs.weighted" %in% summary[j] == TRUE,
                                 paste("abs. weighted se = ",
                                       round(De.stats[i,17], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("serel.weighted" %in% summary[j] == TRUE,
                                 paste("rel. weighted se = ",
                                       round(De.stats[i,18], 2), " %",
                                       " | ",
                                       sep = ""),
                                 "")
                          )
      }

      summary.text <- paste(summary.text, collapse = "")

      label.text[[length(label.text) + 1]]  <- paste(
        "  ",
        summary.text,
        sep = "")
    }

    ## remove outer vertical lines from string
    for(i in 2:length(label.text)) {
      label.text[[i]] <- substr(x = label.text[[i]],
                                start = 3,
                                stop = nchar(label.text[[i]]) - 3)
    }
  }

  ## remove dummy list element
  label.text[[1]] <- NULL

  if(rotate == FALSE) {
    ## convert keywords into summary placement coordinates
    if(missing(summary.pos) == TRUE) {
      summary.pos <- c(limits.x[1], limits.y[2])
      summary.adj <- c(0, 1)
    } else if(length(summary.pos) == 2) {
      summary.pos <- summary.pos
      summary.adj <- c(0, 1)
    } else if(summary.pos[1] == "topleft") {
      summary.pos <- c(limits.x[1], limits.y[2] - par()$cxy[2] * 1)
      summary.adj <- c(0, 1)
    } else if(summary.pos[1] == "top") {
      summary.pos <- c(mean(limits.x), limits.y[2] - par()$cxy[2] * 1)
      summary.adj <- c(0.5, 1)
    } else if(summary.pos[1] == "topright") {
      summary.pos <- c(limits.x[2], limits.y[2] - par()$cxy[2] * 1)
      summary.adj <- c(1, 1)
    }  else if(summary.pos[1] == "left") {
      summary.pos <- c(limits.x[1], mean(limits.y))
      summary.adj <- c(0, 0.5)
    } else if(summary.pos[1] == "center") {
      summary.pos <- c(mean(limits.x), mean(limits.y))
      summary.adj <- c(0.5, 0.5)
    } else if(summary.pos[1] == "right") {
      summary.pos <- c(limits.x[2], mean(limits.y))
      summary.adj <- c(1, 0.5)
    }else if(summary.pos[1] == "bottomleft") {
      summary.pos <- c(limits.x[1], limits.y[1] + par()$cxy[2] * 3.5)
      summary.adj <- c(0, 0)
    } else if(summary.pos[1] == "bottom") {
      summary.pos <- c(mean(limits.x), limits.y[1] + par()$cxy[2] * 3.5)
      summary.adj <- c(0.5, 0)
    } else if(summary.pos[1] == "bottomright") {
      summary.pos <- c(limits.x[2], limits.y[1] + par()$cxy[2] * 3.5)
      summary.adj <- c(1, 0)
    }

    ## convert keywords into legend placement coordinates
    if(missing(legend.pos) == TRUE) {
      legend.pos <- c(limits.x[1], limits.y[2])
      legend.adj <- c(0, 1)
    } else if(length(legend.pos) == 2) {
      legend.pos <- legend.pos
      legend.adj <- c(0, 1)
    } else if(legend.pos[1] == "topleft") {
      legend.pos <- c(limits.x[1], limits.y[2])
      legend.adj <- c(0, 1)
    } else if(legend.pos[1] == "top") {
      legend.pos <- c(mean(limits.x), limits.y[2])
      legend.adj <- c(0.5, 1)
    } else if(legend.pos[1] == "topright") {
      legend.pos <- c(limits.x[2], limits.y[2])
      legend.adj <- c(1, 1)
    } else if(legend.pos[1] == "left") {
      legend.pos <- c(limits.x[1], mean(limits.y))
      legend.adj <- c(0, 0.5)
    } else if(legend.pos[1] == "center") {
      legend.pos <- c(mean(limits.x), mean(limits.y))
      legend.adj <- c(0.5, 0.5)
    } else if(legend.pos[1] == "right") {
      legend.pos <- c(limits.x[2], mean(limits.y))
      legend.adj <- c(1, 0.5)
    } else if(legend.pos[1] == "bottomleft") {
      legend.pos <- c(limits.x[1], limits.y[1])
      legend.adj <- c(0, 0)
    } else if(legend.pos[1] == "bottom") {
      legend.pos <- c(mean(limits.x), limits.y[1])
      legend.adj <- c(0.5, 0)
    } else if(legend.pos[1] == "bottomright") {
      legend.pos <- c(limits.x[2], limits.y[1])
      legend.adj <- c(1, 0)
    }
  } else {
    ## convert keywords into summary placement coordinates
    if(missing(summary.pos) == TRUE) {
      summary.pos <- c(limits.y[1] + par()$cxy[1] * 7.5, limits.x[1])
      summary.adj <- c(0, 0)
    } else if(length(summary.pos) == 2) {
      summary.pos <- summary.pos
      summary.adj <- c(0, 1)
    } else if(summary.pos[1] == "topleft") {
      summary.pos <- c(limits.y[1] + par()$cxy[1] * 7.5, limits.x[2])
      summary.adj <- c(0, 1)
    } else if(summary.pos[1] == "top") {
      summary.pos <- c(mean(limits.y), limits.x[2])
      summary.adj <- c(0.5, 1)
    } else if(summary.pos[1] == "topright") {
      summary.pos <- c(limits.y[2], limits.x[2])
      summary.adj <- c(1, 1)
    }  else if(summary.pos[1] == "left") {
      summary.pos <- c(limits.y[1] + par()$cxy[1] * 7.5, mean(limits.x))
      summary.adj <- c(0, 0.5)
    } else if(summary.pos[1] == "center") {
      summary.pos <- c(mean(limits.y), mean(limits.x))
      summary.adj <- c(0.5, 0.5)
    } else if(summary.pos[1] == "right") {
      summary.pos <- c(limits.y[2], mean(limits.x))
      summary.adj <- c(1, 0.5)
    }else if(summary.pos[1] == "bottomleft") {
      summary.pos <- c(limits.y[1] + par()$cxy[1] * 7.5, limits.x[1])
      summary.adj <- c(0, 0)
    } else if(summary.pos[1] == "bottom") {
      summary.pos <- c(mean(limits.y), limits.x[1])
      summary.adj <- c(0.5, 0)
    } else if(summary.pos[1] == "bottomright") {
      summary.pos <- c(limits.y[2], limits.x[1])
      summary.adj <- c(1, 0)
    }

    ## convert keywords into legend placement coordinates
    if(missing(legend.pos) == TRUE) {
      legend.pos <- c(limits.y[1] + par()$cxy[1] * 7.5, limits.x[1])
      legend.adj <- c(0, 0)
    } else if(length(legend.pos) == 2) {
      legend.pos <- legend.pos
      legend.adj <- c(1, 0)
    } else if(legend.pos[1] == "topleft") {
      legend.pos <- c(limits.y[1] + par()$cxy[1] * 11, limits.x[2])
      legend.adj <- c(1, 0)
    } else if(legend.pos[1] == "top") {
      legend.pos <- c(mean(limits.y), limits.x[2])
      legend.adj <- c(1, 0.5)
    } else if(legend.pos[1] == "topright") {
      legend.pos <- c(limits.y[2], limits.x[2])
      legend.adj <- c(1, 1)
    } else if(legend.pos[1] == "left") {
      legend.pos <- c(limits.y[1] + par()$cxy[1] * 7.5, mean(limits.x))
      legend.adj <- c(0.5, 0)
    } else if(legend.pos[1] == "center") {
      legend.pos <- c(mean(limits.y), mean(limits.x))
      legend.adj <- c(0.5, 0.5)
    } else if(legend.pos[1] == "right") {
      legend.pos <- c(limits.y[2], mean(limits.x))
      legend.adj <- c(0.5, 1)
    } else if(legend.pos[1] == "bottomleft") {
      legend.pos <- c(limits.y[1] + par()$cxy[1] * 7.5, limits.x[1])
      legend.adj <- c(0, 0)
    } else if(legend.pos[1] == "bottom") {
      legend.pos <- c(mean(limits.y), limits.x[1])
      legend.adj <- c(0, 0.5)
    } else if(legend.pos[1] == "bottomright") {
      legend.pos <- c(limits.y[2], limits.x[1])
      legend.adj <- c(0, 1)
    }
  }

  ## define cartesian plot origins
  if(rotate == FALSE) {
    xy.0 <- c(min(ellipse[,1]) * 1.03, min(ellipse[,2]))
  } else {
    xy.0 <- c(min(ellipse[,1]), min(ellipse[,2]) * 1.03)
  }

  ## calculate coordinates for dispersion polygon overlay
  y.max.x <- 2 * limits.x[2] / max(data.global[6])

  polygons <- matrix(nrow = length(data), ncol = 14)
  for(i in 1:length(data)) {

    if(dispersion == "qr") {
      ci.lower <- quantile(data[[i]][,1], 0.25)
      ci.upper <- quantile(data[[i]][,1], 0.75)
    } else if(grepl(x = dispersion, pattern = "p") == TRUE) {
      ci.plot <- as.numeric(strsplit(x = dispersion,
                                     split = "p")[[1]][2])
      ci.plot <- (100 - ci.plot) / 100
      ci.lower <- quantile(data[[i]][,1], ci.plot)
      ci.upper <- quantile(data[[i]][,1], 1 - ci.plot)
    } else if(dispersion == "sd") {
      if(log.z == TRUE) {
        ci.lower <- exp(mean(log(data[[i]][,1])) - sd(log(data[[i]][,1])))
        ci.upper <- exp(mean(log(data[[i]][,1])) + sd(log(data[[i]][,1])))
      } else {
        ci.lower <- mean(data[[i]][,1]) - sd(data[[i]][,1])
        ci.upper <- mean(data[[i]][,1]) + sd(data[[i]][,1])
      }
    } else if(dispersion == "2sd") {
      if(log.z == TRUE) {
        ci.lower <- exp(mean(log(data[[i]][,1])) - 2 * sd(log(data[[i]][,1])))
        ci.upper <- exp(mean(log(data[[i]][,1])) + 2 * sd(log(data[[i]][,1])))
      } else {
        ci.lower <- mean(data[[i]][,1]) - 2 * sd(data[[i]][,1])
        ci.upper <- mean(data[[i]][,1]) + 2 * sd(data[[i]][,1])
      }
    } else {
      stop("Measure of dispersion not supported.")
    }

    if(log.z == TRUE) {
      ci.lower[which(ci.lower < 0)] <- 1
      y.lower <- log(ci.lower)
      y.upper <- log(ci.upper)
    } else {
      y.lower <- ci.lower
      y.upper <- ci.upper
    }

    if(rotate == FALSE) {
      polygons[i,1:7] <- c(limits.x[1],
                           limits.x[2],
                           xy.0[1],
                           par()$usr[2],
                           par()$usr[2],
                           xy.0[1],
                           limits.x[2])
      polygons[i,8:14] <- c(0,
                            (y.upper - z.central.global) *
                              limits.x[2],
                            (y.upper - z.central.global) *
                              xy.0[1],
                            (y.upper - z.central.global) *
                              xy.0[1],
                            (y.lower - z.central.global) *
                              xy.0[1],
                            (y.lower - z.central.global) *
                              xy.0[1],
                            (y.lower - z.central.global) *
                              limits.x[2]
      )
    } else {
      y.max <- par()$usr[4]
      polygons[i,1:7] <- c(limits.x[1],
                           limits.x[2],
                           xy.0[2],
                           y.max,
                           y.max,
                           xy.0[2],
                           limits.x[2])
      polygons[i,8:14] <- c(0,
                            (y.upper - z.central.global) *
                              limits.x[2],
                            (y.upper - z.central.global) *
                              xy.0[2],
                            (y.upper - z.central.global) *
                              xy.0[2],
                            (y.lower - z.central.global) *
                              xy.0[2],
                            (y.lower - z.central.global) *
                              xy.0[2],
                            (y.lower - z.central.global) *
                              limits.x[2]
      )
    }
  }

  ## append information about data in confidence interval
  for(i in 1:length(data)) {
    data.in.2s <- rep(x = FALSE, times = nrow(data[[i]]))
    data.in.2s[data[[i]][,8] > -2 & data[[i]][,8] < 2] <- TRUE
    data[[i]] <- cbind(data[[i]], data.in.2s)
  }

  ## calculate coordinates for 2-sigma bar overlay
  if(bar[1] == TRUE) {
    bars <- matrix(nrow = length(data), ncol = 8)

    for(i in 1:length(data)) {
      bars[i,1:4] <- c(limits.x[1],
                       limits.x[1],
                       ifelse("xlim" %in% names(extraArgs),
                              extraArgs$xlim[2] * 0.95,
                              max(data.global$precision)),
                       ifelse("xlim" %in% names(extraArgs),
                              extraArgs$xlim[2] * 0.95,
                              max(data.global$precision)))

      bars[i,5:8] <- c(-2,
                       2,
                       (data[[i]][1,5] - z.central.global) *
                         bars[i,3] + 2,
                       (data[[i]][1,5] - z.central.global) *
                         bars[i,3] - 2)

    }
  } else {
    bars <- matrix(nrow = length(bar), ncol = 8)

    if(is.numeric(bar) == TRUE & log.z == TRUE) {
      bar <- log(bar)
    }

    for(i in 1:length(bar)) {
      bars[i,1:4] <- c(limits.x[1],
                       limits.x[1],
                       ifelse("xlim" %in% names(extraArgs),
                              extraArgs$xlim[2] * 0.95,
                              max(data.global$precision)),
                       ifelse("xlim" %in% names(extraArgs),
                              extraArgs$xlim[2] * 0.95,
                              max(data.global$precision)))

      bars[i,5:8] <- c(-2,
                       2,
                       (bar[i] - z.central.global) *
                         bars[i,3] + 2,
                       (bar[i] - z.central.global) *
                         bars[i,3] - 2)
    }
  }

  ## calculate error bar coordinates
  if(error.bars == TRUE) {
    arrow.coords <- list(NA)
    for(i in 1:length(data)) {
      arrow.x1 <- data[[i]][,6]
      arrow.x2 <- data[[i]][,6]
      arrow.y1 <- data[[i]][,1] - data[[i]][,2]
      arrow.y2 <- data[[i]][,1] + data[[i]][,2]

      if(log.z == TRUE) {
        arrow.y1 <- log(arrow.y1)
        arrow.y2 <- log(arrow.y2)
      }

      arrow.coords[[length(arrow.coords) + 1]] <- cbind(
        arrow.x1,
        arrow.x2,
        (arrow.y1 - z.central.global) * arrow.x1,
        (arrow.y2 - z.central.global) * arrow.x1)
    }
    arrow.coords[[1]] <- NULL
  }

  ## calculate KDE
  KDE <- list(NA)
  KDE.ext <- 0
  KDE.bw <- numeric(0)

  for(i in 1:length(data)) {
    KDE.i <- density(x = data[[i]][,3],
                     kernel = "gaussian",
                     bw = bw,
                     from = ellipse.values[1],
                     to = ellipse.values[2],
                     weights = data[[i]]$weights)
    KDE.xy <- cbind(KDE.i$x, KDE.i$y)
    KDE.bw <- c(KDE.bw, KDE.i$bw)
    KDE.ext <- ifelse(max(KDE.xy[,2]) < KDE.ext, KDE.ext, max(KDE.xy[,2]))
    KDE.xy <- rbind(c(min(KDE.xy[,1]), 0), KDE.xy, c(max(KDE.xy[,1]), 0))
    KDE[[length(KDE) + 1]] <- cbind(KDE.xy[,1], KDE.xy[,2])
  }
  KDE[1] <- NULL

  ## calculate mean KDE bandwidth
  KDE.bw <- mean(KDE.bw, na.rm = TRUE)

  ## calculate max KDE value for labelling
  KDE.max.plot <- numeric(length(data))

  for(i in 1:length(data)) {
    KDE.plot <- density(x = data[[i]][,1],
                        kernel = "gaussian",
                        bw = bw,
                        from = limits.z[1],
                        to = limits.z[2])
    KDE.max.plot[i] <- max(KDE.plot$y)
  }

  KDE.max.plot <- max(KDE.max.plot, na.rm = TRUE)

  ## calculate histogram data without plotting

  ## create dummy list
  hist.data <- list(NA)

  for(i in 1:length(data)) {
    hist.i <- hist(x = data[[i]][,3],
                   plot = FALSE,
                   breaks = breaks)
    hist.data[[length(hist.data) + 1]] <- hist.i
  }

  ## remove dummy list object
  hist.data[[1]] <- NULL

  ## calculate maximum histogram bar height for normalisation
  hist.max.plot <- numeric(length(data))
  for(i in 1:length(data)) {
    hist.max.plot <- ifelse(max(hist.data[[i]]$counts, na.rm = TRUE) >
                              hist.max.plot, max(hist.data[[i]]$counts,
                                                 na.rm = TRUE), hist.max.plot)
  }
  hist.max.plot <- max(hist.max.plot, na.rm = TRUE)

  ## normalise histogram bar height to KDE dimensions
  for(i in 1:length(data)) {
    hist.data[[i]]$density <- hist.data[[i]]$counts / hist.max.plot *
      KDE.max.plot
  }

  ## calculate line coordinates and further parameters
  if(missing(line) == FALSE) {

    ## check if line parameters are R.Lum-objects
    for(i in 1:length(line)) {
      if(is.list(line) == TRUE) {
        if(is(line[[i]], "RLum.Results")) {
          line[[i]] <- as.numeric(get_RLum(object = line[[i]],
                                                   data.object = "summary")$de)
        }
      } else if(is(object = line, class2 = "RLum.Results")) {
        line <- as.numeric(get_RLum(object = line,
                                            data.object = "summary")$de)
      }
    }

    ## convert list to vector
    if(is.list(line) == TRUE) {
      line <- unlist(line)
    }

    if(log.z == TRUE) {
      line <- log(line)
    }

    line.coords <- list(NA)

    if(rotate == FALSE) {
      for(i in 1:length(line)) {
        line.x <- c(limits.x[1], min(ellipse[,1]), par()$usr[2])
        line.y <- c(0,
                    (line[i] - z.central.global) * min(ellipse[,1]),
                    (line[i] - z.central.global) * min(ellipse[,1]))
        line.coords[[length(line.coords) + 1]] <- rbind(line.x, line.y)
      }
    } else {
      for(i in 1:length(line)) {
        line.x <- c(limits.x[1], min(ellipse[,2]),y.max)
        line.y <- c(0,
                    (line[i] - z.central.global) * min(ellipse[,2]),
                    (line[i] - z.central.global) * min(ellipse[,2]))
        line.coords[[length(line.coords) + 1]] <- rbind(line.x, line.y)
      }
    }

    line.coords[1] <- NULL

    if(missing(line.col) == TRUE) {
      line.col <- seq(from = 1, to = length(line.coords))
    }

    if(missing(line.label) == TRUE) {
      line.label <- rep("", length(line.coords))
    }
  }

  ## calculate rug coordinates
  if(missing(rug) == FALSE) {
    if(log.z == TRUE) {
      rug.values <- log(De.global)
    } else {
      rug.values <- De.global
    }

    rug.coords <- list(NA)

    if(rotate == FALSE) {
      for(i in 1:length(rug.values)) {
        rug.x <- c(xy.0[1] * (1 - 0.013 * (layout$abanico$dimension$rugl / 100)),
                   xy.0[1])
        rug.y <- c((rug.values[i] - z.central.global) * min(ellipse[,1]),
                   (rug.values[i] - z.central.global) * min(ellipse[,1]))
        rug.coords[[length(rug.coords) + 1]] <- rbind(rug.x, rug.y)
      }
    } else {
      for(i in 1:length(rug.values)) {
        rug.x <- c(xy.0[2] * (1 - 0.013 * (layout$abanico$dimension$rugl / 100)),
                   xy.0[2])
        rug.y <- c((rug.values[i] - z.central.global) * min(ellipse[,2]),
                   (rug.values[i] - z.central.global) * min(ellipse[,2]))
        rug.coords[[length(rug.coords) + 1]] <- rbind(rug.x, rug.y)
      }
    }

    rug.coords[1] <- NULL
  }

  ## Generate plot ------------------------------------------------------------

  ## determine number of subheader lines to shift the plot
  if(length(summary) > 0 & summary.pos[1] == "sub") {
    shift.lines <- (length(data) + 1) * layout$abanico$dimension$summary.line/100
  } else {shift.lines <- 1}

  ## extract original plot parameters
  par(bg = layout$abanico$colour$background)
  bg.original <- par()$bg

  if(rotate == FALSE) {
    ## setup plot area
    par(mar = c(4.5, 4.5, shift.lines + 1.5, 7),
        xpd = TRUE,
        cex = cex)

    if(layout$abanico$dimension$figure.width != "auto" |
         layout$abanico$dimension$figure.height != "auto") {
      par(mai = layout$abanico$dimension$margin / 25.4,
          pin = c(layout$abanico$dimension$figure.width / 25.4 -
                    layout$abanico$dimension$margin[2] / 25.4 -
                    layout$abanico$dimension$margin[4] / 25.4,
                  layout$abanico$dimension$figure.height / 25.4 -
                    layout$abanico$dimension$margin[1] / 25.4 -
                    layout$abanico$dimension$margin[3]/25.4))
    }

    ## create empty plot
    par(new = TRUE)
    plot(NA,
         xlim = c(limits.x[1], limits.x[2] * (1 / plot.ratio)),
         ylim = limits.y,
         main = "",
         sub = sub,
         xlab = "",
         ylab = "",
         xaxs = "i",
         yaxs = "i",
         frame.plot = FALSE,
         axes = FALSE)

    ## add y-axis label
    mtext(text = ylab,
          at = mean(x = c(min(ellipse[,2]),
                          max(ellipse[,2])),
                    na.rm = TRUE),
          #        at = 0, ## BUG FROM VERSION 0.4.0, maybe removed in future
          adj = 0.5,
          side = 2,
          line = 3 * layout$abanico$dimension$ylab.line / 100,
          col = layout$abanico$colour$ylab,
          family = layout$abanico$font.type$ylab,
          font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                         layout$abanico$font.deco$ylab],
          cex = cex * layout$abanico$font.size$ylab/12)

    ## calculate upper x-axis label values
    label.x.upper <- if(log.z == TRUE) {
      as.character(round(1/axTicks(side = 1)[-1] * 100, 1))
    } else {
      as.character(round(1/axTicks(side = 1)[-1], 1))
    }

    # optionally, plot 2-sigma-bar
    if(bar[1] != FALSE) {
      for(i in 1:length(bar)) {
        polygon(x = bars[i,1:4],
                y = bars[i,5:8],
                col = bar.fill[i],
                border = bar.line[i])
      }
    }

    ## remove unwanted parts
    polygon(x = c(par()$usr[2],
                  par()$usr[2],
                  par()$usr[2] * 2,
                  par()$usr[2] * 2),
            y = c(min(ellipse[,2]) * 2,
                  max(ellipse[,2]) * 2,
                  max(ellipse[,2]) * 2,
                  min(ellipse[,2]) * 2),
            col = bg.original,
            lty = 0)

    ## optionally, plot dispersion polygon
    if(polygon.fill[1] != "none") {
      for(i in 1:length(data)) {
        polygon(x = polygons[i,1:7],
                y = polygons[i,8:14],
                col = polygon.fill[i],
                border = polygon.line[i])
      }
    }

    ## optionally, add minor grid lines
    if(grid.minor != "none") {

      for(i in 1:length(tick.values.minor)) {
        lines(x = c(limits.x[1], min(ellipse[,1])),
              y = c(0, (tick.values.minor[i] - z.central.global) *
                      min(ellipse[,1])),
              col = grid.minor,
              lwd = 1)
      }

      for(i in 1:length(tick.values.minor)) {
        lines(x = c(xy.0[1], par()$usr[2]),
              y = c((tick.values.minor[i] - z.central.global) *
                      min(ellipse[,1]),
                    (tick.values.minor[i] - z.central.global) *
                      min(ellipse[,1])),
              col = grid.minor,
              lwd = 1)
      }
    }

    ## optionally, add major grid lines
    if(grid.major != "none") {
      for(i in 1:length(tick.values.major)) {
        lines(x = c(limits.x[1], min(ellipse[,1])),
              y = c(0, (tick.values.major[i] - z.central.global) *
                      min(ellipse[,1])),
              col = grid.major,
              lwd = 1)
      }
      for(i in 1:length(tick.values.major)) {
        lines(x = c(xy.0[1], par()$usr[2]),
              y = c((tick.values.major[i] - z.central.global) *
                      min(ellipse[,1]),
                    (tick.values.major[i] - z.central.global) *
                      min(ellipse[,1])),
              col = grid.major,
              lwd = 1)
      }
    }

    ## optionally, plot lines for each bar
    if(lwd[1] > 0 & lty[1] > 0 & bar[1] != FALSE & length(data) == 1) {
      if(bar[1] == TRUE & length(bar) == 1) {
        bar[1] <- z.central.global
      }
      for(i in 1:length(bar)) {
        x2 <- r / sqrt(1 + f^2 * (
          bar[i] - z.central.global)^2)
        y2 <- (bar[i] - z.central.global) * x2
        lines(x = c(limits.x[1], x2, xy.0[1], par()$usr[2]),
              y = c(0, y2, y2, y2),
              lty = lty[i],
              lwd = lwd[i],
              col = centrality.col[i])
      }
    } else if(lwd[1] > 0 & lty[1] > 0 & bar[1] != FALSE) {
      for(i in 1:length(data)) {

        z.line <- ifelse(test = is.numeric(bar[i]) == TRUE,
                         yes = bar[i],
                         no = data[[i]][1,5])

        x2 <- r / sqrt(1 + f^2 * (
          z.line - z.central.global)^2)
        y2 <- (z.line - z.central.global) * x2
        lines(x = c(limits.x[1], x2, xy.0[1], par()$usr[2]),
              y = c(0, y2, y2, y2),
              lty = lty[i],
              lwd = lwd[i],
              col = centrality.col[i])
      }
    }

    ## optionally add further lines
    if(missing(line) == FALSE) {
      for(i in 1:length(line)) {
        lines(x = line.coords[[i]][1,1:3],
              y = line.coords[[i]][2,1:3],
              col = line.col[i])
        text(x = line.coords[[i]][1,3],
             y = line.coords[[i]][2,3] + par()$cxy[2] * 0.3,
             labels = line.label[i],
             pos = 2,
             col = line.col[i],
             cex = cex * 0.9)
      }
    }

    ## add plot title
    cex.old <- par()$cex
    par(cex = layout$abanico$font.size$main / 12)
    title(main = main,
          family = layout$abanico$font.type$main,
          font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                         layout$abanico$font.deco$main],
          col.main = layout$abanico$colour$main,
          line = shift.lines * layout$abanico$dimension$main / 100)
    par(cex = cex.old)

    ## calculate lower x-axis (precision)
    x.axis.ticks <- axTicks(side = 1)
    x.axis.ticks <- x.axis.ticks[c(TRUE, x.axis.ticks <= limits.x[2])]
    x.axis.ticks <- x.axis.ticks[x.axis.ticks <= max(ellipse[,1])]

    ## x-axis with lables and ticks
    axis(side = 1,
         at = x.axis.ticks,
         col = layout$abanico$colour$xtck1,
         col.axis = layout$abanico$colour$xtck1,
         labels = NA,
         tcl = -layout$abanico$dimension$xtcl1 / 200,
         cex = cex)
    axis(side = 1,
         at = x.axis.ticks,
         line = 2 * layout$abanico$dimension$xtck1.line / 100 - 2,
         lwd = 0,
         col = layout$abanico$colour$xtck1,
         family = layout$abanico$font.type$xtck1,
         font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                        layout$abanico$font.deco$xtck1],
         col.axis = layout$abanico$colour$xtck1,
         cex.axis = layout$abanico$font.size$xlab1/12)

    ## extend axis line to right side of the plot
    lines(x = c(max(x.axis.ticks), max(ellipse[,1])),
          y = c(limits.y[1], limits.y[1]),
          col = layout$abanico$colour$xtck1)

    ## draw closing tick on right hand side
    axis(side = 1,
         tcl = -layout$abanico$dimension$xtcl1 / 200,
         lwd = 0,
         lwd.ticks = 1,
         at = limits.x[2],
         labels = FALSE,
         col = layout$abanico$colour$xtck1)

    axis(side = 1,
         tcl = layout$abanico$dimension$xtcl2 / 200,
         lwd = 0,
         lwd.ticks = 1,
         at = limits.x[2],
         labels = FALSE,
         col = layout$abanico$colour$xtck2)

    ## add lower axis label
    mtext(xlab[2],
          at = (limits.x[1] + max(ellipse[,1])) / 2,
          side = 1,
          line = 2.5 * layout$abanico$dimension$xlab1.line / 100,
          col = layout$abanico$colour$xlab1,
          family = layout$abanico$font.type$xlab1,
          font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                         layout$abanico$font.deco$xlab1],
          cex = cex * layout$abanico$font.size$xlab1/12)

    ## add upper axis label
    mtext(xlab[1],
          at = (limits.x[1] + max(ellipse[,1])) / 2,
          side = 1,
          line = -3.5 * layout$abanico$dimension$xlab2.line / 100,
          col = layout$abanico$colour$xlab2,
          family = layout$abanico$font.type$xlab2,
          font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                         layout$abanico$font.deco$xlab2],
          cex = cex * layout$abanico$font.size$xlab2/12)

    ## plot upper x-axis
    axis(side = 1,
         at = x.axis.ticks[-1],
         col = layout$abanico$colour$xtck2,
         col.axis = layout$abanico$colour$xtck2,
         labels = NA,
         tcl = layout$abanico$dimension$xtcl2 / 200,
         cex = cex)

    ## remove first tick label (infinity)
    label.x.upper <- label.x.upper[1:(length(x.axis.ticks) - 1)]

    axis(side = 1,
         at = x.axis.ticks[-1],
         labels = label.x.upper,
         line = -1 * layout$abanico$dimension$xtck2.line / 100 - 2,
         lwd = 0,
         col = layout$abanico$colour$xtck2,
         family = layout$abanico$font.type$xtck2,
         font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                        layout$abanico$font.deco$xtck2],
         col.axis = layout$abanico$colour$xtck2,
         cex.axis = layout$abanico$font.size$xlab2/12)

    ## plot y-axis
    if(y.axis == TRUE) {
      char.height <- par()$cxy[2]
      tick.space <- axisTicks(usr = limits.y, log = FALSE)
      tick.space <- (max(tick.space) - min(tick.space)) / length(tick.space)
      if(tick.space < char.height * 1.7) {
        axis(side = 2,
             tcl = -layout$abanico$dimension$ytcl / 200,
             lwd = 1,
             lwd.ticks = 1,
             at = c(-2, 2),
             labels = c("", ""),
             las = 1,
             col = layout$abanico$colour$ytck)

        axis(side = 2,
             at = 0,
             tcl = 0,
             line = 2 * layout$abanico$dimension$ytck.line / 100 - 2,
             labels = paste("\u00B1", "2"),
             las = 1,
             family = layout$abanico$font.type$ytck,
             font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                            layout$abanico$font.deco$ytck],
             col.axis = layout$abanico$colour$ytck,
             cex.axis = layout$abanico$font.size$ylab/12)
      } else {
        axis(side = 2,
             at = seq(-2, 2, by = 2),
             col = layout$abanico$colour$ytck,
             col.axis = layout$abanico$colour$ytck,
             labels = NA,
             las = 1,
             tcl = -layout$abanico$dimension$ytcl / 200,
             cex = cex)
        axis(side = 2,
             at = seq(-2, 2, by = 2),
             line = 2 * layout$abanico$dimension$ytck.line / 100 - 2,
             lwd = 0,
             las = 1,
             col = layout$abanico$colour$ytck,
             family = layout$abanico$font.type$ytck,
             font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                            layout$abanico$font.deco$ytck],
             col.axis = layout$abanico$colour$ytck,
             cex.axis = layout$abanico$font.size$ylab/12)
      }
    } else {
      axis(side = 2,
           at = 0,
           col = layout$abanico$colour$ytck,
           col.axis = layout$abanico$colour$ytck,
           labels = NA,
           las = 1,
           tcl = -layout$abanico$dimension$ytcl / 200,
           cex = cex)
      axis(side = 2,
           at = 0,
           line = 2 * layout$abanico$dimension$ytck.line / 100 - 2,
           lwd = 0,
           las = 1,
           col = layout$abanico$colour$ytck,
           family = layout$abanico$font.type$ytck,
           font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                          layout$abanico$font.deco$ytck],
           col.axis = layout$abanico$colour$ytck,
           cex.axis = layout$abanico$font.size$ylab/12)
    }

    ## plot minor z-ticks
    for(i in 1:length(tick.values.minor)) {
      lines(x = c(par()$usr[2],
                  (1 + 0.007 * cex * layout$abanico$dimension$ztcl / 100) *
                    par()$usr[2]),
            y = c((tick.values.minor[i] - z.central.global) *
                    min(ellipse[,1]),
                  (tick.values.minor[i] - z.central.global) *
                    min(ellipse[,1])),
            col = layout$abanico$colour$ztck)
    }

    ## plot major z-ticks
    for(i in 1:length(tick.values.major)) {
      lines(x = c(par()$usr[2],
                  (1 + 0.015 * cex * layout$abanico$dimension$ztcl / 100) *
                    par()$usr[2]),
            y = c((tick.values.major[i] - z.central.global) *
                    min(ellipse[,1]),
                  (tick.values.major[i] - z.central.global) *
                    min(ellipse[,1])),
            col = layout$abanico$colour$ztck)
    }

    ## plot z-axes
    lines(ellipse, col = layout$abanico$colour$border)
    lines(rep(par()$usr[2], nrow(ellipse)), ellipse[,2],
          col = layout$abanico$colour$ztck)

    ## plot z-axis text
    text(x = (1 + 0.04 * cex * layout$abanico$dimension$ztcl / 100) *
           par()$usr[2],
         y = (tick.values.major - z.central.global) * min(ellipse[,1]),
         labels = label.z.text,
         adj = 0,
         family = layout$abanico$font.type$ztck,
         font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                        layout$abanico$font.deco$ztck],
         cex = cex * layout$abanico$font.size$ztck/12)

    ## plot z-label
    mtext(text = zlab,
          at = mean(x = c(min(ellipse[,2]),
                          max(ellipse[,2])),
                    na.rm = TRUE),
          #        at = 0, ## BUG from version 0.4.0, maybe removed in future
          side = 4,
          las = 3,
          adj = 0.5,
          line = 5 * layout$abanico$dimension$zlab.line / 100,
          col = layout$abanico$colour$zlab,
          family = layout$abanico$font.type$zlab,
          font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                         layout$abanico$font.deco$zlab],
          cex = cex * layout$abanico$font.size$zlab/12)

    ## plot values and optionally error bars
    if(error.bars == TRUE) {
      for(i in 1:length(data)) {
        arrows(x0 = arrow.coords[[i]][,1],
               x1 = arrow.coords[[i]][,2],
               y0 = arrow.coords[[i]][,3],
               y1 = arrow.coords[[i]][,4],
               length = 0,
               angle = 90,
               code = 3,
               col = value.bar[i])
      }
    }

    for(i in 1:length(data)) {
      points(data[[i]][,6][data[[i]][,6] <= limits.x[2]],
             data[[i]][,8][data[[i]][,6] <= limits.x[2]],
             col = value.dot[i],
             pch = pch[i],
             cex = layout$abanico$dimension$pch / 100)
    }

    ## calculate KDE width
    KDE.max <- 0
    for(i in 1:length(data)) {
      KDE.max <- ifelse(KDE.max < max(KDE[[i]][,2]), max(KDE[[i]][,2]), KDE.max)
    }
    KDE.scale <- (par()$usr[2] - xy.0[1]) / (KDE.max * 1.05)

    ## optionally add KDE plot
    if(kde == TRUE) {

      ## plot KDE lines
      for(i in 1:length(data)) {
        polygon(x = xy.0[1] + KDE[[i]][,2] * KDE.scale,
                y = (KDE[[i]][,1] - z.central.global) * min(ellipse[,1]),
                col = kde.fill[i],
                border = kde.line[i],
                lwd = 1.7)
      }

      ## plot KDE x-axis
      axis(side = 1,
           at = c(xy.0[1], par()$usr[2]),
           col = layout$abanico$colour$xtck3,
           col.axis = layout$abanico$colour$xtck3,
           labels = NA,
           tcl = -layout$abanico$dimension$xtcl3 / 200,
           cex = cex)

      axis(side = 1,
           at = c(xy.0[1], par()$usr[2]),
           labels = as.character(round(c(0, KDE.max.plot), 3)),
           line = 2 * layout$abanico$dimension$xtck3.line / 100 - 2,
           lwd = 0,
           col = layout$abanico$colour$xtck3,
           family = layout$abanico$font.type$xtck3,
           font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                          layout$abanico$font.deco$xtck3],
           col.axis = layout$abanico$colour$xtck3,
           cex.axis = layout$abanico$font.size$xtck3/12)

      mtext(text = paste(xlab[3],
                         " (bw ",
                         round(x = KDE.bw,
                               digits = 3),
                         ")",
                         sep = ""),
            at = (xy.0[1] + par()$usr[2]) / 2,
            side = 1,
            line = 2.5 * layout$abanico$dimension$xlab3.line / 100,
            col = layout$abanico$colour$xlab3,
            family = layout$abanico$font.type$xlab3,
            font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                           layout$abanico$font.deco$xlab3],
            cex = cex * layout$abanico$font.size$xlab3/12)
    }

    ## optionally add histogram or dot plot axis
    if(hist == TRUE) {
      axis(side = 1,
           at = c(xy.0[1], par()$usr[2]),
           labels = as.character(c(0, hist.max.plot)),
           line = -1 * layout$abanico$dimension$xtck3.line / 100 - 2,
           lwd = 0,
           col = layout$abanico$colour$xtck3,
           family = layout$abanico$font.type$xtck3,
           font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                          layout$abanico$font.deco$xtck3],
           col.axis = layout$abanico$colour$xtck3,
           cex.axis = layout$abanico$font.size$xtck3/12)

      ## add label
      mtext(text = "n",
            at = (xy.0[1] + par()$usr[2]) / 2,
            side = 1,
            line = -3.5 * layout$abanico$dimension$xlab2.line / 100,
            col = layout$abanico$colour$xlab2,
            family = layout$abanico$font.type$xlab2,
            font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                           layout$abanico$font.deco$xlab2],
            cex = cex * layout$abanico$font.size$xlab2/12)

      ## plot ticks
      axis(side = 1,
           at = c(xy.0[1], par()$usr[2]),
           col = layout$abanico$colour$xtck2,
           col.axis = layout$abanico$colour$xtck2,
           labels = NA,
           tcl = layout$abanico$dimension$xtcl2 / 200,
           cex = cex)

      ## calculate scaling factor for histogram bar heights
      hist.scale <- (par()$usr[2] - xy.0[1]) / (KDE.max.plot * 1.05)

      ## draw each bar for each data set
      for(i in 1:length(data)) {
        for(j in 1:length(hist.data[[i]]$density)) {
          ## calculate x-coordinates
          hist.x.i <- c(xy.0[1],
                        xy.0[1],
                        xy.0[1] + hist.data[[i]]$density[j] * hist.scale,
                        xy.0[1] + hist.data[[i]]$density[j] * hist.scale)

          ## calculate y-coordinates
          hist.y.i <- c((hist.data[[i]]$breaks[j] - z.central.global) *
                          min(ellipse[,1]),
                        (hist.data[[i]]$breaks[j + 1] - z.central.global) *
                          min(ellipse[,1]),
                        (hist.data[[i]]$breaks[j + 1] - z.central.global) *
                          min(ellipse[,1]),
                        (hist.data[[i]]$breaks[j] - z.central.global) *
                          min(ellipse[,1]))

          ## remove data out of z-axis range
          hist.y.i <- ifelse(hist.y.i < min(ellipse[,2]),
                             min(ellipse[,2]),
                             hist.y.i)
          hist.y.i <- ifelse(hist.y.i > max(ellipse[,2]),
                             max(ellipse[,2]),
                             hist.y.i)

          ## draw the bars
          polygon(x = hist.x.i,
                  y = hist.y.i,
                  col = kde.fill[i],
                  border = kde.line[i])
        }
      }
    }

    ## optionally add dot plot
    if(dots == TRUE) {
      for(i in 1:length(data)) {
        for(j in 1:length(hist.data[[i]]$counts)) {

          ## calculate scaling factor for histogram bar heights
          dots.distance <- (par()$usr[2] - (xy.0[1] + par()$cxy[1] * 0.4)) / hist.max.plot

          dots.x.i <- seq(from = xy.0[1] + par()$cxy[1] * 0.4,
                          by = dots.distance,
                          length.out = hist.data[[i]]$counts[j])

          dots.y.i <- rep((hist.data[[i]]$mids[j] - z.central.global) *
                            min(ellipse[,1]), length(dots.x.i))

          ## remove data out of z-axis range
          dots.x.i <- dots.x.i[dots.y.i >= min(ellipse[,2]) &
                                 dots.y.i <= max(ellipse[,2])]
          dots.y.i <- dots.y.i[dots.y.i >= min(ellipse[,2]) &
                                 dots.y.i <= max(ellipse[,2])]

          if(max(c(0, dots.x.i), na.rm = TRUE) >= (par()$usr[2] -
                                                     par()$cxy[1] * 0.4)) {
            dots.y.i <- dots.y.i[dots.x.i < (par()$usr[2] - par()$cxy[1] * 0.4)]
            dots.x.i <- dots.x.i[dots.x.i < (par()$usr[2] - par()$cxy[1] * 0.4)]
            pch.dots <- c(rep(20, length(dots.x.i) - 1), 15)
          } else {
            pch.dots <- rep(20, length(dots.x.i))
          }

          ## plot points
          points(x = dots.x.i,
                 y = dots.y.i,
                 pch = "|",
                 cex = 0.7 * cex,
                 col = kde.line[i])

        }
      }
    }

    ## optionally add stats, i.e. min, max, median sample text
    if(length(stats) > 0) {
      text(x = stats.data[,1],
           y = stats.data[,2],
           pos = 2,
           labels = round(stats.data[,3], 1),
           family = layout$abanico$font.type$stats,
           font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                          layout$abanico$font.deco$stats],
           cex = cex * layout$abanico$font.size$stats/12,
           col = layout$abanico$colour$stats)
    }

    ## optionally add rug
    if(rug == TRUE) {
      for(i in 1:length(rug.coords)) {
        lines(x = rug.coords[[i]][1,],
              y = rug.coords[[i]][2,],
              col = value.rug[data.global[i,10]])
      }
    }

    ## plot KDE base line
    lines(x = c(xy.0[1], xy.0[1]),
          y = c(min(ellipse[,2]), max(ellipse[,2])),
          col = layout$abanico$colour$border)

    ## draw border around plot
    if(frame == 1) {
      polygon(x = c(limits.x[1], min(ellipse[,1]), par()$usr[2],
                    par()$usr[2], min(ellipse[,1])),
              y = c(0, max(ellipse[,2]), max(ellipse[,2]),
                    min(ellipse[,2]), min(ellipse[,2])),
              border = layout$abanico$colour$border,
              lwd = 0.8)
    } else if(frame == 2) {
      polygon(x = c(limits.x[1], min(ellipse[,1]), par()$usr[2],
                    par()$usr[2], min(ellipse[,1]), limits.x[1]),
              y = c(2, max(ellipse[,2]), max(ellipse[,2]),
                    min(ellipse[,2]), min(ellipse[,2]), -2),
              border = layout$abanico$colour$border,
              lwd = 0.8)
    } else if(frame == 3) {
      polygon(x = c(limits.x[1], par()$usr[2],
                    par()$usr[2], limits.x[1]),
              y = c(max(ellipse[,2]), max(ellipse[,2]),
                    min(ellipse[,2]), min(ellipse[,2])),
              border = layout$abanico$colour$border,
              lwd = 0.8)
    }

    ## optionally add legend content
    if(missing(legend) == FALSE) {
      ## store and change font familiy
      par.family <- par()$family
      par(family = layout$abanico$font.type$legend)

      legend(x = legend.pos[1],
             y = 0.8 * legend.pos[2],
             xjust = legend.adj[1],
             yjust = legend.adj[2],
             legend = legend,
             pch = pch,
             col = value.dot,
             text.col = value.dot,
             text.font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                                 layout$abanico$font.deco$legend],
             cex = cex * layout$abanico$font.size$legend/12,
             bty = "n")

      ## restore font family
      par(family = par.family)
    }

    ## optionally add subheader text
    mtext(text = mtext,
          side = 3,
          line = (shift.lines - 2) * layout$abanico$dimension$mtext / 100,
          col = layout$abanico$colour$mtext,
          family = layout$abanico$font.type$mtext,
          font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                         layout$abanico$font.deco$mtext],
          cex = cex * layout$abanico$font.size$mtext / 12)

    ## add summary content
    for(i in 1:length(data)) {
      if(summary.pos[1] != "sub") {
        text(x = summary.pos[1],
             y = summary.pos[2],
             adj = summary.adj,
             labels = label.text[[i]],
             col = summary.col[i],
             family = layout$abanico$font.type$summary,
             font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                            layout$abanico$font.deco$summary],
             cex = cex * layout$abanico$font.size$summary / 12)
      } else {
        if(mtext == "") {
          mtext(side = 3,
                line = (shift.lines- 1 - i) *
                  layout$abanico$dimension$summary / 100 ,
                text = label.text[[i]],
                col = summary.col[i],
                family = layout$abanico$font.type$summary,
                font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                               layout$abanico$font.deco$summary],
                cex = cex * layout$abanico$font.size$summary / 12)
        }
      }
    }
  } else {
    ## setup plot area
    par(mar = c(4, 4, shift.lines + 5, 4),
        xpd = TRUE,
        cex = cex)

    if(layout$abanico$dimension$figure.width != "auto" |
         layout$abanico$dimension$figure.height != "auto") {
      par(mai = layout$abanico$dimension$margin / 25.4,
          pin = c(layout$abanico$dimension$figure.width / 25.4 -
                    layout$abanico$dimension$margin[2] / 25.4 -
                    layout$abanico$dimension$margin[4] / 25.4,
                  layout$abanico$dimension$figure.height / 25.4 -
                    layout$abanico$dimension$margin[1] / 25.4 -
                    layout$abanico$dimension$margin[3]/25.4))
    }

    ## create empty plot
    par(new = TRUE)
    plot(NA,
         xlim = limits.y,
         ylim = c(limits.x[1], limits.x[2] * (1 / plot.ratio)),
         main = "",
         sub = sub,
         xlab = "",
         ylab = "",
         xaxs = "i",
         yaxs = "i",
         frame.plot = FALSE,
         axes = FALSE)

    ## add y-axis label
    mtext(text = ylab,
          at = 0,
          adj = 0.5,
          side = 1,
          line = 3 * layout$abanico$dimension$ylab.line / 100,
          col = layout$abanico$colour$ylab,
          family = layout$abanico$font.type$ylab,
          font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                         layout$abanico$font.deco$ylab],
          cex = cex * layout$abanico$font.size$ylab/12)

    ## calculate upper x-axis label values
    label.x.upper <- if(log.z == TRUE) {
      as.character(round(1/axTicks(side = 2)[-1] * 100, 1))
    } else {
      as.character(round(1/axTicks(side = 2)[-1], 1))
    }

#     ## optionally, plot 2-sigma-bar
#     if(bar.fill[1] != "none") {
#
#       if(is.numeric(centrality) == TRUE & length(centrality) > length(data)) {
#         for(i in 1:length(centrality)) {
#           polygon(x = bars[i,1:4],
#                   y = bars[i,5:8],
#                   col = bar.fill[i],
#                   border = bar.line[i])
#         }
#       } else {
#         for(i in 1:length(data)) {
#           polygon(y = bars[i,1:4],
#                   x = bars[i,5:8],
#                   col = bar.fill[i],
#                   border = bar.line[i])
#         }
#       }
#     }

    ## remove unwanted parts
    polygon(y = c(par()$usr[2],
                  par()$usr[2],
                  par()$usr[2] * 2,
                  par()$usr[2] * 2),
            x = c(min(ellipse[,2]) * 2,
                  max(ellipse[,2]) * 2,
                  max(ellipse[,2]) * 2,
                  min(ellipse[,2]) * 2),
            col = bg.original,
            lty = 0)

    ## optionally, plot dispersion polygon
    if(polygon.fill[1] != "none") {
      for(i in 1:length(data)) {
        polygon(x = polygons[i,8:14],
                y = polygons[i,1:7],
                col = polygon.fill[i],
                border = polygon.line[i])
      }
    }

    ## optionally, add minor grid lines
    if(grid.minor != "none") {
      for(i in 1:length(tick.values.minor)) {
        lines(y = c(limits.x[1], min(ellipse[,1])),
              x = c(0, (tick.values.minor[i] - z.central.global) * min(ellipse[,1])),
              col = grid.minor,
              lwd = 1)
      }
      for(i in 1:length(tick.values.minor)) {
        lines(y = c(xy.0[2], par()$usr[2]),
              x = c((tick.values.minor[i] - z.central.global) * min(ellipse[,1]),
                    (tick.values.minor[i] - z.central.global) * min(ellipse[,1])),
              col = grid.minor,
              lwd = 1)
      }
    }

    ## optionally, add major grid lines
    if(grid.major != "none") {
      for(i in 1:length(tick.values.major)) {
        lines(y = c(limits.x[1], min(ellipse[,2])),
              x = c(0, (tick.values.major[i] - z.central.global) * min(ellipse[,2])),
              col = grid.major,
              lwd = 1)
      }
      for(i in 1:length(tick.values.major)) {
        lines(y = c(xy.0[2],y.max),
              x = c((tick.values.major[i] - z.central.global) * min(ellipse[,2]),
                    (tick.values.major[i] - z.central.global) * min(ellipse[,2])),
              col = grid.major,
              lwd = 1)
      }
    }

#     ## optionally, plot central value lines
#     if(lwd[1] > 0 & lty[1] > 0) {
#       for(i in 1:length(data)) {
#         x2 <- r / sqrt(1 + f^2 * (
#           data[[i]][1,5] - z.central.global)^2)
#         y2 <- (data[[i]][1,5] - z.central.global) * x2
#         lines(y = c(limits.x[1], x2, xy.0[2],y.max),
#               x = c(0, y2, y2, y2),
#               lty = lty[i],
#               lwd = lwd[i],
#               col = centrality.col[i])
#       }
#     }

    ## optionally add further lines
    if(missing(line) == FALSE) {
      for(i in 1:length(line)) {
        lines(y = line.coords[[i]][1,1:3],
              x = line.coords[[i]][2,1:3],
              col = line.col[i])
        text(y = line.coords[[i]][1,3],
             x = line.coords[[i]][2,3] + par()$cxy[2] * 0.3,
             labels = line.label[i],
             pos = 2,
             col = line.col[i],
             cex = cex * 0.9)
      }
    }

    ## add plot title
    cex.old <- par()$cex
    par(cex = layout$abanico$font.size$main / 12)
    title(main = main,
          family = layout$abanico$font.type$main,
          font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                         layout$abanico$font.deco$main],
          col.main = layout$abanico$colour$main,
          line = (shift.lines + 3.5) * layout$abanico$dimension$main / 100)
    par(cex = cex.old)

    ## calculate lower x-axis (precision)
    x.axis.ticks <- axTicks(side = 2)
    x.axis.ticks <- x.axis.ticks[c(TRUE, x.axis.ticks <= limits.x[2])]
    x.axis.ticks <- x.axis.ticks[x.axis.ticks <= max(ellipse[,2])]

    ## x-axis with lables and ticks
    axis(side = 2,
         at = x.axis.ticks,
         col = layout$abanico$colour$xtck1,
         col.axis = layout$abanico$colour$xtck1,
         labels = NA,
         tcl = -layout$abanico$dimension$xtcl1 / 200,
         cex = cex)
    axis(side = 2,
         at = x.axis.ticks,
         line = 2 * layout$abanico$dimension$xtck1.line / 100 - 2,
         lwd = 0,
         col = layout$abanico$colour$xtck1,
         family = layout$abanico$font.type$xtck1,
         font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                        layout$abanico$font.deco$xtck1],
         col.axis = layout$abanico$colour$xtck1,
         cex.axis = layout$abanico$font.size$xlab1/12)

    ## extend axis line to right side of the plot
    lines(y = c(max(x.axis.ticks), max(ellipse[,2])),
          x = c(limits.y[1], limits.y[1]),
          col = layout$abanico$colour$xtck1)

    ## draw closing tick on right hand side
    axis(side = 2,
         tcl = -layout$abanico$dimension$xtcl1 / 200,
         lwd = 0,
         lwd.ticks = 1,
         at = limits.x[2],
         labels = FALSE,
         col = layout$abanico$colour$xtck1)

    axis(side = 2,
         tcl = layout$abanico$dimension$xtcl2 / 200,
         lwd = 0,
         lwd.ticks = 1,
         at = limits.x[2],
         labels = FALSE,
         col = layout$abanico$colour$xtck2)

    ## add lower axis label
    mtext(xlab[2],
          at = (limits.x[1] + max(ellipse[,2])) / 2,
          side = 2,
          line = 2.5 * layout$abanico$dimension$xlab1.line / 100,
          col = layout$abanico$colour$xlab1,
          family = layout$abanico$font.type$xlab1,
          font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                         layout$abanico$font.deco$xlab1],
          cex = cex * layout$abanico$font.size$xlab1/12)

    ## add upper axis label
    mtext(xlab[1],
          at = (limits.x[1] + max(ellipse[,2])) / 2,
          side = 2,
          line = -3.5 * layout$abanico$dimension$xlab2.line / 100,
          col = layout$abanico$colour$xlab2,
          family = layout$abanico$font.type$xlab2,
          font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                         layout$abanico$font.deco$xlab2],
          cex = cex * layout$abanico$font.size$xlab2/12)

    ## plot upper x-axis
    axis(side = 2,
         at = x.axis.ticks[-1],
         col = layout$abanico$colour$xtck2,
         col.axis = layout$abanico$colour$xtck2,
         labels = NA,
         tcl = layout$abanico$dimension$xtcl2 / 200,
         cex = cex)

    ## remove first tick label (infinity)
    label.x.upper <- label.x.upper[1:(length(x.axis.ticks) - 1)]

    axis(side = 2,
         at = x.axis.ticks[-1],
         labels = label.x.upper,
         line = -1 * layout$abanico$dimension$xtck2.line / 100 - 2,
         lwd = 0,
         col = layout$abanico$colour$xtck2,
         family = layout$abanico$font.type$xtck2,
         font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                        layout$abanico$font.deco$xtck2],
         col.axis = layout$abanico$colour$xtck2,
         cex.axis = layout$abanico$font.size$xlab2/12)

    ## plot y-axis
    if(y.axis == TRUE) {
      char.height <- par()$cxy[2]
      tick.space <- axisTicks(usr = limits.y, log = FALSE)
      tick.space <- (max(tick.space) - min(tick.space)) / length(tick.space)
      if(tick.space < char.height * 1.7) {
        axis(side = 1,
             tcl = -layout$abanico$dimension$ytcl / 200,
             lwd = 1,
             lwd.ticks = 1,
             at = c(-2, 2),
             labels = c("", ""),
             las = 1,
             col = layout$abanico$colour$ytck)

        axis(side = 1,
             at = 0,
             tcl = 0,
             line = 2 * layout$abanico$dimension$ytck.line / 100 - 2,
             labels = paste("\u00B1", "2"),
             las = 1,
             family = layout$abanico$font.type$ytck,
             font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                            layout$abanico$font.deco$ytck],
             col.axis = layout$abanico$colour$ytck,
             cex.axis = layout$abanico$font.size$ylab/12)
      } else {
        axis(side = 1,
             at = seq(-2, 2, by = 2),
             col = layout$abanico$colour$ytck,
             col.axis = layout$abanico$colour$ytck,
             labels = NA,
             las = 1,
             tcl = -layout$abanico$dimension$ytcl / 200,
             cex = cex)
        axis(side = 1,
             at = seq(-2, 2, by = 2),
             line = 2 * layout$abanico$dimension$ytck.line / 100 - 2,
             lwd = 0,
             las = 1,
             col = layout$abanico$colour$ytck,
             family = layout$abanico$font.type$ytck,
             font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                            layout$abanico$font.deco$ytck],
             col.axis = layout$abanico$colour$ytck,
             cex.axis = layout$abanico$font.size$ylab/12)
      }
    } else {
      axis(side = 1,
           at = 0,
           col = layout$abanico$colour$ytck,
           col.axis = layout$abanico$colour$ytck,
           labels = NA,
           las = 1,
           tcl = -layout$abanico$dimension$ytcl / 200,
           cex = cex)
      axis(side = 1,
           at = 0,
           line = 2 * layout$abanico$dimension$ytck.line / 100 - 2,
           lwd = 0,
           las = 1,
           col = layout$abanico$colour$ytck,
           family = layout$abanico$font.type$ytck,
           font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                          layout$abanico$font.deco$ytck],
           col.axis = layout$abanico$colour$ytck,
           cex.axis = layout$abanico$font.size$ylab/12)
    }

    ## plot minor z-ticks
    for(i in 1:length(tick.values.minor)) {
      lines(y = c(par()$usr[4],
                  (1 + 0.015 * cex * layout$abanico$dimension$ztcl / 100) *
                    y.max),
            x = c((tick.values.minor[i] - z.central.global) *
                    min(ellipse[,2]),
                  (tick.values.minor[i] - z.central.global) *
                    min(ellipse[,2])),
            col = layout$abanico$colour$ztck)
    }

    ## plot major z-ticks
    for(i in 1:length(tick.values.major)) {
      lines(y = c(par()$usr[4],
                  (1 + 0.03 * cex * layout$abanico$dimension$ztcl / 100) *
                    y.max),
            x = c((tick.values.major[i] - z.central.global) *
                    min(ellipse[,2]),
                  (tick.values.major[i] - z.central.global) *
                    min(ellipse[,2])),
            col = layout$abanico$colour$ztck)
    }

    ## plot z-axes
    lines(ellipse, col = layout$abanico$colour$border)
    lines(y = rep(par()$usr[4], nrow(ellipse)),
          x = ellipse[,1],
          col = layout$abanico$colour$ztck)

    ## plot z-axis text
    text(y = (1 + 0.06 * cex * layout$abanico$dimension$ztcl / 100) *
           y.max,
         x = (tick.values.major - z.central.global) * min(ellipse[,2]),
         labels = label.z.text,
         adj = 0.5,
         family = layout$abanico$font.type$ztck,
         font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                        layout$abanico$font.deco$ztck],
         cex = cex * layout$abanico$font.size$ztck/12)

    ## plot z-label
    mtext(text = zlab,
          at = 0,
          side = 3,
          las = 1,
          adj = 0.5,
          line = 2.5 * layout$abanico$dimension$zlab.line / 100,
          col = layout$abanico$colour$zlab,
          family = layout$abanico$font.type$zlab,
          font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                         layout$abanico$font.deco$zlab],
          cex = cex * layout$abanico$font.size$zlab/12)

    ## plot values and optionally error bars
    if(error.bars == TRUE) {
      for(i in 1:length(data)) {
        arrows(y0 = arrow.coords[[i]][,1],
               y1 = arrow.coords[[i]][,2],
               x0 = arrow.coords[[i]][,3],
               x1 = arrow.coords[[i]][,4],
               length = 0,
               angle = 90,
               code = 3,
               col = value.bar[i])
      }
    }

    for(i in 1:length(data)) {
      points(y = data[[i]][,6][data[[i]][,6] <= limits.x[2]],
             x = data[[i]][,8][data[[i]][,6] <= limits.x[2]],
             col = value.dot[i],
             pch = pch[i],
             cex = layout$abanico$dimension$pch / 100)
    }

    ## calculate KDE width
    KDE.max <- 0
    for(i in 1:length(data)) {
      KDE.max <- ifelse(KDE.max < max(KDE[[i]][,2]), max(KDE[[i]][,2]), KDE.max)
    }
    KDE.scale <- (par()$usr[4] - xy.0[2]) / (KDE.max * 1.05)


    ## optionally add KDE plot
    if(kde == TRUE) {

      ## plot KDE lines
      for(i in 1:length(data)) {
        polygon(y = xy.0[2] + KDE[[i]][,2] * KDE.scale,
                x = (KDE[[i]][,1] - z.central.global) * min(ellipse[,2]),
                col = kde.fill[i],
                border = kde.line[i],
                lwd = 1.7)
      }

      ## plot KDE x-axis
      axis(side = 2,
           at = c(xy.0[2], y.max),
           col = layout$abanico$colour$xtck3,
           col.axis = layout$abanico$colour$xtck3,
           labels = NA,
           tcl = -layout$abanico$dimension$xtcl3 / 200,
           cex = cex)

      axis(side = 2,
           at = c(xy.0[2], y.max),
           labels = as.character(round(c(0, KDE.max.plot), 3)),
           line = 2 * layout$abanico$dimension$xtck3.line / 100 - 2,
           lwd = 0,
           col = layout$abanico$colour$xtck3,
           family = layout$abanico$font.type$xtck3,
           font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                          layout$abanico$font.deco$xtck3],
           col.axis = layout$abanico$colour$xtck3,
           cex.axis = layout$abanico$font.size$xtck3/12)

      mtext(text = paste(xlab[3],
                         " (bw ",
                         round(x = KDE.bw,
                               digits = 3),
                         ")",
                         sep = ""),
            at = (xy.0[2] + y.max) / 2,
            side = 2,
            line = 2.5 * layout$abanico$dimension$xlab3.line / 100,
            col = layout$abanico$colour$xlab3,
            family = layout$abanico$font.type$xlab3,
            font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                           layout$abanico$font.deco$xlab3],
            cex = cex * layout$abanico$font.size$xlab3/12)
    }

    ## optionally add histogram or dot plot axis
    if(hist == TRUE) {
      axis(side = 2,
           at = c(xy.0[2], y.max),
           labels = as.character(c(0, hist.max.plot)),
           line = -1 * layout$abanico$dimension$xtck3.line / 100 - 2,
           lwd = 0,
           col = layout$abanico$colour$xtck3,
           family = layout$abanico$font.type$xtck3,
           font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                          layout$abanico$font.deco$xtck3],
           col.axis = layout$abanico$colour$xtck3,
           cex.axis = layout$abanico$font.size$xtck3/12)

      ## add label
      mtext(text = "n",
            at = (xy.0[2] + y.max) / 2,
            side = 2,
            line = -3.5 * layout$abanico$dimension$xlab2.line / 100,
            col = layout$abanico$colour$xlab2,
            family = layout$abanico$font.type$xlab2,
            font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                           layout$abanico$font.deco$xlab2],
            cex = cex * layout$abanico$font.size$xlab2/12)

      ## plot ticks
      axis(side = 2,
           at = c(xy.0[2], y.max),
           col = layout$abanico$colour$xtck2,
           col.axis = layout$abanico$colour$xtck2,
           labels = NA,
           tcl = layout$abanico$dimension$xtcl2 / 200,
           cex = cex)

      ## calculate scaling factor for histogram bar heights
      hist.scale <- (par()$usr[4] - xy.0[2]) / (KDE.max.plot * 1.05)

      ## draw each bar for each data set
      for(i in 1:length(data)) {
        for(j in 1:length(hist.data[[i]]$density)) {
          ## calculate x-coordinates
          hist.x.i <- c(xy.0[2],
                        xy.0[2],
                        xy.0[2] + hist.data[[i]]$density[j] * hist.scale,
                        xy.0[2] + hist.data[[i]]$density[j] * hist.scale)

          ## calculate y-coordinates
          hist.y.i <- c((hist.data[[i]]$breaks[j] - z.central.global) *
                          min(ellipse[,2]),
                        (hist.data[[i]]$breaks[j + 1] - z.central.global) *
                          min(ellipse[,2]),
                        (hist.data[[i]]$breaks[j + 1] - z.central.global) *
                          min(ellipse[,2]),
                        (hist.data[[i]]$breaks[j] - z.central.global) *
                          min(ellipse[,2]))

          ## remove data out of z-axis range
          hist.y.i <- ifelse(hist.y.i < min(ellipse[,1]),
                             min(ellipse[,1]),
                             hist.y.i)
          hist.y.i <- ifelse(hist.y.i > max(ellipse[,1]),
                             max(ellipse[,1]),
                             hist.y.i)

          ## draw the bars
          polygon(y = hist.x.i,
                  x = hist.y.i,
                  col = kde.fill[i],
                  border = kde.line[i])
        }
      }
    }

    ## optionally add dot plot
    if(dots == TRUE) {
      for(i in 1:length(data)) {
        for(j in 1:length(hist.data[[i]]$counts)) {

          ## calculate scaling factor for histogram bar heights
          dots.distance <- (par()$usr[4] - (xy.0[2] + par()$cxy[1] * 0.4)) / hist.max.plot

          dots.x.i <- seq(from = xy.0[2] + par()$cxy[2] * 0.4,
                          by = dots.distance,
                          length.out = hist.data[[i]]$counts[j])

          dots.y.i <- rep((hist.data[[i]]$mids[j] - z.central.global) *
                            min(ellipse[,2]), length(dots.x.i))

          ## remove data out of z-axis range
          dots.x.i <- dots.x.i[dots.y.i >= min(ellipse[,1]) &
                                 dots.y.i <= max(ellipse[,1])]
          dots.y.i <- dots.y.i[dots.y.i >= min(ellipse[,1]) &
                                 dots.y.i <= max(ellipse[,1])]

          if(max(c(0, dots.x.i), na.rm = TRUE) >= (par()$usr[4] -
                                                     par()$cxy[2] * 0.4)) {
            dots.y.i <- dots.y.i[dots.x.i < (par()$usr[4] - par()$cxy[2] * 0.4)]
            dots.x.i <- dots.x.i[dots.x.i < (par()$usr[4] - par()$cxy[2] * 0.4)]
            pch.dots <- c(rep(20, length(dots.x.i) - 1), 15)
          } else {
            pch.dots <- rep(20, length(dots.x.i))
          }

          ## plot points
          points(y = dots.x.i,
                 x = dots.y.i,
                 pch = "-",
                 cex = 0.7 * cex,
                 col = kde.line[i])
        }
      }
    }

    ## optionally add stats, i.e. min, max, median sample text
    if(length(stats) > 0) {
      text(y = stats.data[,1],
           x = stats.data[,2],
           pos = 2,
           labels = round(stats.data[,3], 1),
           family = layout$abanico$font.type$stats,
           font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                          layout$abanico$font.deco$stats],
           cex = cex * layout$abanico$font.size$stats/12,
           col = layout$abanico$colour$stats)
    }

    ## optionally add rug
    if(rug == TRUE) {
      for(i in 1:length(rug.coords)) {
        lines(y = rug.coords[[i]][1,],
              x = rug.coords[[i]][2,],
              col = value.rug[data.global[i,10]])
      }
    }

    ## plot KDE base line
    lines(y = c(xy.0[2], xy.0[2]),
          x = c(min(ellipse[,1]), max(ellipse[,1])),
          col = layout$abanico$colour$border)

    ## draw border around plot
    polygon(y = c(limits.x[1], min(ellipse[,2]), y.max,
                  y.max, min(ellipse[,2])),
            x = c(0, max(ellipse[,1]), max(ellipse[,1]),
                  min(ellipse[,1]), min(ellipse[,1])),
            border = layout$abanico$colour$border,
            lwd = 0.8)

    ## optionally add legend content
    if(missing(legend) == FALSE) {
      ## store and change font familiy
      par.family <- par()$family
      par(family = layout$abanico$font.type$legend)

      legend(y = legend.pos[2],
             x = 0.8 * legend.pos[1],
             xjust = legend.adj[2],
             yjust = legend.adj[1],
             legend = legend,
             pch = pch,
             col = value.dot,
             text.col = value.dot,
             text.font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                                 layout$abanico$font.deco$legend],
             cex = cex * layout$abanico$font.size$legend/12,
             bty = "n")

      ## restore font family
      par(family = par.family)
    }

    ## optionally add subheader text
    mtext(text = mtext,
          side = 3,
          line = (shift.lines - 2 + 3.5) * layout$abanico$dimension$mtext / 100,
          col = layout$abanico$colour$mtext,
          family = layout$abanico$font.type$mtext,
          font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                         layout$abanico$font.deco$mtext],
          cex = cex * layout$abanico$font.size$mtext / 12)

    ## add summary content
    for(i in 1:length(data)) {
      if(summary.pos[1] != "sub") {
        text(x = summary.pos[1],
             y = summary.pos[2],
             adj = summary.adj,
             labels = label.text[[i]],
             col = summary.col[i],
             family = layout$abanico$font.type$summary,
             font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                            layout$abanico$font.deco$summary],
             cex = cex * layout$abanico$font.size$summary / 12)
      } else {
        if(mtext == "") {
          mtext(side = 3,
                line = (shift.lines - 1 + 3.5 - i) *
                  layout$abanico$dimension$summary / 100 ,
                text = label.text[[i]],
                col = summary.col[i],
                family = layout$abanico$font.type$summary,
                font = (1:4)[c("plain", "bold", "italic", "bold italic") ==
                               layout$abanico$font.deco$summary],
                cex = cex * layout$abanico$font.size$summary / 12)
        }
      }
    }
  }

  ##sTeve
  if(fun){sTeve()}


  ## create and resturn numeric output
  if(output == TRUE) {
    return(list(xlim = limits.x,
                ylim = limits.y,
                zlim = limits.z,
                polar.box = c(limits.x[1],
                              limits.x[2],
                              min(ellipse[,2]),
                              max(ellipse[,2])),
                cartesian.box = c(xy.0[1],
                                  par()$usr[2],
                                  xy.0[2],
                                  max(ellipse[,2])),
                plot.ratio = plot.ratio,
                data = data,
                data.global = data.global,
                KDE = KDE))
  }
}
