#' Create a violin plot
#'
#' Draws a kernal densiy plot in combination with a boxplot in its middle. The shape of the violin
#' is constructed using a mirrored density curve. This plot is especially designed for cases
#' where the individual errors are zero or to small to be visualised. The idea for this plot is
#' based on the the 'volcano plot' in the ggplot2 package by Hadely Wickham and Winston Chang.
#' The general idea for the Violin Plot seems to be introduced by Hintze and Nelson (1998).
#'
#' The function is passing several arguments to the function \code{\link{plot}},
#' \code{\link[stats]{density}}, \code{\link[graphics]{boxplot}}: Supported arguments are: \code{xlim}, \code{main}, \code{xlab},
#' \code{ylab}, \code{col.violin}, \code{col.boxplot}, \code{mtext}, \code{cex}, \code{mtext}
#'
#' \bold{\code{Valid summary keywords}}\cr
#'
#' 'n', 'mean', 'median', 'sd.abs', 'sd.rel', 'se.abs', 'se.rel', 'skewness', 'kurtosis' \cr
#'
#' @param data \code{\link{numeric}} or \code{\linkS4class{RLum.Results}}
#' object (required): input data for plotting. Alternatively a \code{\link{data.frame}} or
#' a \code{\link{matrix}} can be provided, but only the first column will be considered by the
#' function
#'
#' @param boxplot \code{\link{logical}} (with default): enable or disable boxplot
#'
#' @param rug \code{\link{logical}} (with default): enable or disable rug
#'
#' @param summary \code{\link{character}} (optional): add statistic measures of
#' centrality and dispersion to the plot. Can be one or more of several
#' keywords. See details for available keywords.
#'
#' @param summary.pos \code{\link{numeric}} or \code{\link{character}} (with
#' default): optional position keywords (cf., \code{\link{legend}})
#' for the statistical summary. Alternatively, the keyword \code{"sub"} may be
#' specified to place the summary below the plot header. However, this latter
#' option in only possible if \code{mtext} is not used.
#'
#' @param na.rm \code{\link{logical}} (with default): exclude NA values
#' from the data set prior to any further operations.
#'
#' @param \dots further arguments and graphical parameters passed to
#' \code{\link{plot.default}}, \code{\link[stats]{density}} and \code{\link{boxplot}}. See details for
#' further information
#'
#' @note Although the code for this function was developed independently and just the idea for the plot
#' was based on the 'ggplot2' package plot type 'volcano', it should be mentioned that, beyond this,
#' two other R packages exist providing a possibility to produces this kind of plot, namely:
#' 'vioplot' and 'violinmplot' (see References for details).
#'
#' @section Function version: 0.1.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne (France)
#'
#' @references
#'
#' Daniel Adler (2005). vioplot: A violin plot is a combination of a box plot and a kernel density plot.
#' R package version 0.2 http://CRAN.R-project.org/package=violplot
#'
#' Hintze, J.L., Nelson, R.D., 1998. A Box Plot-Density Trace Synergism. The American Statistician 52, 181-184.
#'
#' Raphael W. Majeed (2012). violinmplot: Combination of violin plot with mean and standard deviation.
#' R package version 0.2.1. http://CRAN.R-project.org/package=violinmplot
#'
#' Wickham. H (2009). ggplot2: elegant graphics for data analysis. Springer New York.
#'
#' @seealso \code{\link[stats]{density}}, \code{\link{plot}}, \code{\link{boxplot}}, \code{\link{rug}},
#' \code{\link{calc_Statistics}}
#'
#' @examples
#' ## read example data set
#' data(ExampleData.DeValues, envir = environment())
#' ExampleData.DeValues <- Second2Gray(ExampleData.DeValues$BT998, c(0.0438,0.0019))
#'
#' ## create plot straightforward
#' plot_ViolinPlot(data = ExampleData.DeValues)
#'
#' @export
plot_ViolinPlot <- function(
  data,
  boxplot = TRUE,
  rug = TRUE,
  summary = NULL,
  summary.pos = "sub",
  na.rm = FALSE,
  ...
) {


  # Integrity tests and conversion --------------------------------------------------------------

    ##Prechecks

    if(missing(data)){
      stop("[plot_ViolinPlot()] I don't know what to do, data input needed." )

    }else{

      ##check for RLum.Results object
      if(is(data, "RLum.Results")){
        data <- get_RLum(data)

      }

      ##if data.frame or matrix
      if(is(data, "data.frame") | is(data, "matrix")){
        data <- data[,1]

      }

    }

    ##Remove NA values
    if(na.rm){
      data <- na.exclude(data)

    }

    #Further checks
    if(!is(summary.pos, "character")){
      stop("[plot_ViolinPlot()] argument 'summary.pos' needs to be of type character!")

    }

  # Pre-calculations ----------------------------------------------------------------------------

  ##density for the violin
  density <-
    density(x = data,
            bw = ifelse("bw" %in% names(list(...)),list(...)$bw,"nrd0"))

  ##some statistical parameter, get rid of the weighted statistics
  stat.summary <- suppressWarnings(calc_Statistics(as.data.frame(data), digits = 2))[[-1]]

    ##make valid summary string
    if(is.null(summary)){
      summary <- c("n","median")

    }

    ##at least show a warning for invalid keywords
    if(!all(summary %in% names(stat.summary))){
      warning(paste0("[plot_ViolinePlot()] At least one 'summary' keyword is invalid. Valid keywords are: ",
                     paste(names(stat.summary), collapse = ", ")), call. = FALSE)
    }

    ##make sure that only valid keywords make it
    summary <- summary[(summary %in% names(stat.summary))]

    stat.text <-
      paste(names(stat.summary[summary]), " = ", stat.summary[summary], collapse = " \n")

    stat.mtext <-
      paste(names(stat.summary[summary]), " = ", stat.summary[summary], collapse = " | ")





  # Plot settings -------------------------------------------------------------------------------

  ##set default values
  plot.settings <- list(
    xlim = range(density$x),
    main = "Violin Plot",
    xlab = expression(paste(D[e], "/(a.u.)")),
    ylab = "Density",
    col.violin = rgb(0,0,0,0.2),
    col.boxplot = NULL,
    mtext = ifelse(summary.pos != 'sub', "", stat.mtext),
    cex = 1
  )

  ##modify list accordingly
  plot.settings <- modifyList(plot.settings, val = list(...))


  # Plot ----------------------------------------------------------------------------------------

  ##open empty plot area
  plot(
    NA,NA,
    xlim = plot.settings$xlim,
    ylim = c(0.2,1.8),
    xlab = plot.settings$xlab,
    ylab = plot.settings$ylab,
    yaxt = "n",
    main = plot.settings$main,
    cex = plot.settings$cex
  )

  ##add polygon ... the violin
  polygon(
    x = c(density$x, rev(density$x)),
    y = c(1 + density$y / max(density$y) * 0.5,
          rev(1 - density$y / max(density$y) * 0.5)),
    col = plot.settings$col.violin,
    border = plot.settings$col.violin
  )

  ##add the boxplot
  if(boxplot){
    boxplot(
      data,
      outline = TRUE,
      boxwex = 0.4,
      horizontal = TRUE,
      axes = FALSE,
      add = TRUE,
      col = plot.settings$col.boxplot
    )

  }

  ##add rug
  if(rug){
    rug(x = data)

  }

  ##add mtext
  if(!is.null(plot.settings$mtext)){
    mtext(side = 3, text = plot.settings$mtext)

  }

  ##add stat.text
  if (summary.pos != "sub") {

    valid_keywords <-
      c(
        "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center"
      )

    if (any(
      summary.pos %in% valid_keywords
    )) {
      legend(summary.pos, legend = stat.text, bty = "n")

    }else{
      warning_text <- paste0("Value provided for 'summary.pos' is not a valid keyword, valid keywords are:",
                             paste(valid_keywords, collapse = ", "))
      warning(warning_text)

    }

  }

}
