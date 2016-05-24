#' Elemental Graphic for Display of Dependent Sample Data
#'
#' Plots dependent sample data beginning from a scatterplot for the X,Y pairs;
#' proceeds to display difference scores as point projections; also X and Y
#' means, as well as the mean of the difference scores.
#'
#' Paired X and Y values are plotted as scatterplot. The identity reference line
#' (for Y = X) is drawn. Parallel projections of data points to (a lower-left)
#' line segment show how each point relates to its X-Y = D difference;
#' semitransparent "shadow" points are used to display the distribution of
#' difference scores, with thin grey lines leading from each raw datapoint to
#' its shadow projection on the difference distribution. The range of that
#' difference score distribution is drawn as a blue line beneath the shadow
#' points and the mean difference is displayed as a heavy dashed purple line,
#' parallel to the identity reference line. Means for X and Y are also plotted
#' (as thin dashed vertical and horizontal lines), and rug plots are shown for
#' the distributions of X (at the top of graphic) and Y (on the right side). The
#' 95\% confidence interval for the population mean difference is also shown
#' graphically as a green band, perpendicular to the mean treatment effect line.
#' Because all data points are plotted relative to the identity line, and
#' summary results are shown graphically, clusters, data trends, outliers, and
#' possible uses of transformations are readily seen, possibly to be
#' accommodated.
#'
#' In summary, the graphic shows all initial data points relative to the
#' identity line, adds projections (to the 'north' and 'east') showing the
#' marginal distributions of X and Y, as well as projections to the 'southwest'
#' where the difference scores for each point are drawn. Means for all three
#' distributions are shown using straight lines; the confidence interval for the
#' population mean difference score is also shown. Summary statistics are
#' printed as side effects of running the function for the dependent sample
#' analysis.
#'
#' @param data is an n X 2 dataframe or matrix. First column defines X
#'   (intially for horzontal axis), the second defines Y.
#' @param main optional main title (as character); can be supplied by user. The default value is
#'   \code{"default_granova_title"}, which will print a generic title for the graphic.
#' @param revc reverses X,Y specifications
#' @param xlab optional label (as character) for horizontal axis. If not
#'   defined, axis labels are taken from colnames of data.
#' @param ylab optional label (as character) for vertical axis. If not
#'   defined, axis labels are taken from colnames of data.
#' @param conf.level The confidence level at which to perform a dependent sample t-test.
#'   Defaults to \code{0.95} (95\% Confidence)
#' @param plot.theme argument indicating a ggplot2 theme to apply to the
#'   graphic; defaults to a customized theme created for the dependent sample graphic
#' @param northeast.padding (numeric) extends axes toward lower left,
#'   effectively moving data points to the southwest. Defaults to zero padding.
#' @param southwest.padding (numeric) extends axes toward upper right,
#'   effectively moving data points to the southwest. Defaults to zero padding.
#'   Making both southwest and northeast padding smaller moves points farther apart,
#'   while making both larger moves data points closer together.
#' @param ... Optional arguments to/from other functions
#' @return Returns a plot object of class \code{ggplot}.
#'
#' @author Brian A. Danielak \email{brian@@briandk.com}\cr
#'   Robert M. Pruzek \email{RMPruzek@@yahoo.com}
#'
#' with contributions by:\cr
#'   William E. J. Doane \email{wil@@drdoane.com}\cr
#'   James E. Helmreich \email{James.Helmreich@@Marist.edu}\cr
#'   Jason Bryer \email{jason@@bryer.org}
#'
#' @seealso \code{\link{granovagg.1w}},
#'   \code{\link{granovagg.ds}}, \code{\link{granovaGG}}
#'
#' @example demo/granovagg.ds.R
#' @import stats
#' @import utils
#' @export
#' @references Pruzek, R. M., & Helmreich, J. E. (2009). Enhancing Dependent Sample Analyses with Graphics. Journal of Statistics Education, 17(1), 21.
#' @references Wickham, H. (2009). Ggplot2: Elegant Graphics for Data Analysis. New York: Springer.
#' @references Wilkinson, L. (1999). The Grammar of Graphics. Statistics and computing. New York: Springer.
granovagg.ds <- function(data       = NULL,
                         revc       = FALSE,
                         main       = "default_granova_title",
                         xlab       = NULL,
                         ylab       = NULL,
                         conf.level = 0.95,
                         plot.theme = "theme_granova_ds",
                         northeast.padding = 0,
                         southwest.padding = 0,
                         ...
                )

{

  GetData <- function(data) {
    data <- CheckData(data)
    data <- ReverseXAndY(data)
    data <- EnsureDataHasColumnNames(data)
    data <- EnsureDataIsADataFrame(data)
    return(data)
  }

  FormatDataForPlotting <- function(dsp) {
    return(
      data.frame(
        x_values = dsp$data[, 2],
        y_values = dsp$data[, 1]
      )
    )
  }

  CheckData <- function(data) {
    IsDataNull(data)
    IsDataInTwoColumnFormat(data)
    return(data)
  }

  ReverseXAndY <- function(data) {
    output <- data
    if (revc) {
      output[, 1] <- data[, 2]
      output[, 2] <- data[, 1]
      names(output) <- names(data)[2:1]
    }
    return(output)
  }

  IsDataNull <- function(data) {
    if (is.null(length(data))) {
      stop("It looks like you didn't pass any data to granovagg.ds")
    }
  }

  IsDataInTwoColumnFormat <- function(data) {
    message <- "It looks like the data you handed in isn't in two-column (n x 2) format. granovagg.ds needs n x 2 data to work."
    if (is.null(dim(data))) {
      stop(message)
    }

    if (dim(data)[2] != 2) {
      stop(message)
    }
  }

  EnsureDataIsADataFrame <- function(data) {
    output <- data
    if (!is.data.frame(data)) {
      output <- as.data.frame(output)
    }
    return(output)
  }

  EnsureDataHasColumnNames <- function(data) {
    output <- data
    if (is.null(colnames(data))) {
      colnames(data) <- c("x", "y")
    }

    return(data)
  }

  GetXs <- function(data) {
    return(data[, 1])
  }

  GetYs <- function(data) {
    return(data[, 2])
  }

  GetEffect <- function(dsp) {
    return(GetXs(dsp$plotting_data) - GetYs(dsp$plotting_data))
  }

  GetTtest <- function(data, conf.level) {
    return(t.test(data[, 1],
                  data[, 2],
                  paired     = TRUE,
                  conf.level = conf.level
                 )
          )
  }

  GetStats <- function(dsp, conf.level) {
    ttest <- GetTtest(dsp$data, conf.level)
    return(data.frame(lower.treatment.effect = as.numeric(ttest$conf.int[1]),
                      mean.treatment.effect  = as.numeric(ttest$estimate[1]),
                      upper.treatment.effect = as.numeric(ttest$conf.int[2]),
                      t.statistic            = as.numeric(ttest$statistic[1])
                     )
          )
  }

  GetGraphicsParams <- function(dsp) {
    .aggregate.data.range  <- c(
      range(dsp$plotting_data$x_values),
      range(dsp$plotting_data$y_values)
    )
    .extrema               <- c(max(.aggregate.data.range), min(.aggregate.data.range))
    .square.data.range     <- max(.extrema) - min(.extrema)
    .southwest.padding     <- (65/100) * .square.data.range
    .northeast.padding     <- (15/100) * .square.data.range
    .lower.graphical.bound <- min(.extrema) - .southwest.padding
    .upper.graphical.bound <- max(.extrema) + .northeast.padding
    .bounds                <- c(.lower.graphical.bound, .upper.graphical.bound)
    .center                <- mean(.bounds)
    .crossbow.anchor       <- mean(.bounds) + min(.bounds)
    .shadow.offset         <- (1/100)*.square.data.range

    return(list(square.data.range = .square.data.range,
                bounds            = .bounds,
                shadow.offset     = .shadow.offset,
                anchor            = .crossbow.anchor,
                point.size        = I(2.5),
                mean.line.size    = I(1/2)
               )
          )
  }

  GetShadows <- function(dsp) {
    x.shadow <- (dsp$effect / 2) +
                (3 * dsp$params$bounds[1] + dsp$params$bounds[2]) / 4 +
                (4 * dsp$params$shadow.offset)
    y.shadow <- x.shadow - dsp$effect
    return(data.frame(x.shadow, y.shadow))
  }

  GetCrossbow <- function(dsp) {
    return(data.frame(x     = min(dsp$shadows$x.shadow) - (2 * dsp$params$shadow.offset),
                      y     = max(dsp$shadows$y.shadow) - (2 * dsp$params$shadow.offset),
                      x.end = max(dsp$shadows$x.shadow) - (2 * dsp$params$shadow.offset),
                      y.end = min(dsp$shadows$y.shadow) - (2 * dsp$params$shadow.offset)
                     )
          )
  }

  GetCIBand <- function(dsp) {
    return(data.frame(x.end = ((dsp$params$anchor - dsp$stats$lower.treatment.effect) / 2) - (3 * (dsp$params$shadow.offset)),
                      y.end = ((dsp$params$anchor + dsp$stats$lower.treatment.effect) / 2) - (3 * (dsp$params$shadow.offset)),
                      x     = ((dsp$params$anchor - dsp$stats$upper.treatment.effect) / 2) - (3 * (dsp$params$shadow.offset)),
                      y     = ((dsp$params$anchor + dsp$stats$upper.treatment.effect) / 2) - (3 * (dsp$params$shadow.offset)),
                      color = factor(paste(100 * conf.level, "% CI", " (t = ", round(dsp$stats$t.statistic, digits = 2), ")", sep =""))
                     )
          )

  }

  GetTreatmentLine <- function(dsp) {
    return(data.frame(intercept = dsp$stats$mean.treatment.effect,
                      slope     = 1,
                      color     = factor(paste("Mean Diff. =", round(dsp$stats$mean.treatment.effect, digits = 2)))
                     )
          )
  }

  GetCrossElementCoordinates <- function(dsp) {
    color <- NULL # to appease R CMD check
    crossbow <- dsp$crossbow
    ci.band  <- subset(dsp$CIBand, select = -color)
    output   <- rbind(crossbow, ci.band)
    return(output)
  }

  EnsureCrossElementsAppearInVisualBounds <- function(dsp) {
    minimum <- min(dsp$params$bounds,
                   dsp$cross.elements$x,
                   dsp$cross.elements$y.end
               )
    return(c(minimum - dsp$params$shadow.offset, max(dsp$params$bounds)))
  }

  GetTrails <- function(dsp) {
    return(data.frame(x.trail.start = dsp$plotting_data$x_values,
                      y.trail.start = dsp$plotting_data$y_values,
                      x.trail.end   = GetXs(dsp$shadow),
                      y.trail.end   = GetYs(dsp$shadow)
                     )
          )
  }

  GetColors <- function(dsp) {
    return(list(treatment.line = "#542570",
                rugplot        = "black",
                mean.line      = "#542570",
                CIBand         = "#33A02C",
                crossbow       = "#377EB8"
               )
          )
  }

  PrintSummary <- function(dsp) {
    summary <- GetPrintedSummary(dsp)
    summary <- RenamePrintedSummaryRows(summary)
    summary <- round(summary, digits = 3)

    print(summary)
  }

  GetPrintedSummary <- function(dsp) {
    n <- dim(dsp$data)[1]
    mean.1 <- mean(dsp$data[, 1])
    mean.2 <- mean(dsp$data[, 2])
    mean.d <- mean.1 - mean.2
    standard.deviation.d <- sd(dsp$data[, 1] - dsp$data[, 2])
    effect.size <- (mean.d / standard.deviation.d)
    r.xy <- cor(dsp$data[, 1], dsp$data[, 2])
    r.x.plus.y.d <- cor((dsp$data[, 1] + dsp$data[, 2]), (dsp$data[, 1] - dsp$data[, 2]))
    lower.treatment.confidence <- dsp$stats$upper.treatment.effect
    upper.treatment.confidence <- dsp$stats$lower.treatment.effect
    t.value <- dsp$t.test$statistic
    degrees.of.freedom <- dsp$t.test$parameter
    p.value <- dsp$t.test$p.value

    return(matrix(c(n,
                    mean.1,
                    mean.2,
                    mean.d,
                    standard.deviation.d,
                    effect.size,
                    r.xy,
                    r.x.plus.y.d,
                    lower.treatment.confidence,
                    upper.treatment.confidence,
                    t.value,
                    degrees.of.freedom,
                    p.value
                  ),
                 ncol = 1
           )
    )
  }

  RenamePrintedSummaryRows <- function(summary) {
    dimnames(summary) <- list(
                           c("n",
                             paste(colnames(dsp$data)[1], "mean"),
                             paste(colnames(dsp$data)[2], "mean"),
                             paste("mean(D = ", colnames(dsp$data)[1], " - ", colnames(dsp$data)[2], ")",  sep = ""),
                             paste("SD(D)"),
                             paste("Effect Size"),
                             paste("r(", colnames(dsp$data)[1], ", ", colnames(dsp$data)[2], ")", sep = ""),
                             paste("r(", colnames(dsp$data)[1], " + ", colnames(dsp$data)[2], ", D)", sep = ""),
                             paste("Lower ", (100 * conf.level), "% ", "Confidence Interval", sep = ""),
                             paste("Upper ", (100 * conf.level), "% ", "Confidence Interval", sep = ""),
                             paste("t (D-bar)"),
                             paste("df.t"),
                             paste("p-value (t-statistic)")
                           ), "Summary Statistics")

    return(summary)
  }

  dsp                <- list(data = GetData(data))
  dsp$plotting_data  <- FormatDataForPlotting(dsp)
  dsp$effect         <- GetEffect(dsp)
  dsp$stats          <- GetStats(dsp, conf.level)
  dsp$t.test         <- GetTtest(dsp$data, conf.level)
  dsp$params         <- GetGraphicsParams(dsp)
  dsp$shadows        <- GetShadows(dsp)
  dsp$crossbow       <- GetCrossbow(dsp)
  dsp$CIBand         <- GetCIBand(dsp)
  dsp$cross.elements <- GetCrossElementCoordinates(dsp)
  dsp$params$bounds  <- EnsureCrossElementsAppearInVisualBounds(dsp)
  dsp$treatment.line <- GetTreatmentLine(dsp)
  dsp$trails         <- GetTrails(dsp)
  dsp$colors         <- GetColors(dsp)
  PrintSummary(dsp)



  # Because of the way ggplot2 creates plot objects, layers can be
  # added to a plot p simply by calling "p <- p + newLayer"


  InitializeGgplot <- function(dsp) {
    return(
      ggplot(
        aes_string(
          x = "x_values",
          y = "y_values"
        ),
        data = dsp$plotting_data
      )
    )
  }

  TreatmentLine <- function(dsp) {
    return(
      geom_abline(
        aes_string(
          intercept = "intercept",
          slope     = "slope",
          color     = "color"
        ),
        alpha    = 0.5,
        size     = I(1),
        linetype = "dashed",
        data     = dsp$treatment.line
      )
    )
  }

  RawData <- function(dsp) {
    return(
      geom_point(
        aes_string(
          x = "x_values",
          y = "y_values"
        ),
        data = dsp$plotting_data,
        size = dsp$params$point.size
      )
    )
  }

  IdentityLine <- function() {
    return(
      geom_abline(
        slope     = 1,
        intercept = 0,
        alpha     = 0.75,
        size      = 1
      )
    )
  }

  ScaleX <- function(dsp) {
    return(
      scale_x_continuous(
        limits = dsp$params$bounds
      )
    )
  }

  ScaleY <- function(dsp) {
    return(
      scale_y_continuous(
        limits = dsp$params$bounds
      )
    )
  }

  PadViewingWindow <- function(params) {
    ne.offset = params$square.data.range * southwest.padding
    sw.offset = params$square.data.range * northeast.padding
    padded.window = c(params$bounds[1] - sw.offset, params$bounds[2] + ne.offset)

    return(
      coord_cartesian(
        xlim = padded.window,
        ylim = padded.window
      )
    )
  }

  RugPlot <- function(dsp) {
    return(
      geom_rug(
        size  = 1/2,
        alpha = 1/3,
        color = dsp$colors$rugplot,
        sides = "tr", # top and right sides
        data  = dsp$plotting_data
      )
    )
  }

  XMeanLine <- function(dsp) {
    return(
      geom_vline(
        xintercept = mean(dsp$data[, 2]),
        color      = dsp$colors$mean.line,
        size       = dsp$params$mean.line.size,
        linetype   = "dashed",
        alpha      = I(1/2)
      )
    )
  }

  YMeanLine <- function(dsp)  {
    return(
      geom_hline(
        yintercept = mean(dsp$data[, 1]),
        color      = dsp$colors$mean.line,
        size       = dsp$params$mean.line.size,
        linetype   = "dashed",
        alpha      = I(1/2)
      )
    )
  }

  Crossbow <- function(dsp) {
    return(
      geom_segment(
        aes_string(x    = "x",
                   y    = "y",
                   xend = "x.end",
                   yend = "y.end"
        ),
        size  = 3/4,
        alpha = 3/4,
        color = dsp$colors$crossbow,
        data  = dsp$crossbow
      )
    )
  }

  CIBand <- function(dsp) {
    return(
      geom_segment(
        aes_string(
          x     = "x",
          y     = "y",
          xend  = "x.end",
          yend  = "y.end",
          color = "color"
        ),
        size = 2,
        data = dsp$CIBand
      )
    )
  }

  Shadows <- function(dsp) {
    return(
      geom_point(
        aes_string(
          x = "x.shadow",
          y = "y.shadow"
        ),
        data  = dsp$shadow,
        size  = dsp$params$point.size,
        shape = 16,
        alpha = 1/2
      )
    )
  }

  Trails <- function(dsp) {
    return(
      geom_segment(
        aes_string(
          x    = "x.trail.start",
          y    = "y.trail.start",
          xend = "x.trail.end",
          yend = "y.trail.end"
        ),
        data     = dsp$trails,
        size     = 1/3,
        color    = "black",
        linetype = 1,
        alpha    = 1/10
      )
    )
  }

  ColorScale <- function(dsp) {
    colors <- c(dsp$colors$treatment.line, dsp$colors$CIBand)

    return(scale_color_manual(values = colors, name = ""))
  }

  XLabel <- function(dsp) {
    result <- colnames(dsp$data)[1]
    if(!is.null(xlab)) {
      result <- xlab
    }

    return(xlab(result))
  }

  YLabel <- function(dsp) {
    result <- colnames(dsp$data)[2]
    if(!is.null(ylab)) {
      result <- ylab
    }

    return(ylab(result))
  }

  Title <- function(main) {
    output.title <- "Dependent Sample Assessment Plot"
    if (main != "default_granova_title") {
      output.title <- main
    }

    return(
      ggtitle(output.title)
    )

  }

  ForceCoordinateAxesToBeEqual <- function() {
    return(coord_fixed(ratio = 1))
  }


  p <- InitializeGgplot(dsp)
  p <- p + TreatmentLine(dsp)
  p <- p + XMeanLine(dsp) + YMeanLine(dsp)
  p <- p + Shadows(dsp)
  p <- p + Trails(dsp)
  p <- p + RawData(dsp)
  # p <- p + Theme(plot.theme)
  p <- p + IdentityLine()
  p <- p + RugPlot(dsp)
  p <- p + Crossbow(dsp)
  p <- p + CIBand(dsp)
  p <- p + ColorScale(dsp)
  p <- p + ScaleX(dsp) + ScaleY(dsp)
  p <- p + PadViewingWindow(dsp$params)
  p <- p + ForceCoordinateAxesToBeEqual()
  p <- p + Title(main)
  p <- p + XLabel(dsp)
  p <- p + YLabel(dsp)

  return(p)

}
