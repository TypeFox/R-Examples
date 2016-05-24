#' Elemental Graphic Display for Contrast Effect of ANOVA
#'
#' Provides graphic displays that shows data and effects for a priori contrasts
#' in ANOVA contexts; also corresponding numerical results.
#'
#' Function provides graphic displays of contrast effects for prespecified
#' contrasts in ANOVA. Data points are displayed as relevant for each contrast
#' based on comparing groups according to the positive and negative contrast
#' coefficients for each contrast on the horizontal axis, against response
#' values on the vertical axis. Data points corresponding to groups not being
#' compared in any contrast (coefficients of zero) are ignored. For each
#' contrast (generally as part of a 2 x 2 panel) a line segment is given that
#' compares the (weighted) mean of the response variable for the negative
#' coefficients versus the positive coefficients. Standardized contrasts are
#' used, wherein the sum of (magnitudes) of negative coefficients is unity; and
#' the same for positive coefficients. If a line is `notably' different from
#' horizontal (i.e. slope of zero), a `notable' effect has been identified;
#' however, the question of statistical significance generally depends on a
#' sound context-based estimate of standard error for the corresponding effect.
#' This means that while summary aov numerical results and test statistics are
#' presented (see below), the appropriateness of the default standard error
#' generally requires the analyst's judgment. The response values are to be
#' input in (a stacked) form, i.e. as a vector, for all cells (cf. arg. ylab).
#' The matrix of contrast vectors \code{contrasts} must have G rows (the number
#' of groups), and a number of columns equal to the number of prespecified
#' contrasts, at most G-1. If the number of columns of \code{contrasts} is G-1,
#' then the number per group, or cell size, is taken to be
#' \code{length(data)/G}, where \code{G = nrow(contrasts)}.
#'
#' If the number of columns of \code{contrasts} is less than G-1 then the user
#' must stipulate \code{npg}, the number in each group or cell.  The function
#' is designed for the case when all cell sizes are the same, and may be most
#' helpful when the a priori contrasts are mutually orthogonal (e.g., in power
#' of 2 designs, or their fractional counterparts; also when specific row or
#' column comparisons, or their interactions (see the example below based on
#' rat weight gain data)). It is not essential that contrasts be mutually
#' orthogonal; but mutual linear independence is required. (When factor levels
#' correspond to some underlying continuum a standard application might use
#' \code{con = contr.poly(G)}, for G the number of groups; consider also
#' \code{contr.helmert(G)}.)  The final plot in each application shows the data
#' for all groups or cells in the design, where groups are simply numbered from
#' 1:G, for G the number of groups, on the horizontal axis, versus the response
#' values on the vertical axis.
#'
#' @param data Vector of scores for all equally sized groups, or a data.fame or
#'   matrix where each column represents a group.
#' @param contrasts Matrix of column contrasts with dimensions (number of
#'   groups [G]) x (number of contrasts) [generally (G x G-1)].
#' @param ylab Character; y axis label. Defaults to a generic granova title.
#' @param plot.theme argument indicating a ggplot2 theme to apply to the
#'   graphic; defaults to a customized theme created for the contrast graphic
#' @param print.four.plots.per.page If \code{TRUE}, the function lays out four plots per page and sends
#'   each page to the graphics device. When running R interactively, you'll have an opportunity to review each page
#'   before seeing the next page. Also, when \code{print.four.plots.per.page} is \code{TRUE}, the function won't
#'   return any plot objects as output. When \code{print.four.plots.per.page} is set to \code{FALSE},
#'   the function returns a list of ggplot objects, one element per plot.
#' @param jj Numeric; controls \code{\link{jitter}} and allows you to control the
#'   degree of jitter in the contrast plots. \code{jj} is divided by 100 and passed as the \code{amount}
#'   parameter to \code{\link{jitter}}.
#' @param ... Optional arguments to/from other functions.
#' @return If \code{print.four.plots.per.page} is set to \code{FALSE}, the function returns
#'   a list of ggplot objects, one element per plot. That allows you to access any individual plot
#'   or plots, then modify them as you wish (with ggplot2 commands, for example).
#'   When \code{print.four.plots.per.page} is set to \code{TRUE}
#'   (the default), the function prints four plots per page on a graphical device
#'   but returns \code{NULL}.
#'
#'   The function also provides printed output:
#'   \item{Weighted Means}{Table showing the (weighted) means for positive
#'      and negative coefficients for each (row) contrast, and for each row, the
#'      difference between these means, and the standardized effect size in the
#'      final column.}
#'   \item{summary.lm}{Summary results for a linear
#'      model analysis based on the R function \code{lm} (When effects are simple,
#'      as in an equal n's power of 2 design, mean differences will generally
#'      correspond to the linear regression coefficients as seen in the \code{lm}
#'      summary results.)}
#'   \item{Contrasts}{The contrast matrix you specified.}
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
#' @keywords hplot
#' @example demo/granovagg.contr.R
#' @references Wickham, H. (2009). Ggplot2: Elegant Graphics for Data Analysis. New York: Springer.
#' @references Wilkinson, L. (1999). The Grammar of Graphics. Statistics and computing. New York: Springer.
#' @import gridExtra
#' @import reshape2
#' @import stats
#' @import utils
#' @export
granovagg.contr <- function(data,
                            contrasts,
                            ylab       = "default_y_label",
                            plot.theme = "theme_granova_contr",
                            print.four.plots.per.page = TRUE,
                            jj = 1,
                            ...
                   )
{

  # Plots responses by contrasts.
  # 'data' must be vector of scores for all equal size groups.
  # 'con' must be matrix of column contrasts with dimensions (number of groups) x (number of contrasts)
  # [generally n X n-1].  The number of rows = number 'cells' or groups.
  # Basic lm (regression) results are provided; orthogonal contrasts are ideal (but not essential).
  # 'jj' controls jitter.

  # 'ctr' is shorthand for the ConTRast data object that will hold all the information for plotting
  FormatResponseData <- function(data) {
    is.data.one.dimensional <- is.null(dim(data)[2])
    if (is.data.one.dimensional) {
      return(data)
    }

    return(stack(as.data.frame(data))[, 1])
  }

  GetDegreeOfJitter <- function(jj) {
    result <- jj / 100
    return(result)
  }

  std.contr <- function(contrasts, tolerance = sqrt(.Machine$double.eps)^0.6) {
      if (!is.matrix(contrasts)) {
          contrasts <- as.matrix(contrasts)
      }
      if (sum(abs(colMeans(contrasts))) > tolerance) {
          stop("Input vector/matrix must have mean zero (for each column)")
      }
      if (ncol(contrasts) == 1) {
          contrasts <- matrix(contrasts, ncol = 1)
      }
      dg <- apply(abs(contrasts), 2, sum)
      if (length(dg) == 1) {
          dg <- as.matrix(dg)
      }
      standardized.contrasts <- round(2 * contrasts %*% diag(1/dg), 3)

      return(standardized.contrasts)
  }

  indic <- function(xx) {
             mm <- matrix(0, length(xx), length(unique(xx)))
             indx <- ifelse(xx == col(mm), 1, 0)
             indx
  }

  AdaptVariablesFromGranovaComputations <- function () {

    response            <- FormatResponseData(data)
    contrasts           <- as.matrix(contrasts)
    number.of.groups    <- nrow(contrasts)
    responses.per.group <- length(response)/number.of.groups
    group.identifiers   <- rep(1:number.of.groups, ea = responses.per.group)
    indicator.matrix    <- indic(group.identifiers)
    indicated.contrasts <- indicator.matrix %*% contrasts
    standardized.contrasts <- std.contr(indicated.contrasts)

    return(
        list(
          response                      = response,
          contrast.matrix               = contrasts,
          scaled.standardized.contrasts = standardized.contrasts * responses.per.group,
          number.of.contrasts           = dim(standardized.contrasts)[2],
          number.of.groups              = number.of.groups,
          responses.per.group           = responses.per.group
        )
    )
  }

  GetContrastPlotData <- function (ctr) {
    return(
      lapply(
        X         = 1:ctr$number.of.contrasts,
        FUN       = ExtractDataForPlot,
        contrasts = ctr$scaled.standardized.contrasts,
        response  = ctr$response
      )
    )
  }

  GetLinearModel <- function(ctr) {
    Contrast <- ctr$scaled.standardized.contrasts
    Response <- ctr$response

    return(lm(Response ~ Contrast))
  }

  ExtractDataForPlot <- function (contrasts, response, index) {
      non.zero.indicators <- contrasts[, index] != 0
      x.values <- contrasts[, index][non.zero.indicators]
      y.values <- response[non.zero.indicators]
      raw.data <- data.frame(x.values, y.values)

    return(
        list(
          raw.data     = raw.data,
          summary.data = GetSummary(raw.data)
        )
    )
  }

  GetSummary <- function(data) {
    # To appease R CMD check
    x.values <- NULL
    y.values <- NULL

    return(
      ddply(data, .(x.values > 0), summarise,
        contrasts          = mean(x.values),
        responses          = mean(y.values)
      )
    )
  }

  GetContrastPlots <- function (ctr) {
    return(
      lapply(
        X             = 1:ctr$number.of.contrasts,
        FUN           = ComposeContrastPlot,
        plot.data     = ctr$contrast.plot.data
      )
    )
  }

  ComposeContrastPlot <- function(plot.data, index) {
    p <- ggplot()
    p <- p + MeanResponse(plot.data[[index]]$raw.data$y.values)
    p <- p + JitteredResponsesByContrast(plot.data[[index]]$raw.data)
    p <- p + EffectsOfContrasts(plot.data[[index]]$summary.data)
    p <- p + ConnectEffectMeans(plot.data[[index]]$summary.data)
    p <- p + Theme(plot.theme)
    p <- p + ContrastPlotTitle(ctr, index)
    p <- p + ContrastPlotXLabel(ctr, index)
    p <- p + YLabel()

    return(p)
  }

  MeanResponse <- function(response) {
    return(
      geom_hline(
        aes_string(yintercept = mean(response)),
        color = brewer.pal(8, "Set1")[1],
        data  = as.data.frame(data),
        alpha = 0.5,
        size  = 0.3
      )
    )
  }

  JitteredResponsesByContrast <- function (data) {
    return(
      geom_point(
        aes_string(
          x = "x.values",
          y = "y.values"
        ),
        data     = data,
        position = position_jitter(height = 0, width = GetDegreeOfJitter(jj))
      )
    )
  }

  EffectsOfContrasts <- function(data) {
    return(
      geom_point(
        aes_string(
          x = "contrasts",
          y = "responses"
        ),
        data  = data,
        color = brewer.pal(8, "Set1")[2],
        size  = 3,
        alpha = 0.75
      )
    )
  }

  ConnectEffectMeans <- function(data) {
    return(
      geom_line(
        aes_string(
          x = "contrasts",
          y = "responses"
        ),
        data  = data,
        color = brewer.pal(8, "Set1")[2],
        alpha = 1
      )
    )
  }

  ContrastPlotTitle <- function(ctr, index) {
    plot.title <- paste("Coefficients vs. Response\n", GetContrastName(ctr$contrast.matrix, index))
    return(
        ggtitle(plot.title)
    )
  }

  ContrastPlotXLabel <- function(ctr, index) {
    return(
        xlab(paste(GetContrastName(ctr$contrast.matrix, index)))
    )
  }

  YLabel <- function() {
    result <- ylab
    if (ylab == "default_y_label") {
      result <- "Outcome (Response)"
    }

    return(ylab(paste(result)))
  }

  GetSummaryPlotData <- function(ctr) {
    raw.data <- as.data.frame(
                   matrix(ctr$response, ncol = ctr$number.of.groups)
                 )
    raw.data <- RenameSummaryColumnNames(raw.data)
    raw.data <- melt(raw.data)
    raw.data$variable <- as.numeric(raw.data$variable)
    summary.data <- GetGroupSummary(raw.data)

    return(list(
                raw.data     = raw.data,
                summary.data = summary.data
          )
    )
  }

  RenameSummaryColumnNames <- function(data) {
    colnames(data) <- sapply(
                        1:ncol(data),
                        function(index) {paste(index)}
                      )

    return(data)
  }

  GetGroupSummary <- function(data) {
    # Appease R CMD Check
    variable <- NULL
    value <- NULL
    standard.deviation <- NULL

    output <- ddply(
      data,
      .(variable),
      summarise,
      group      = unique(variable),
      group.mean = mean(value),
      standard.deviation = sd(value)
    )
    output <- transform(
      output,
      pooled.standard.deviation = mean(standard.deviation^2)^0.5
    )
    return(subset(output, select = -variable))
  }

  ComposeSummaryPlot <- function(plot.data) {
    p <- ggplot()
    p <- p + MeanResponse(plot.data$raw.data$value)
    p <- p + RawScoresByGroup(plot.data$raw.data)
    p <- p + MeansByGroup(plot.data$summary.data)
    p <- p + ConnectGroupResponseMeans(plot.data$summary.data)
    p <- p + Theme(plot.theme)
    p <- p + GroupSummaryPlotTitle(ctr)
    p <- p + GroupSummaryXLabel()
    p <- p + YLabel()
    return(p)
  }

  RawScoresByGroup <- function(data) {
    return(
      geom_point(
        aes_string(
          x = "as.factor(variable)",
          y = "value"
        ),
        data = data,
        position = position_jitter(height = 0, width = 3 * GetDegreeOfJitter(jj))
      )
    )
  }

  MeansByGroup <- function(data) {
    return(
      geom_point(
        aes_string(
          x = "group",
          y = "group.mean"
        ),
        data  = data,
        color = brewer.pal(8, "Set1")[2],
        size  = I(3),
        alpha = 1
      )
    )
  }

  ConnectGroupResponseMeans <- function(data) {
    return(
      geom_line(
        aes_string(
          x = "group",
          y = "group.mean"
        ),
        data  = data,
        color = brewer.pal(8, "Set1")[2],
        alpha = 1
      )
    )
  }

  GroupSummaryPlotTitle <- function(ctr) {
    plot.title <- paste("Responses for all groups\n", "each n = ", ctr$responses.per.group)
    return(
      ggtitle(plot.title)
    )
  }

  GroupSummaryXLabel <- function() {
    return(xlab("Group Indicator"))
  }

  CollateOutputPlots <- function(ctr) {
    output <- list(NULL)

    for (i in 1:ctr$number.of.contrasts) {
      output[[i]] <- ctr$contrast.plots[[i]]
    }
    output[[ctr$number.of.contrasts + 1]] <- ctr$summary.plot

    return(output)
  }

  OptionalPlotPrinting <- function(output) {
    if (print.four.plots.per.page) {
      LayoutFourPlotsPerPage(output)
    }
  }

  PrintOutput <- function() {
    PrintLinearModelSummary(ctr$linear.model)
    PrintSummaryDataByContrast(ctr)
    PrintSummaryDataByGroup(ctr)
    PrintContrasts(ctr)
  }

  PrintLinearModelSummary <- function(model) {
    model.summary <- summary(model)
    message("\nLinear Model Summary")
    print(model.summary)
  }

  GetSummaryDataByContrast <- function(x, pooled.standard.deviation) {
    ExtractData <- function(x) {
      summary.data <- x$summary.data
      neg <- summary.data$responses[summary.data$contrasts <= 0]
      pos <- summary.data$responses[summary.data$contrasts > 0]
      diff <- pos - neg
      stEftSze <- (pos - neg) / pooled.standard.deviation
      return(data.frame(neg, pos, diff, stEftSze))
    }
    output <- ldply(x, .fun = ExtractData)
    output <- ForceRowNamesToBeContrastNumbers(output)
    return(output)
  }

  ForceRowNamesToBeContrastNumbers <- function(x) {
    dimnames(x)[[1]] <- sapply(1:(dim(x)[[1]]), function(x) {paste("Contrast", x, sep="")})

    return(x)
  }

  PrintSummaryDataByContrast <- function(ctr) {
    message("\n(Weighted) means, mean differences, and standardized effect size")
    print(
      GetSummaryDataByContrast(ctr$contrast.plot.data, ctr$summary.plot.data$summary.data$pooled.standard.deviation[1]), digits = 3)
  }

  PrintSummaryDataByGroup <- function(ctr) {
    message("\nSummary statistics by group")
    print(ctr$summary.plot.data$summary.data, digits = 4)
  }

  PrintContrasts <- function(ctr) {
    message("\nThe contrasts you specified")
    print(ctr$contrast.matrix, digits = 3)
  }

  GetOutput <- function(ctr) {
    four.plot.message <- paste("Since you elected to print four plots per page\n",
                               "granovagg.contr won't return any plot objects.", sep = ""
                         )
    if (print.four.plots.per.page) {
      message(four.plot.message)
      LayoutFourPlotsPerPage(ctr$output)
      output <- NULL
    }

    else {
      output <- ctr$output
    }
    return(output)
  }

  ctr                        <- AdaptVariablesFromGranovaComputations()
  ctr$linear.model           <- GetLinearModel(ctr)
  ctr$contrast.plot.data     <- GetContrastPlotData(ctr)
  ctr$contrast.plots         <- GetContrastPlots(ctr)
  ctr$summary.plot.data      <- GetSummaryPlotData(ctr)
  ctr$summary.plot           <- ComposeSummaryPlot(ctr$summary.plot.data)
  ctr$output                 <- CollateOutputPlots(ctr)
  PrintOutput()

  return(GetOutput(ctr))

}
