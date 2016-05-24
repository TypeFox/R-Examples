#' @title Plot residual-based predictiveness (RBP) curve.
#'
#' @description plots the RBP curve
#' 
#' @template arg_obj
#' @param main [\code{character(1)}]\cr
#' An overall title for the plot.
#' @param xlab [\code{character(1)}]\cr
#'   Label for X-axis.
#'   Default is \dQuote{Cumulative Percentage}.
#' @param ylab [\code{character(1)}]\cr
#'   Label for Y-axis.
#'   Default is \dQuote{Estimated Residuals}.
#' @param type [\code{character(1)}]\cr
#'   The plot type that should be drawn, see \code{\link{plot}} for all possible types.
#'   Default is \code{type = "l"} for \bold{l}ines.
#' @param ylim [\code{numeric(2)}]\cr
#'   Limits for Y-axis.
#'   Default is \code{c(-1, 1.1)}.
#' @param x.adj [\code{numeric(2)}]\cr
#'   Adjustment for the X-axis.
#' @param y.adj [\code{numeric(2)}]\cr
#'   Adjustment for the Y-axis.
#' @param cond.axis [\code{logical(1)}]\cr
#'   Should an additional axis be plotted reflecting residuals conditional on y?
#'   Default is \code{FALSE}.
#' @param title.line [\code{integer(1)}]\cr
#'   Where to plot the title, see \code{\link{title}}.
#' @param add [\code{logical(1)}]\cr
#'   Should RBP plot be added to current plot?
#'   Default is \code{FALSE}.
#' @param ... [any]\cr
#'   Passed to \code{\link{plot}} or \code{\link{lines}}, depending on \code{add}.
#' @export
#' @import mlr
#' @examples
#' 
#' # Download data
#' mydata = read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
#' head(mydata)
#' 
#' # Build logit model and plot RBP curve
#' mylogit <- glm(admit ~ ., data = mydata, family = "binomial")
#' y = mydata$admit
#' pred1 = predict(mylogit, type="response")
#' obj1 = makeRBPObj(pred1, y)
#' plotRBPCurve(obj1, cond.axis = TRUE, type = "b")
#' 
#' \dontrun{
#' # Build logit model using mlr and plot RBP curve
#' task = makeClassifTask(data = mydata, target = "admit")
#' lrn = makeLearner("classif.logreg", predict.type = "prob")
#' tr = train(lrn, task)
#' pred2 = getPredictionProbabilities(predict(tr, task))
#' obj2 = makeRBPObj(pred2, y)
#' plotRBPCurve(obj2, cond.axis = TRUE, type = "b", col = 2)
#' }

plotRBPCurve = function (obj,
  main = "RBP Curve",
  xlab = "Cumulative Percentage",
  ylab = "Estimated Residuals",
  type = "l",
  ylim = c(-1, 1.2),
  x.adj = c(NA, -0.5),
  y.adj = c(NA, NA),
  cond.axis = FALSE,
  title.line = ifelse(cond.axis, 3, 2),
  add = FALSE,
  ...) {

  # argument checks
  assertClass(obj, "RBPObj")
  #assertString(main)
  #assertString(xlab)
  #assertString(ylab)
  assertString(type)
  assertNumeric(ylim, len = 2L)
  assertNumeric(x.adj, len = 2L)
  assertNumeric(y.adj, len = 2L)
  assertFlag(cond.axis)
  assertNumber(title.line)
  assertFlag(add)

  # plot or add RBP curve
  if (add) {
    lines(x = obj$axis.x, y = obj$axis.y, ...)
  } else {
    plot(x = obj$axis.x, y = obj$axis.y,
      xlab = xlab, ylab = ylab, ylim = ylim,
      main = "", type = type, yaxt = "n", xaxt = "n", ...)
    axis(1L, hadj = x.adj[1], padj = x.adj[2])
    axis(2L, las = 2L, hadj = y.adj[1], padj = y.adj[2])
    #abline(h = 0L, col = "grey")
  }

  # add conditional axis
  omp = obj$one.min.prev
  xAxis = seq(0, 1, by = 0.2)
  if (cond.axis) {
    abline(v = omp, col = "grey")
    axis(side = 1L, at = xAxis*omp, labels = xAxis,
      padj = -0.5, hadj = 0.6, pos = par()$usr[4L])
    axis(side = 3L, at = omp + xAxis*(1 - omp),
      padj = 0.5, hadj = 0.4, labels = xAxis)
  }

  title(main, line = title.line)

  return(invisible(NULL))
}

