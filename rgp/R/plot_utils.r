## plot_utils.R
##   - Utilities for plotting several types of RGP objects
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Show an overlayed plot of multiple functions
##'
##' Creates and shows and overlayed plot of one or more functions of one variable \eqn{y = f(x)}.
##'
##' @param funcs A list of functions of one variable to plot.
##' @param from The left bound of the plot, i.e. the minimum \eqn{x} value to plot.
##' @param to The right bound of the plot, i.e. the maximum \eqn{x} value to plot.
##' @param steps The number of steps, or samples, to plot.
##' @param type The plot type (e.g. l = line) as passed on to \code{\link{matplot}}.
##' @param lty The line types as passed on to \code{\link{matplot}}.
##' @param lwd The line widths as passed on to \code{\link{matplot}}.
##' @param lend The line end cap types as passed on to \code{\link{matplot}}.
##' @param pch The plot chars as passed on to \code{\link{matplot}}.
##' @param col The plot colors as passed on to \code{\link{matplot}}.
##' @param cex The character expansion sizes as passed on to \code{\link{matplot}}.
##' @param bg The background (fill) colors as passed on to \code{\link{matplot}}.
##' @param xlab The x axis label as passed on to \code{\link{matplot}}.
##' @param ylab The y axis label as passed on to \code{\link{matplot}}.
##' @param legendpos The position of the legend, passed as the \code{x} parameter to
##'   \code{\link{legend}}. 
##' @param bty The box type parameter of the legend, passed as the \code{bty} parameter to
##'   \code{\link{legend}}.
##' @param ... Graphic parameters for \code{\link{par}} and further arguments to \code{plot}.
##'   For example, use the \code{main} parameter to set a title.
##'
##' @examples
##' plotFunctions(list(function(x) sin(x),
##'                    function(x) cos(x),
##'                    function(x) 0.5*sin(2*x)+1),
##'               -pi, pi, 256)
##'
##' @export
plotFunctions <- function(funcs, from = 0, to = 1, steps = 1024,
                          type = "l", lty = 1:5, lwd = 1, lend = par("lend"),
                          pch = NULL,
                          col = 1:6, cex = NULL, bg = NA,
                          xlab = "x", ylab = "y",
                          legendpos = "bottomright", bty = "n", ...) {
  xs <- seq(from, to, length = steps)
  yslists <- Map(function(f) as.vector(Map(f, xs), mode = "numeric"), funcs)
  ys <- c(); for (ylist in yslists) ys <- c(ys, ylist)
  ysmatrix <- matrix(ys, ncol = length(funcs))

  legendBuilder <- function(f) {
    fformals <- names(formals(f))
    fbody <- exprToPlotmathExpr(body(f))
    as.expression(bquote(f(.(fformals[[1]])) == .(fbody)))
  }
  legendContent <- sapply(funcs, legendBuilder)

  matplot(xs, ysmatrix, type = type, xlab = xlab, ylab = ylab,
          lty = lty, lwd = lwd, lend = lend, pch = pch,
          col = col, cex = cex, bg = bg, ...)
  legend(legendpos, legendContent, bty = bty,
         col = col, bg = bg, lty = lty, lwd = lwd, pch = pch)
  invisible()
}

##' Plot a 2D function as a 3D surface
##'
##' Creates and shows and perspective plot of a 2D function of either the form
##' \eqn{z = f(x, y)} or \eqn{z = f(xv)}, where \eqn{xv} is a numeric of length 2.
##'
##' @param func A 2D function to plot.
##' @param lo A vector of lower limits of the plot (one entry for each dimension). 
##' @param up A vector of upper limits of the plot (one entry for each dimension). 
##' @param samples The number of samples in each dimension.
##' @param palette The color palette, use \code{NULL} to disable.
##' @param ... Graphic parameters for \code{\link{persp}}.
##' 
##' @export
plotFunction3d <- function (func = function(x) sum(x^2),
                            lo = c(0, 0), up = c(1, 1), samples = 10,
                            palette = gray.colors(256), ...) {
  x <- seq(lo[1], up[1], length = samples)
  y <- seq(lo[2], up[2], length = samples)
  z <- if (length(formals(func)) == 1) {
    fn <- function(a, b) apply(cbind(a, b), 1, func)
    outer(x, y, fn)
  } else if (length(formals(func)) == 2) {
    outer(x, y, func)
  } else {
    stop("plotFunction3d: function arity must be either 1 or 2")
  }

  if (is.null(palette)) {
    persp(x, y, z, ...)
  } else {
    colInd <- cut(z, length(palette))
    persp(x, y, z, col = palette[colInd], ...)
  }
}

##' Convert a function to an expression plottable by plotmath
##'
##' Tries to convert a function \code{func} to an expression plottable by \code{\link{plotmath}}
##' by replacing arithmetic operators and "standard" functions by plottable counterparts.
##'
##' @param func The function to convert.
##' @return An expression plottable by \code{\link{plotmath}}.
##'
##' @seealso \code{\link{funcToIgraph}}
##' @export
funcToPlotmathExpr <- function(func)
  bquote(f(.(names(formals(func)))) == .(exprToPlotmathExpr(body(func))))

##' Convert any expression to an expression that is plottable by plotmath
##'
##' Tries to convert a GP-generated expression \code{expr} to an expression plottable by
##' \code{\link{plotmath}} by replacing GP variants of arithmetic operators by their standard
##' counterparts.
##'
##' @param expr The GP-generated expression to convert.
##' @return An expression plottable by \code{\link{plotmath}}.
##'
##' @export
exprToPlotmathExpr <- function(expr)
  if (is.call(expr)) {
    expl <- as.list(expr)
    func <- expl[1]
    if (func == as.name("`/`") || func == as.name("safeDivide"))
      bquote(frac(.(exprToPlotmathExpr(expl[[2]])), .(exprToPlotmathExpr(expl[[3]]))))
    else if (func == as.name("`*`"))
      bquote(.(exprToPlotmathExpr(expl[[2]])) %.% .(exprToPlotmathExpr(expl[[3]])))
    else if (func == as.name("sqrt") || func == as.name("safeSqroot"))
      bquote(sqrt(.(exprToPlotmathExpr(expl[[2]]))))
    else if (func == as.name("loge") || func == as.name("safeLogn"))
      bquote(ln(.(exprToPlotmathExpr(expl[[2]]))))
    else
      as.call(c(func, Map(exprToPlotmathExpr, rest(expl))))
  } else expr

##' Visualization of functions and expressions as trees
##'
##' The following functions plot R expressions and functions as trees. The igraph package
##' is required for most of these functions.
##' \code{exprToGraph} transforms an R expression into a graph given as a character vector
##' of vertices V and a even-sized numeric vector of edges E. Two elements i and i+1 in E
##' encode a directed edge from V[i] to V[i+1]. 
##' \code{funcToIgraph} and \code{exprToIgraph} return an igraph graph object for an R
##' function or an R expression.
##' 
##' @param func An R function.
##' @param expr An R expression.
##' @return The result (see the details section).
##'
##' @seealso \code{\link{funcToPlotmathExpr}}
##' @rdname gpIndividualVisualization
##' @export
funcToIgraph <- function(func) exprToIgraph(body(func))

##' @rdname gpIndividualVisualization
##' @export
exprToIgraph <- function(expr) {
  if (!require("igraph")) stop("exprToIgraph: Package 'igraph' not installed.")
  exprGraph <- exprToGraph(expr)
  exprIgraph <- graph(exprGraph$edges - 1, n = length(exprGraph$vertices)) # igraph vertexes are counted from zero
  V(exprIgraph)$label <- exprGraph$vertices
  exprIgraph
}

##' @rdname gpIndividualVisualization
##' @export
exprToGraph <- function(expr) {
  vertices <- character()
  edges <- numeric()
  exprToGraphRecursive <- function(expr, cv = 1, c = 1)
  if (is.call(expr)) {
    f <- expr[[1]]
    vertices <<- c(vertices, as.character(f))
    edges <<- c(edges, cv, c)
    args <- as.list(expr)[-1]
    exprToGraphRecursive(args, cv = c, c = c + 1)
  } else if (is.list(expr) && !identical(expr, list())) {
    rc <- exprToGraphRecursive(expr[[1]], cv = cv, c = c) # first
    exprToGraphRecursive(expr[-1], cv = cv, c = rc) # rest
  } else if (is.list(expr)) {
    # empty argument list, do nothing
    c
  } else {
    vertices <<- c(vertices, as.character(expr))
    edges <<- c(edges, cv, c)
    c + 1
  }
  exprToGraphRecursive(expr)
  list(vertices = vertices,
       edges = edges[-2:-1]) # the first edge is the reflexive root edge
}

##' Plot a GP Pareto Front 
##'
##' Plots fitness/complexity/age Pareto fronts for multi-objective 
##' GP. The z-coordinate represents individual age and is shown in form
##' of a color scale, where younger individuals are bright green, individuals
##' with age \code{maxZ} are black. Individuals not on the first Pareto
##' front are shown as small gray circles, regardless of age. 
##'
##' @param x A vector of type \code{numeric} representing individual fitness.
##' @param y A vector of type \code{numeric} representing individual complexity. 
##' @param z A vector of type \code{integer} representing individual age.
##' @param indicesToMark A index vector of points to mark with red crosses.
##' @param maxZ The individual age at the large end of the age color scale.
##' @param main The plot's title.
##' @param ... Graphic parameters for \code{\link{par}} and further arguments to \code{plot}.
##'   For example, use the \code{main} parameter to set a title.
##'
##' @seealso \code{\link{funcToIgraph}}
##' @import emoa
##' @export
plotParetoFront <- function(x, y, z, indicesToMark = integer(), maxZ = 50,
                            main = sprintf("Population Pareto Front Plot (% Individuals)", length(x)), ...) {
  ranks <- nds_rank(rbind(x, y, z))
  zColorScale <- colorRamp(c("#00FF00", "#006600","#0000FF", "#000000"))
  rankOneXs <- x[ranks == 1]; rankOneYs <- y[ranks == 1]; rankOneAges <- z[ranks == 1]
  zColorScaleIndex <- pmin(rankOneAges / maxZ, 1.0)
  rankOneAgeColors <- rgb(zColorScale(zColorScaleIndex), maxColorValue = 255)

  plot(rankOneXs, rankOneYs,
       col = rankOneAgeColors, pch = 19, main = main, ...)
  points(x[ranks > 1], y[ranks > 1], col = "gray", pch = 1)
  points(x[indicesToMark], y[indicesToMark], col = "red", pch = 4)
}

