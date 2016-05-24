#' Mosaic plots of expected and realized population sizes
#'
#' Create mosaic plots of expected (i.e., estimated) and realized (i.e.,
#' simulated) population sizes.
#'
#' If \code{method} is \code{"split"}, the two tables of expected and realized
#' population sizes are combined into a single table, with an additional
#' conditioning variable indicating expected and realized values. A conditional
#' plot of this table is then produced using \code{\link[vcd]{cotabplot}}.
#'
#' @name spMosaic
#' @param x An object of class \code{"spTable"} created using function
#' \code{\link{spTable}}.
#' @param method A character string specifying the plot method. Possible values
#' are \code{"split"} to plot the expected population sizes on the left hand
#' side and the realized population sizes on the right hand side, and
#' \code{"color"}
#' @param \dots if \code{method} is \code{"split"}, further arguments to be
#' passed to \code{\link[vcd]{cotabplot}}.  If \code{method} is \code{"color"},
#' further arguments to be passed to \code{\link[vcd]{strucplot}}
#' @author Andreas Alfons and Bernhard Meindl
#' @seealso \code{\link{spTable}}, \code{\link[vcd]{cotabplot}},
#' \code{\link[vcd]{strucplot}}
#' @keywords hplot
#' @export
#' @examples
#' set.seed(1234)  # for reproducibility
#' data(eusilcS)   # load sample data
#' samp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize",
#'   strata="db040", weight="db090")
#' eusilcP <- simStructure(data=samp, method="direct", basicHHvars=c("age","rb090"))
#' abb <- c("B","LA","Vi","C","St","UA","Sa","T","Vo")
#' tab <- spTable(eusilcP, select=c("rb090", "db040", "hsize"))
#'
#' # expected and realized population sizes
#' spMosaic(tab, method = "split",
#'   labeling=labeling_border(abbreviate=c(db040=TRUE)))
#'
#' # realized population sizes colored according to relative
#' # differences with expected population sizes
#' spMosaic(tab, method = "color",
#'   labeling=labeling_border(abbreviate=c(db040=TRUE)))
#'
spMosaic <- function(x, method = c("split", "color"), ...) {
  if ( !class(x) == "spTable" ) {
    stop("input argument 'x' must be of class 'spTable'!\n")
  }
  method <- match.arg(method)
  # define local version of 'cotabplot'
  if ( method == "split" ) {
    ## split the plot with the sample on the left and the simulated population
    ## on the right
    # combine the tables for the sample and the simulated population
    tab <- as.table(x)
    dn <- dimnames(tab)
    dn[[length(dn)]] <- c("Sample", "Population")
    names(dn)[length(dn)] <- "Data"
    dimnames(tab) <- dn
    # define a local wrapper for cotabplot()
    localCotabplot <- function(x, ..., cond, panel) {
      cotabplot(x, cond="Data", ...)
    }
    # produce the plot
    localCotabplot(tab, ...)
  } else {
    ## plot the simulated population with cells colored according to the
    ## differences with the sample
    # adjust the expected frequencies from the sample such that the total is
    # the same as in the simulated population
    tab <- list(expected=x$expected * sum(x$realized) / sum(x$expected),
                realized=x$realized)
    # define a local wrapper for strucplot()
    localStrucplot <- function(x, ..., gp_args = list(), legend_args = list(),
                               # the following arguments are ignored
                               residuals, expected, df, condvars, shade,
                               gp, legend) {
      # compute residuals relative to sample
      residuals <- (x$realized - x$expected) / ifelse(x$expected > 0, x$expected, 1)
      # change default number of tick marks in legend
      if(is.null(legend_args$ticks)) legend_args$ticks <- 5
      # call strucplot()
      strucplot(x$realized, residuals=residuals, expected=x$expected,
                shade=TRUE, gp=spShading, gp_args=gp_args,
                legend=legend_resbased, legend_args=legend_args, ...)
    }
    # produce the plot
    localStrucplot(tab, ...)
  }
}

# shading function to color cells in strucplot()
spShading <- function(observed, residuals, expected, df = NULL, steps = 201,
                      h = NULL, c = NULL, l = NULL, power = 1.3, lty = 1,
                      eps = NULL, line_col = "black", ...) {
  ## set defaults
  if(is.null(h)) h <- c(260, 0)
  if(is.null(c)) c <- 100
  if(is.null(l)) l <- c(90, 50)

  ## get h/c/l and lty
  steps <- rep_len(steps, 1)  # number of steps for color gradient
  my.h <- rep_len(h, 2)       # positive and negative hue
  my.c <- rep_len(c, 1)       # maximum chroma
  my.l <- rep_len(l, 2)       # maximum and minimum luminance
  lty <- rep_len(lty, 2)      # positive and negative lty

  ## find range of residuals
  r <- range(residuals)

  ## obtain color information
  col.bins <- seq(r[1], r[2], length.out=steps)
  ## TODO: replace with heat_hcl: check how it looks like
  legend.col <- rev(heat_hcl(steps, h=my.h, c.=my.c, l=rev(my.l),
                                power=power))

  ## store lty information for legend
  lty.bins <- 0
  legend.lty <- lty[2:1]
  legend <- list(lty=legend.lty, lty.bins=lty.bins)

  ## set up function that computes color/lty from residuals
  rval <- function(x) {
    res <- as.vector(x)

    bin <- cut(res, breaks=col.bins, labels=FALSE, include.lowest=TRUE)
    fill <- legend.col[bin]
    dim(fill) <- dim(x)

    col <- rep_len(line_col, length(res))
    if(!is.null(eps)) {
      eps <- abs(eps)
      col[res > eps] <- legend.col[1]
      col[res < -eps] <- legend.col[length(legend.col)]
    }
    dim(col) <- dim(x)

    lty <- ifelse(x > 0, lty[1], lty[2])
    dim(lty) <- dim(x)

    return(structure(list(col = col, fill = fill, lty = lty), class = "gpar"))
  }
  attr(rval, "legend") <- legend
  attr(rval, "p.value") <- NULL
  return(rval)
}
class(spShading) <- "grapcon_generator"
