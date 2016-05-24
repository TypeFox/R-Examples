#' Plot the result of Precision-Recall curve evaluation
#'
#' The \code{plot_eval_results} function validates Precision-Recall curves
#'    and creates a plot.
#'
#' @param object An S3 object that contains evaluation results of
#'     Precision-Recall curves.
#'
#' @param base_plot A Boolean value to specify whether the base points are
#'     plotted.
#'
#' @param ncol An integer used for the column size of multiple panes.
#'
#' @param nrow An integer used for the row size of multiple panes.
#'
#' @param ret_grob A Boolean value to specify whether the function returns a
#'     grob object.
#'
#' @param ... Not used by this function.
#'
#' @return A data frame with validation results.
#'
#' @examples
#' library(ggplot2)
#'
#' ## Plot evaluation results on test datasets r1, r2, and r3
#' testset <- create_testset("curve", c("c1", "c2", "c3"))
#' toolset <- create_toolset(set_names = "crv5")
#' eres1 <- run_evalcurve(testset, toolset)
#' autoplot(eres1)
#'
#' @rdname autoplot
#' @export
autoplot.evalcurve <- function(object, base_plot = TRUE, ret_grob = FALSE,
                               ncol = NULL, nrow = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Validate arguments
  new_args <- .validate_autoplot_evalcurve_args(object, base_plot, ret_grob,
                                                ncol, nrow, ...)

  # Create plots
  plots <- .create_plots(new_args$object)
  if (new_args$base_plot) {
    bplot <- .plot_base(object$basepoints)
    plots <- c(list(bplot), plots)
  }

  # Combine multiple plots
  ncolrow = .get_row_col(new_args$ncol, new_args$nrow, length(plots))

  ggrob <- .combine_plots(plots, ncolrow$ncol, ncolrow$nrow)

  if (new_args$ret_grob) {
    return(ggrob)
  } else {
    .plot_grob(ggrob)
  }
}

#
# Get the numbers of rows and columns for arrangeGrob
#
.get_row_col <- function (ncol, nrow, nplot) {
  if (is.null(ncol)) {
    if (nplot >= 4) {
      ncol <- 2
    } else {
      ncol <- nplot
    }
  }

  if (is.null(nrow)) {
    if (nplot >= 5) {
      nrow <- 3
    } else if (nplot == 4) {
      nrow <- 2
    } else {
      nrow <- 1
    }
  }

  list(ncol = ncol, nrow = nrow)
}

#
# Create a list of ggplot objects
#
.create_plots <- function(evalcurve) {
  preds <- evalcurve$predictions
  tscores <- evalcurve$testsum
  uniqnames <- unique(paste(tscores$toolset, tscores$toolname))

  plotfunc <- function(uname) {
    unamevec <- strsplit(uname, " ")[[1]]
    toolset <- unamevec[1]
    toolname <- unamevec[2]

    pcrows <- preds$toolset == toolset & preds$toolname == toolname
    pcurves <- preds[pcrows, ]

    tsrows <- tscores$toolset == toolset & tscores$toolname == toolname
    tscore <- tscores[tsrows, ]

    if (toolset == toolname || length(unique(tscores$toolset)) == 1) {
      tname <- toolname
    } else {
      tname <- paste(toolset, toolname, sep = ":")
    }

    .plot_curves(evalcurve$basepoints, pcurves, tscore, tname)
  }

  plots <- lapply(uniqnames, plotfunc)
}

#
# Plot base points
#
.plot_base <- function(basepoints, title = "Base points", yintercept = 0.5) {
  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_hline(yintercept = yintercept, colour = "grey",
                               linetype = 3)
  p <- p + ggplot2::geom_point(data = basepoints,
                               ggplot2::aes_string(x = "x", y = "y",
                                                   colour = "testset",
                                                   shape = "testset"),
                               size = 2)
  p <- p + ggplot2::scale_shape(solid = FALSE)
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::ggtitle(title)
  p <- p + ggplot2::xlab("recall")
  p <- p + ggplot2::ylab("precision")
  p <- p + ggplot2::theme(legend.title = ggplot2::element_blank())
  p <- p + ggplot2::theme(legend.position = "none")
  p <- p + ggplot2::scale_y_continuous(limits = c(0.0, 1.0))
  p <- p + ggplot2::coord_fixed(ratio = 1)
}

#
# Plot curves for a specified tool
#
.plot_curves <- function(basepoints, pcurves, tscore, toolname,
                         yintercept = 0.5) {

  p <- .plot_base(basepoints, toolname, yintercept)
  p <- p + ggplot2::geom_line(data = pcurves,
                              ggplot2::aes_string(x = "x", y = "y",
                                                  colour = "testset"))
  p <- p + ggplot2::geom_text(data = tscore,
                              ggplot2::aes_string(x = "lbl_pos_x",
                                                  y = "lbl_pos_y",
                                                  label = "label",
                                                  colour = "testset"))
}

#
# Plot curves for a specified tool
#
.combine_plots <- function(plots, ncol, nrow) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra needed for this function to work. Please install it.",
         call. = FALSE)
  }

  gfunc <- function(...) {gridExtra::arrangeGrob(..., ncol = ncol, nrow = nrow)}
  grob <- suppressWarnings(do.call(gfunc, plots))
}

#
# Plot grob object
#
.plot_grob <- function(grob) {
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("grid needed for this function to work. Please install it.",
         call. = FALSE)
  }
  graphics::plot.new()
  grid::grid.draw(grob)
}

#
# Validate arguments and return updated arguments
#
.validate_autoplot_evalcurve_args <- function(object, base_plot, ret_grob, ncol,
                                              nrow, ...) {

  if (!methods::is(object, "evalcurve")) {
    stop("Ivalid object type", call. = FALSE)
  }

  assertthat::assert_that(assertthat::is.flag(base_plot))
  assertthat::assert_that(assertthat::is.flag(ret_grob))

  if (!is.null(ncol) && !is.null(nrow)) {
    assertthat::assert_that(assertthat::is.number(ncol))
    assertthat::assert_that(ncol > 0)

    assertthat::assert_that(assertthat::is.number(nrow))
    assertthat::assert_that(nrow > 0)
  } else if ((!is.null(ncol) && is.null(nrow))
             || (is.null(ncol) && !is.null(nrow))) {
    stop("Both ncol and nrow must be set", call. = FALSE)
  }

  list(object = object, base_plot = base_plot, ret_grob = ret_grob, ncol = ncol,
       nrow = nrow)
}
