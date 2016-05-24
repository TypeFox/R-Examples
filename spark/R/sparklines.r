
#' Spark line of a numeric vector on the terminal
#'
#' For marking the minumum/maximum, the \code{crayon} package
#' is needed.
#'
#' @inheritParams spark
#' @param data The data to visualize. It can be a numeric
#'   vector, or anything that can be cut into intervals
#'   with \code{cut}. Infinite values in numeric data are ignored,
#'   and a black character is plotted instead of them.
#' @param min If not NULL, then a crayon style to mark the
#'   minimum value.
#' @param max If not NULL, then a crayon style to mark the
#'   maximum value.
#' @param ... Not used, it is an error if given.
#' @return Character scalar containing the spark line.
#'
#' @family spark
#' @method spark default
#' @export
#' @examples
#' ## Annual number of Lynx trappings
#' spark(lynx, width = "auto")
#'
#' ## Luteinizing Hormone in Blood Samples,
#' ## in blue, if the terminal supports it
#' cat(crayon::blue(spark(lh, print = FALSE)), "\n")

spark.default <- function(data, width = c("data", "auto", "screen"),
                          print = TRUE, min = NULL, max = NULL, ...) {

  if (length(list(...))) stop("Extra arguments are invalid")

  width <- match.arg(width)

  win_size <- getOption("width")
  if (width == "auto") {
    width <- if (width <= win_size) "data" else "screen"
  }

  if (width == "screen") data <- scale_to(data, win_size)

  if (is.numeric(data)) data[!is.finite(data)] <- NA

  code <- cut(data, breaks = length(spark_ticks)) %>%
    as.integer()

  res <- ifelse(is.na(code), ' ', spark_ticks[code]) %>%
    mark_max(data = data, multiplier = -1, mark = min) %>%
    mark_max(data = data, multiplier =  1, mark = max) %>%
    paste(collapse = "")

  if (print) {
    cat(res, "\n")
    invisible(res)
  } else {
    res
  }
}


mark_max <- function(data, ticks, multiplier, mark) {

  if (is.null(mark)) return(ticks)
  if (!is_installed("crayon")) {
    warning("The 'crayon' package is needed for min/max marks")
    return(ticks)
  }

  wm <- which_max(multiplier * data)
  ticks[wm] <- import("crayon", "style")(ticks[wm], as = mark)
  ticks
}
