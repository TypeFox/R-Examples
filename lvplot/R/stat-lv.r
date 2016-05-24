
"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}


#' @rdname geom_lv
#' @param conf confidence level
#' @param percent numeric value: percent of data in outliers
#' @param k number of letter values shown
#' @param na.rm If \code{FALSE} (the default), removes missing values with
#'    a warning.  If \code{TRUE} silently removes missing values.
#' @section Computed/reported variables:
#' \describe{
#'   \item{k}{Number of Letter Values used for the display}
#'   \item{LV}{Name of the Letter Value}
#'   \item{width}{width of the interquartile box}
#' }
#' @export
stat_lv <- function(mapping = NULL, data = NULL, geom = "lv",
  position = "dodge", na.rm = TRUE, conf = 0.95, percent = NULL, k = NULL, show.legend = NA,
  inherit.aes = TRUE, ...)
{
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatLv,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      conf = conf,
      k = k,
      percent = percent,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname geom_lv
#' @export
StatLv <- ggplot2::ggproto("StatLv", ggplot2::Stat,
  required_aes = c("x", "y"),
  non_missing_aes = "weight",

  setup_params = function(data, params) {
    params$width <- params$width %||% resolution(data$x) * 0.75
    params
  },

  compute_group = function(data, scales, width = NULL, na.rm = FALSE, k = NULL, conf=0.95, percent=NULL) {
    n <- nrow(data)

    k <- determineDepth(n, k, alpha=conf, percent)
    # compute letter values and outliers
    stats <- calcLV(data$y, k)
    outliers <- data$y < min(stats) | data$y > max(stats)
    res <- outputLVplot(data$y, stats, k, outliers, alpha=conf)

    df <- data.frame(res$letter.val)
    df$k <- k
    df$LV <- factor(row.names(df), levels=c("M", LETTERS[6:1], LETTERS[c(26:14, 12:7)]))
    df$ci <- list(res$conf.int)
    df$outliers <- list(data$y[res$outliers])

    df$x <- if (is.factor(data$x)) data$x[1] else mean(range(data$x))
    df$width <- width
    df$relvarwidth <- sqrt(n)
    yrange <- range(data$y)
    if (any(outliers)) {
      yrange <- range(c(res$letter.val[k,-1], data$y[!outliers]), na.rm = TRUE)
    }
    df$ymin <- yrange[1]
    df$ymax <- yrange[2]

    df
  }
)
