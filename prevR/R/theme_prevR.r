#' prevR themes for ggplot2
#' 
#' Two custom themes for ggplot2 graphs, hidding axis.
#' 
#' @param base_size base font size
#' 
#' @seealso \code{\link[ggplot2]{ggtheme}}\{\pkg{ggplot2}\}
#' 
#' @export

theme_prevR <- function (base_size = 12) {
  '%+replace%'(theme_grey(base_size),
    theme(
      axis.title        = element_blank(),
      axis.text         = element_blank(),
      axis.ticks.length = unit(0, "cm"),
      axis.ticks.margin = unit(0, "lines"),
      plot.margin       = unit(c(0, 0, 0, 0), "lines"),
      complete          = TRUE
    )
  )
}

#' @export
#' @rdname theme_prevR

theme_prevR_light <- function (base_size = 12) {
  '%+replace%'(theme_grey(base_size),
    theme(
      axis.title        = element_blank(),
      axis.text         = element_blank(),
      panel.background  = element_blank(),
      panel.grid        = element_blank(),
      axis.ticks.length = unit(0, "cm"),
      axis.ticks.margin = unit(0, "lines"),
      plot.margin       = unit(c(0, 0, 0, 0), "lines"),
      complete          = TRUE
    )
  )
}
