#' A Theme for No of At Risk in Kaplan-Meier Plots
#'
#' These functions generate graphics themes for no of at risk in Kaplan-Meier Plots.
#'
#' @param base_size base font size
#' @param base_family font family name
#' @param base_theme base theme name
#' @seealso \code{\link[ggplot2:ggplot]{ggplot}}
#'
#' @rdname theme-natrisk
#' @keywords color
#' @importFrom ggplot2 %+replace%
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @export
theme_natrisk <- function(base_theme, base_size = 20, base_family = "") {
  
  plotbg   <- base_theme()$plot.background$fill

  natheme <- base_theme(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line        = element_blank(),
      axis.text.x      = element_blank(),
      axis.title.x     = element_blank(),
      legend.position  = "none",
      panel.background = element_rect(fill = "transparent", colour = "transparent"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      plot.background  = element_rect(fill = "transparent", colour = "transparent"),
      plot.title       = element_blank()
    )
  natheme$axis.title.y$colour <- "transparent"
  natheme$axis.text.y$colour <- "transparent"
  natheme$axis.ticks$colour <- "transparent"
  natheme$plot.margin[1] <- 0
  natheme$plot.margin[3] <- 0
  if (!is.null(natheme$axis.ticks.x)) natheme$axis.ticks.x$colour <- "transparent"
  if (!is.null(natheme$axis.ticks.y)) natheme$axis.ticks.y$colour <- "transparent"
  
  natheme
  
}