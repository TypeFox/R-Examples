#' A Simple Black and White Theme
#'
#' @param base_size base font size
#' @param base_family font family name
#' @seealso \code{\link[ggplot2:ggplot]{ggplot}}
#'
#' @rdname theme-simple
#' @keywords color
#' @importFrom ggplot2 %+replace%
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @export
theme_simple <- function(base_size = 16, base_family = "") {

  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
    strip.background = element_rect(fill = "grey90", colour = "grey90"),
    panel.border     = element_rect(fill = "transparent", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

}
