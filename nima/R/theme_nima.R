#' Nima's ggplot2 theme
#'
#' Nima's ggplot2 theme: black border and no ticks
#'
#' @param base_size Base font size
#' @param base_family Base font family
#' @param ... Passed to \code{\link[ggplot2]{theme}}
#'
#' @importFrom ggplot2 ggplot
#'
#' @export
#'
#' @return An object as returned by \code{\link[ggplot2]{theme}}
#'
#' @examples
#' library(ggplot2)
#' mtcars$cyl <- factor(mtcars$cyl)
#' p <- ggplot(mtcars, aes(y=mpg, x = disp, color = cyl)) 
#' p <- p + geom_point() + theme_nima()
#' p
#'
#' @seealso
#' \code{\link[ggplot2]{theme}}


theme_nima <- function(base_size = 12, base_family = "", ...) {
    ggplot2::"%+replace%"(
      ggplot2::theme_grey(base_size = base_size, base_family = base_family),
      ggplot2::theme(axis.ticks.length = grid::unit(0, "cm"),
                     panel.border = ggplot2::element_rect(fill = NA,
                                                          color = "black"),
                     strip.background=ggplot2::element_rect(fill="gray80",
                                                            color="black"), ...)
      )
}

#' @rdname theme_nima
#' @export
nima_theme <- theme_nima
