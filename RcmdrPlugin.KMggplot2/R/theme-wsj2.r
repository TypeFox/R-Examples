#' A Wrapper Theme for ggthemes's theme_wsj
#'
#' @param base_size base font size
#' @param base_family font family name
#' @seealso \code{\link[ggplot2:ggplot]{ggplot}} \code{\link[ggthemes:theme_wsj]{theme_wsj}}
#'
#' @rdname theme-wsj2
#' @keywords color
#' @importFrom ggplot2 %+replace%
#' @importFrom ggplot2 theme
#' @importFrom ggthemes theme_wsj
#' @export
theme_wsj2 <- function(base_size = 16, base_family = "") {

  theme_wsj(base_size = base_size, base_family = base_family, title_family = base_family) %+replace%
    theme(
      legend.direction = NULL,
      legend.justification = "center", 
      legend.box = NULL
    )

}
