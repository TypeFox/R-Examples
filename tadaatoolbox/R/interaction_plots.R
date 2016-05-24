#' Interaction plots
#'
#' Easily generate interaction plots of two nominal grouping
#' variables and a numeric response variable.
#' @param data A \code{data.frame}
#' @param response Response variable
#' @param group1 First grouping variable
#' @param group2 Second grouping variable
#' @param brewer_palette The name of the \link[RColorBrewer]{RColorBrewer} palette to use, defaults to \code{Set1}
#'
#' @return Invisible: A list with two ggplot2 objects named \code{p1} and \code{p2}.
#' Printed: The two ggplot2 objects.
#' @export
#' @family Tadaa-functions
#' @import ggplot2
#' @examples
#' tadaa_int(ngo, stunzahl, jahrgang, geschl)
tadaa_int <- function(data, response, group1, group2, brewer_palette = "Set1"){

  sdots <- lazyeval::interp(~mean(variable, na.rm = T), variable = substitute(response))

  data <- dplyr::group_by_(data, substitute(group1), substitute(group2))
  data <- dplyr::summarize_(data, .dots = list(mw = sdots))


  p1 <- ggplot(data = data, aes_string(x = substitute(group1), y = "mw", colour = substitute(group2))) +
          geom_point(shape = 23) +
          geom_line(aes_string(group = substitute(group2))) +
          scale_colour_brewer(palette = brewer_palette) +
          labs(title = paste("Interaction of ", substitute(group1), "and",
                             substitute(group2)), y = "Mean")

  p2 <- ggplot(data = data, aes_string(x = substitute(group2), y = "mw", colour = substitute(group1))) +
          geom_point(shape = 23) +
          geom_line(aes_string(group = substitute(group1))) +
          scale_colour_brewer(palette = brewer_palette) +
          labs(title = paste("Interaction of ", substitute(group2), "and",
                             substitute(group1)), y = "Mean")

  print(p1)
  print(p2)

  invisible(list(p1 = p1, p2 = p2))
}
