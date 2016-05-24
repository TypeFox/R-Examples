#' Emphasize a certain age in Lexis grid
#' 
#' Add a coloured rectangle to an existing Lexis grid to highlight a certain age in that Lexis grid.
#' 
#' @param lg, an existing object originally created with \code{lexis.grid()}.
#' @param age numeric, set the age to highlight.
#' @param fill character, set colour to fill the rectangle. Default is \code{"yellow"}.
#' @param alpha numeric, set alpha, the level of transparency for \code{fill}. Default is \code{0.5}.
#' @details Takes an existing Lexis grid and adds a coloured rectangle that highlights all triangles belonging to a certain age.
#' @return A ggplot2 object.
#' @author Philipp Ottolinger
#' @import ggplot2
#' @importFrom utils tail
#' @export lexis.age
#' @examples 
#' library(LexisPlotR)
#' lexis <- lexis.grid(year.start = 1900, year.end = 1905, age.start = 0, age.end = 5)
#' lexis <- lexis.age(lg = lexis, age = 3)

lexis.age <- function(lg, age, fill = "yellow", alpha = 0.5) {
  age <- as.numeric(age)
  year.start <- as.Date(ggplot_build(lg)$data[[1]][1,1], origin="1970-01-01")
  year.end <- as.Date(tail(ggplot_build(lg)$data[[1]]$xend,1), origin = "1970-01-01")
  age.start <- ggplot_build(lg)$data[[1]][1,3]
  age.end <- tail(ggplot_build(lg)$data[[1]]$yend,1)
  if (age > age.end) { stop("Out of bounds.") }
  if (age < age.start) { stop("Out of bounds.") }
  x <- NULL
  y <- NULL
  df <- data.frame(x = c(year.start, year.end, year.end, year.start),
                   y = c(age, age, age+1, age+1))
  lg <- lg + geom_polygon(data = df, aes(x,y), fill = fill, alpha = alpha, colour = NA)
  return(lg)
}