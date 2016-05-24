#' Emphasize a certain year in Lexis grid.
#' 
#' Takes an existing Lexis grid and adds a coloured rectangle to highlight a certain age.
#' 
#' @param lg, an existing object originally created with \code{lexis.grid()}.
#' @param year numeric, set the year to highlight.
#' @param fill character, set the colour of the rectangle. Default is \code{"blue"}.
#' @param alpha numeric, set the transparency of the rectangle. Default is \code{0.5}.
#' @details Takes an existing Lexis grid and adds a coloured rectangle to the plot. The rectangle will highlight a certain year in the grid.
#' @return A ggplot2 object.
#' @author Philipp Ottoliner
#' @import ggplot2
#' @importFrom utils tail
#' @export lexis.year
#' @examples 
#' lg <- lexis.grid(year.start = 1900, year.end = 1905, age.start = 0, age.end = 5)
#' lexis.year(lg = lg, year = 1902)

lexis.year <- function(lg, year, fill = "blue", alpha = 0.5) {
  year.start <- as.Date(ggplot_build(lg)$data[[1]][1,1], origin="1970-01-01")
  year_start <- as.numeric(substr(year.start,1,4))
  year.end <- as.Date(tail(ggplot_build(lg)$data[[1]]$xend,1), origin = "1970-01-01")
  year_end <- as.numeric(substr(year.end,1,4))
  age.start <- ggplot_build(lg)$data[[1]][1,3]
  age.end <- tail(ggplot_build(lg)$data[[1]]$yend,1)
  if (year > year_end - 1) { stop("Out of bounds.") }
  if (year < year_start) { stop("Out of bounds.") }
  year1 <- year + 1
  year <- as.Date(paste(year, "-01-01", sep = ""), origin = "1970-01-01")
  year1 <- as.Date(paste(year1, "-01-01", sep = ""), origin = "1970-01-01")
  x <- NULL
  y <- NULL
  df <- data.frame(x = c(year, year1, year1, year),
                   y = c(age.start, age.start, age.end, age.end))
  lg <- lg + geom_polygon(data = df, aes(x,y), fill = fill, alpha = alpha, colour = NA)
  return(lg)
}