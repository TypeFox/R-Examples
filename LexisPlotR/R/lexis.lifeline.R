#' Plot lifelines into a Lexis grid
#' 
#' Add lifelines to an existing Lexis grid.
#' 
#' @param lg, an existing object originally created with \code{lexis.grid()}.
#' @param entry character, set the entry or birth date of an individual in format \code{"YYYY-MM-DD"}.
#' @param exit character, set the exit or death date of an individual in format \code{"YYYY-MM-DD"}. Default is \code{NA} (no exit or death observed).
#' @param lineends logical, if \code{TRUE} lineends will be marked. Default is \code{FALSE}.
#' @param colour character, set the colour of the lifelines. Default is \code{"red"}.
#' @param alpha numeric, set the transparency of the lifelines. Default is \code{1} (no transparency).
#' @param lwd numeric, set the linewidth of the lifelines. Default is \code{0.5}.
#' @details Takes an existing Lexis grid and adds lifelines to the grid. Input can be a single dates or dates from a vector.
#' @return A ggplot2 object.
#' @author Philipp Ottolinger
#' @import ggplot2
#' @importFrom utils tail
#' @export lexis.lifeline
#' @examples 
#' lg <- lexis.grid(year.start = 1900, year.end = 1905, age.start = 0, age.end = 5)
#' lexis.lifeline(lg = lg, entry = "1901-09-23")
#' lexis.lifeline(lg = lg, entry = "1901-09-23", exit = "1904-03-03")
#' 
lexis.lifeline <- function(lg, entry, exit = NA, lineends = F, colour = "red", alpha = 1, lwd = 0.5) {
  if (!is.ggplot(lg)) { stop("No valid ggplot object.") }
  entry <- as.Date(entry, origin = "1970-01-01")
  exit <- as.Date(exit, origin = "1970-01-01")
  year.start <- as.Date(ggplot_build(lg)$data[[1]][1,1], origin="1970-01-01")
  year.end <- as.Date(tail(ggplot_build(lg)$data[[1]]$xend,1), origin = "1970-01-01")
  age.start <- ggplot_build(lg)$data[[1]][1,3]
  age.end <- tail(ggplot_build(lg)$data[[1]]$yend,1)
  x <- NULL
  y <- NULL
  xend <- NULL
  yend <- NULL
  case <- data.frame(entry, exit)
  case$x <- entry
  case$xend <- ifelse(is.na(exit), year.end, exit)
  case$xend <- as.Date(case$xend, origin = "1970-01-01")
  case$y <- 0
  case$yend <- ifelse(is.na(case$exit), how.old(case$entry, year.end), how.old(case$entry, case$exit))
  lg <- lg + geom_segment(data = case, aes(x=x,xend=xend,y=y,yend=yend), colour = colour, alpha = alpha, lwd = lwd)
  if (lineends == TRUE) {
    lg <- lg + geom_point(data = case[!is.na(case$exit),], aes(x=xend, y=yend), size=2, shape = 3)
  }
return(lg)
}

