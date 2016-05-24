#' Fill Lexis triangles by HMD data
#' 
#' The function opens an existing Lexis grid and fill the triangles according to data from the Human Mortality Database.
#' 
#' @param lg, an existing object originally created with \code{lexis.grid()}.
#' @param hmd.data, a data.frame created with \code{prepare.hmd()}.
#' @param column character, the name of the column of \code{hmd.data} the triangles shall be filled with.
#' @details The function creates a subset of \code{hmd.data} that fits in the dimensions of the existing Lexis grid.
#' The triangles will be filled according to the data in \code{column}.
#' @author Philipp Ottolinger
#' @import ggplot2
#' @importFrom utils tail
#' @export lexis.hmd
#' @examples 
#' library(LexisPlotR)
#' lg <- lexis.grid(year.start = 1980, year.end = 1985, age.start = 0, age.end = 5)
#' # Load sample data
#' path <- system.file("extdata", "Deaths_lexis_sample.txt", package = "LexisPlotR")
#' deaths.triangles <- prepare.hmd(path)
#' lexis.hmd(lg = lg, hmd.data = deaths.triangles, column = "Total")
#' 
#' ### Plot data not explicitly present in HMD data
#' deaths.triangles$RatioMale <- deaths.triangles$Male / deaths.triangles$Total
#' lexis.hmd(lg, deaths.triangles, "RatioMale")

lexis.hmd <- function(lg, hmd.data, column) {
  year.start <- as.Date(ggplot_build(lg)$data[[1]][1,1], origin="1970-01-01")
  year_start <- as.numeric(substr(year.start, 1, 4))
  year.end <- as.Date(tail(ggplot_build(lg)$data[[1]]$xend,1), origin = "1970-01-01")
  year_end <- as.numeric(substr(year.end, 1, 4))
  age.start <- ggplot_build(lg)$data[[1]][1,3]
  age.end <- tail(ggplot_build(lg)$data[[1]]$yend,1)
  filterYear <- year_start:(year_end - 1)
  filterAge <- age.start:(age.end - 1)
  data <- hmd.data[hmd.data$Year %in% filterYear & hmd.data$Age %in% filterAge,]
  n <- dim(data)[1]
  for (i in 1:n) {
    xx <- c(data[i,"x1"],data[i,"x2"],data[i,"x3"])
    yy <- c(data[i,"y1"],data[i,"y2"],data[i,"y3"])
    zz <- data[i, column]
    df <- data.frame(xx,yy,zz)
    lg <- lg + geom_polygon(data = df, aes(x = xx, y = yy, fill=zz))
  }
  lg <- lg + labs(fill = column)
  return(lg)
}