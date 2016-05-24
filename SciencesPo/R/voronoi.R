#' @title Voronoi diagram
#'
#' @description Computes the voronoi diagram.
#'
#' @param p An integer for the size of the
#' @param n An integer for the size of
#' @param dim The dimension of the image.
#' @param plot Logical, if \code{TRUE}, the plot is returned, else, the \code{data.frame} is returned.
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' @details
#' \url{https://en.wikipedia.org/wiki/Voronoi_diagram}
#' @importFrom data.table := setDT dcast.data.table
#' @export
#' @examples
#' \dontrun{ voronoi(p=2, n=20, dim=1000) }
#'
`voronoi` <- function(p, n=100, dim=1000, plot=TRUE){
  s1 <- Sys.time()
  dim.image <- dim
  colors <- grDevices::rainbow(n)
  points <- data.frame(id = 1:n, x0 = runif(n) * dim, y0 = runif(n) * dim)
  tmp <- data.table::data.table(expand.grid(x = 1:dim.image,
                                y = 1:dim.image, id = 1:n), key = "id")
  tmp <- merge(tmp, points, by = "id")

  `.distance` <- function(a, b, c, d, p){
    (abs(a-c)^p + abs(b-d)^p)^(1/p)}
 # R check barked on global variable on distance.
distance <-NULL
  tmp$distance <- .distance(tmp$x, tmp$y, tmp$x0, tmp$y0, p)

  tmp[, rank := rank(distance, ties.method = "random"), by = c("x", "y")]

  frame <- tmp[tmp$rank == 1,]

if(plot==TRUE){
   frame[,4:7] <-NULL
     frame$color <- colors[frame$id]
     imagen <- as.matrix(data.table::dcast.data.table(data.table::setDT(frame), x ~ y, value.var = "color")[,-1, with=FALSE])
     graphics::frame()
 grid::grid.raster(imagen)

  s2 <- Sys.time()
  timediff <- c( s2 - s1 )
  cat("\n")
  cat("Date of Analysis: ",format(Sys.time(), "%a %b %d %Y"), "\n", "Computation time: ",timediff,sep="","\n")
  cat("-----------------------------------\n")

  }else{
    s2 <- Sys.time()
    timediff <- c( s2 - s1 )
    cat("\n")
    cat("Date of Analysis: ",format(Sys.time(), "%a %b %d %Y"), "\n", "Computation time: ",timediff,sep="","\n")
    cat("-----------------------------------\n")
return(frame)
}
}### end -- voronoi function
NULL
