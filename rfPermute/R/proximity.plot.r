#' @title Plot Random Forest Proximity Scores.
#' @description Create a plot of Random Forest proximity scores using 
#'   multi-dimensional scaling.
#' 
#' @param rf A \code{randomForest} object.
#' @param dim.x,dim.y Numeric values giving x and y dimensions to plot from 
#'   multidimensional scaling of proximity scores.
#' @param legend.loc Character keyword specifying location of legend. 
#'   See \link{legend}.
#' @param grp.cols Character vector specifying colors for classes.
#' @param circle.size Size of circles around correctly classified 
#'   points as argument to 'cex'. Set to NULL for no circles.
#' 
#' @details Produces a scatter plot of proximity scores for \code{dim.x} and 
#'   \code{dim.y} dimensions from a multidimensional scale (MDS) conversion of 
#'   proximity scores from a \code{randomForest} object. For classification 
#'   models, a convex hull is drawn around the a-priori classes with points 
#'   colored according to original (inner) and predicted (outer) class.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(iris)
#' iris.rf <- randomForest(Species ~ ., data = iris, importance = TRUE, proximity = TRUE)
#' iris.rf
#' proximity.plot(iris.rf, legend.loc = "topleft")
#' 
#' @importFrom stats cmdscale
#' @importFrom grDevices chull rainbow
#' @importFrom graphics plot lines legend points
#' @export
#' 
proximity.plot <- function(rf, dim.x = 1, dim.y = 2, legend.loc = NULL, 
                           grp.cols = NULL, circle.size = 4) {
  if(is.null(rf$proximity)) {
    stop("'rf' has no 'proximity' element. rerun with 'proximity = TRUE'")
  }
  
  prox.cmd <- cmdscale(1 - rf$proximity, k = max(c(dim.x, dim.y)))
  prox.cmd <- prox.cmd[, c(dim.x, dim.y)]
  if(is.null(grp.cols)) grp.cols <- rainbow(length(levels(rf$y)))
  bg.col <- grp.cols[as.numeric(rf$y)]
  plot(prox.cmd[, 1], prox.cmd[, 2], bg = bg.col, pch = 21, 
       xlab = paste("Dimension", dim.x), ylab = paste("Dimension", dim.y)
  )
  
  if(rf$type != "regression") {
    loc.hull <- tapply(1:nrow(prox.cmd), rf$y, function(i) {
      ch <- chull(prox.cmd[i, 1], prox.cmd[i, 2])
      c(i[ch], i[ch[1]])
    })
    for(ch in 1:length(loc.hull)) {
      lines(prox.cmd[loc.hull[[ch]], 1:2], col = grp.cols[ch])
    }
    
    if(!is.null(legend.loc)) {
      legend(legend.loc, legend = levels(rf$y), pch = 21, pt.bg = grp.cols)
    }
    
    if(!is.null(circle.size) & is.numeric(circle.size)) {
      pt.col <- grp.cols[as.numeric(rf$predicted)]
      points(prox.cmd[, 1], prox.cmd[, 2], col = pt.col, pch = 1, cex = circle.size)
    }
  }
  
  invisible(prox.cmd)
}