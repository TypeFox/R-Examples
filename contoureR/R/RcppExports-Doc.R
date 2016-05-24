#' Contour Walker Function, Rcpp Interface to C++ Routine
#' 
#' This function is the R interface to the C++ core contour-walker function, totally essential for this package.
#' 
#' The underlying C++ routine establishes a pointer-driven adjacency-network within each of the triangles (Dels) in the 
#' supplied Delaunay mesh, that is to say that Two Del's are deemed to be 'networked' with each other if they share 
#' adjacent edges (Edges), which are drawn between two points (Nodes). Dels are deemed as 'fully-networked' if they hold 
#' reciprocating pointers with exactly three (3) other Dels. Dels can be partially networked if one or two of the 
#' available pointers remain unassigned, which may be the case if the particular Del is a participating member of the 
#' convex hull.
#'  
#' Because C++ uses zero-based indexing, whilst R uses 1-based indexing, should the Delaunay Mesh (\code{dm}) be provided 
#' where the minimum value is 1, then it will be deemed to be a 1-based set and reduced accordingly. 
#' 
#' Prior to traversing the mesh, in order to prevent degeneracies, should any Dels in the mesh have vertices which are 
#' of the same \code{z} value, and/or equal to one of the intended contouring levels, then these nodes will be 
#' pertubated (along the \code{z} direction) by an infitesimally small amount. The consequences of this approach is 
#' that when interpolating along the edges of the Del, the path will always leave at some point (if even trivially small) 
#' inbetween two (2) nodes -- the path will NEVER leave directly through a node, which would otherwise lead to potential 
#' confusion as to the appropriate recipient of the path under such circumstance.
#' 
#' This function is not particularly convenient, and a more convenient wrapper has been produced, with all the usual
#' checks and balances, for further information, see the \code{\link{getContourLines}} function.
#' 
#' @param dm the n x 3 matrix representing the indexes of the vertices of the triangles in the delaunay mesh.
#' No values should be greater than the number of rows in \code{xyz}.
#' @param xyz the m x 3 matrix of xyz coordinates for all the points in the data.
#' @param levels a numeric vector of levels to contour.
#' @param maximumPertubation the maximum pertubation amount (positive or negative) as a percentage.
#' @return \code{matrix} with 6 columns: \code{LevelID, GroupID, PathID, x, y} and \code{z}
#' @rdname contourWalker
#' @name contourWalker
#' @seealso \code{\link{getContourLines}}
#' @examples
#' n = 100
#' x = runif(n)
#' y = runif(n)
#' df = expand.grid(x,y); 
#' colnames(df) = c("x","y")
#' df$z = df$x^2 + df$y^2
#' dm  = getDelaunayMesh(df$x,df$y)
#' res = contourWalker(dm,as.matrix(df),levels=pretty(df$z,n=20))
#' res = data.frame(res); colnames(res) = c('LID','GID','PID','x','y','z')
#' res$Group = interaction(res$LID,res$GID)
#' library(ggplot2)
#' ggplot(res,aes(x,y,group=Group,colour=z)) + geom_path()
NULL

#' Convex Hull via Andrews Monotone, Rcpp Interface to C++ Routine
#' 
#' This function is the R interface to the C++ implementation of Andrews Monotone, a well known algorithm for solving
#' the convex hull in \code{O(nlogn)} time complexity.
#' @param x NumericVector of x values
#' @param y NumericVector of y values
#' @param includeColinear whether to inlude points that line \strong{ON} the hull, by default this is set to FALSE, as this is
#' the true definition of the convex hull.
#' @param  zeroBased Whether the return indexes should be zero based (true, for use in C++), 
#' or One-Based (false, for use in R).
#' @usage 
#' convexHullAM_Indexes(x, y, includeColinear=FALSE,zeroBased = TRUE)
#' convexHullAM_Points(x, y,includeColinear=FALSE)
#' @return \code{convexHullAM_Indexes} returns an integer vector of the indexes of the points, 
#' whilst \code{convexHullAM_Points} returns an \code{n x 2} matrix of the points themselves.
#' @examples
#' library(contoureR)
#' library(ggplot2)
#' set.seed(1)
#' x  = runif(100)
#' y  = runif(100)
#' ch = convexHullAM_Indexes(x,y,includeColinear=FALSE,zeroBased = FALSE)
#' ggplot(data.frame(x,y),aes(x,y)) + 
#'  geom_point() + 
#'  geom_path(data=data.frame(x,y)[ch,],colour="red")
#' @name convexHullAM
#' @rdname convexHullAM
#' @aliases convexHullAM_Indexes convexHullAM_Points
NULL
