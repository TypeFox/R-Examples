#' Get Delaunay Mesh
#' 
#' Retrieves the Delaunay Mesh for a series of x and y points in 2D. 
#' With the exception of a few brief checks, is almost a direct wrapper to the \code{\link[geometry]{delaunayn}} 
#' function as part of the geometry package.
#' @param x numeric vector of x values
#' @param y numeric vector of y values of same length as x
#' @return \code{matrix} object having three columns that represent the (1-based) indexes of each vertex 
#' relative to the data in the \code{x} and \code{y} input parameters.
#' @examples
#' #Generate a sample Delaunay Mesh
#' set.seed(1)
#' x  = runif(100)
#' y  = runif(100)
#' dm = getDelaunayMesh(x,y)
#' 
#' #To demonstrate, Lets view the mesh
#' library(ggplot2)
#' library(reshape)
#' df = as.data.frame(dm); df$id = 1:nrow(df); df = melt(df,id="id")
#' df = cbind(df,data.frame(x,y)[df$value,])
#' ggplot(data=df,aes(x,y,group=id)) + 
#'  geom_polygon(aes(fill=id),color="gray")
#' @export
getDelaunayMesh <- function(x,y){
  if(!all(is.numeric(x),is.numeric(y)))stop('x and y must both be numeric')
  if(length(x) != length(y)) stop('x and y must be of the same length')
  tryCatch({ 
    dm = suppressMessages(delaunayn(matrix(c(x,y),nrow=length(x),ncol=2),options="Qt"))
  },error=function(e){
    writeLines(as.character(e))
    dm <<-matrix(0,0,3)
  })
  dm
}