#' Order Points Clockwise or Counter-Clockwise
#' 
#' Returns the indexes of supplied points, \code{x} and \code{y}, ordered either clockwise or 
#' anticlockwise about another point, which by default is taken to be the non-weighted midpoint 
#' of the supplied data
#' @inheritParams getDelaunayMesh
#' @param xm the x value of the reference point
#' @param ym the y value of the reference point
#' @param ... not used
#' @param clockwise order in clockwise or anticlockwise manner
#' @examples
#' #Generate a random set of points and put them clockwise order
#' set.seed(1)
#' x  = runif(100)
#' y  = runif(100)
#' op    = orderPoints(x,y)
#' 
#' #To demonstrate, Lets view the points in order
#' library(ggplot2)
#' df    = data.frame(x,y)
#' df    = df[op,]; 
#' df$id = 1:nrow(df)
#' ggplot(data=df,aes(x,y,colour=id)) + 
#'     geom_path() + geom_point() +
#'     scale_colour_gradient(low="green",high="red")
#' @export
orderPoints = function(x,y,...,xm=mean(range(x)),ym=mean(range(y)),clockwise=TRUE){
  if(!all(c(is.numeric(x),is.numeric(y),is.numeric(xm),is.numeric(ym)))) stop('x,y,xm and ym must be numeric')
  if(length(x) != length(y)) stop('x and y vectors must be the same length')
  order( (if(clockwise){-1}else{1})*atan2(y-ym[1],x-xm[1] ) )
}