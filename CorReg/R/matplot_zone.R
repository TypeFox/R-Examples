#' Matplot with curves comparison by background colors.
#'@description Plot the columns of one matrix against the columns of another, with conditionnal background for easier comparison of curves.
#' @param x the abscisses
#' @param y matrix of the curves (columns)
#' @param col list of colors (like in matplot)
#' @param what a function to choose a winner. Takes y as an input and must return a vector of colors (can be positive integers) of size length(x) 
#' @param alpha parameter for transparency of the background
#' @param ylim ranges of y axe
#' @param xlim ranges of x axe
#' @param type character string (length 1 vector) or vector of 1-character strings indicating the type of plot for each column of y. The first character of type defines the first plot, the second character the second, etc. Characters in type are cycled through; e.g., "pl" alternately plots points and lines.
#' @param xlab title for x axe
#' @param ylab title for y axe
#' @param ... Other graphical parameters 
# ' @param main the main title (like in matplot)
#' @export
#' @examples
#'     \dontrun{

#' require(CorReg)
#' n=15
#' x=1:n
#' y=cbind(c(rnorm(5,0,1),rnorm(5,1,1),rnorm(5,2,1)),
#'          c(rnorm(5,0,1),rnorm(5,1,1),rnorm(5,4,1)),
#'          c(rnorm(5,1,3),rnorm(5,1,2),rnorm(5,1,1)))
#' matplot_zone(x,y,type="l",what=which.max,main="Highest curve")
#' #background color follows color of the highest curve
#' matplot_zone(x,y,type="l",what=which.min,main="Lowest curve")
#' #background color follows color of the lowest curve
#' 
#' }
matplot_zone<-function(x=x,y=y,col=1:6,alpha=0.2,what=which.min,ylim=NULL,xlim=NULL,type="p",xlab=NULL,ylab=NULL,...){
   matplot(x,y,ylim=ylim,type=type,xlab=xlab,ylab=ylab,xlim=xlim,...)
   if(missing(x)){x<-seq_len(NROW(y))}
   victory_int(x=x,y=y,col=col,what=what,alpha=alpha)
   matplot(x,y,ylim=ylim,add=TRUE,col=col,type=type,xlab=xlab,ylab=ylab,xlim=xlim,...)
}