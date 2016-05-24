#' Add lines 
#' 
#' Add lines to a network plot
#' 
#' @usage addlines(Y,X,col="lightblue",alength=0,...)
#' @param Y a sociomatrix 
#' @param X coordinates of nodes
#' @param col color of lines. Can be a vector of length equal to the number of edges to be drawn
#' @param alength length of arrows to be drawn
#' @param \ldots additional plotting parameters
#' @author Peter Hoff
#' @examples
#'
#' data(addhealthc3) 
#' Y<-addhealthc3$Y
#' X<-xnet(Y) 
#' netplot(Y,X) 
#' addlines(Y,X,col=Y[Y!=0]) 
#' 
#' @export addlines
addlines<-function(Y,X,col="lightblue",alength=0,...)
{
  # add links between nodes of a sociomatrix 
  links<-which(Y!=0,arr.ind=TRUE)
  suppressWarnings(
  arrows(X[links[,1],1],X[links[,1],2],X[links[,2],1],X[links[,2],2],
           col=col,length=alength,...)
                   )
}



