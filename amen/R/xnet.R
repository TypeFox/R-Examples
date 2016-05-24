#' @title Network embedding
#' 
#' @description 
#' Compute an embedding of a sociomatrix into a two-dimensional space. 
#' 
#' @param Y (square matrix) The sociomatrix. 
#' @param fm (logical scalar) If TRUE, the Fruchterman-Reingold layout will 
#' be used (requires the network package). 
#' @param seed (integer) The random seed (the FR layout is random).  
#' 
#' @details 
#' Coordinates are obtained using the Fruchterman-Reingold layout if the 
#' package \code{network} is installed, and otherwise uses the first two 
#' eigenvectors the sociomatrix. 
#' 
#' @return (matrix) A matrix of two-dimensional coordinates. 
#' 
#' @author Peter Hoff
#' 
#' @examples
#' data(addhealthc3) 
#' Y<-addhealthc3$Y
#' X<-xnet(Y) 
#' netplot(Y,X) 
#' 
#' @export xnet
xnet<-function(Y,fm=suppressWarnings(require("network")),seed=1)
{
  if(!is.null(seed)) { set.seed(seed) }
  if(fm)
  {
    x<-network::as.network(Y)
    n <- network::network.size(x)
    d <- network::as.matrix.network(x, matrix.type = "adjacency")
    d[is.na(d)] <- 0
    d <- network::as.network(matrix(as.numeric(d > 0), n, n))
    U<-network::network.layout.fruchtermanreingold(d,layout.par=NULL)
  }
  if(!fm)
  {
    Y0<-Y ; diag(Y0)<-1
    U<-Re(eigen(Y0)$vec[,1:2])
  }
  U<-sweep(U,2,apply(U,2,mean)) ; U<-sweep(U,2,apply(U,2,sd),"/")
  U
}



