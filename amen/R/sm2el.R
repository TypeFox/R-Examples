#' Sociomatrix to edgelist
#' 
#' Construction of an edgelist from a sociomatrix
#' 
#' @usage sm2el(sm,directed=TRUE)
#' @param sm a sociomatrix with possibly valued relations
#' @param directed if TRUE, only use the upper triangular part of the matrix to enumerate edges 
#' @return an edglist 
#' @author Peter Hoff
#' @examples
#' 
#' Y<-matrix(rpois(10*10,.5),10,10) ; diag(Y)<-NA
#' E<-sm2el(Y) 
#' el2sm(E) - Y 
#' 
#' @export sm2el
sm2el<-function(sm,directed=TRUE)
{
  if(!directed){ sm[lower.tri(sm)]<-0 }
  el<-which(sm!=0,arr.ind=TRUE)
  w<-sm[el]
  if(var(w)>0) { el<-cbind(el,w) }
  el[order(el[,1]),]
}

