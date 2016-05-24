#' Edgelist to sociomatrix 
#' 
#' Construction of a sociomatrix from an edgelist
#' 
#' @usage el2sm(el,directed=TRUE,nadiag=all(el[,1]!=el[,2]))
#' @param el a matrix in which each row contains the indices of an edge and possibly the weight for the edge
#' @param directed if FALSE, then a relation is placed in both entry ij and ji of the sociomatrix, for each edge ij (or ji)
#' @param nadiag put NAs on the diagonal 
#' @return a sociomatrix 
#' @author Peter Hoff
#' @examples
#' 
#' Y<-matrix(rpois(10*10,.5),10,10) ; diag(Y)<-NA
#' E<-sm2el(Y) 
#' el2sm(E) - Y 
#' 
#' @export el2sm
el2sm<-function(el,directed=TRUE,nadiag=all(el[,1]!=el[,2]))
{ 
  w<-rep(1,nrow(el))
  if(ncol(el)>2){ w<-el[,3] }

  if( is.numeric(el) && all(round(el[,1:2])==el[,1:2])  ) { nodes<-1:max(el) }
  if(!(is.numeric(el) && all(round(el[,1:2])==el[,1:2])))
  {
    nodes<-sort(unique(c(el[,1:2])))
  }

  el<-cbind( match(el[,1],nodes) ,  match(el[,2],nodes) )

  n<-max(el[,1:2])
  sm <- matrix(0,n,n)            # construct sociomatrix 
  sm[el[,1:2]]<-w                # fill in 
  if(nadiag) { diag(sm) <- NA  } # set diagonal to NA 
  if(!directed){ sm<-sm+t(sm) }
  dimnames(sm)[[1]]<-dimnames(sm)[[2]]<-nodes
  sm
}

