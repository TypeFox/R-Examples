#' @export
#' @name getnns
#' @title Nearest neighbor of individuals based on a dissimilarity matrix
#' @description Get nearest neighbor of individuals and the dissimilarity between individuals and his nearest neighbor
#' @param diss a dissimilarity matrix
#' @param flag a vector of size \code{n} which indicates if we want to compute nearest neighbor of individual i (flag[i]=1) or not (flag[i]=0)
#' @return \item{nn}{a vector of size n which indicates who is the nearest neighbor}
#' @return \item{nndiss}{value of the dissimilarity between individuals an their nearest neighbors}
#'@keywords internal

getnns<-function(diss,flag){
  #function applied to a row vector. Will be used on the whole matrix with "apply"
  mini.vec.flag<-function(vect_flag){
    L<-length(vect_flag)
    flagi<-vect_flag[L]
    vect<-vect_flag[-L]
    if (flagi==0){res<-c(0,0)}
    else{
      search.min<-which(vect==min(vect[vect>0]))[1]
      minobs<-search.min
      mindis<-vect[search.min]
      res<-c(round(minobs,1),mindis)
    }
    return(res)                          
  }  
  
  res<-apply(cbind(diss,flag),1,FUN=mini.vec.flag)
  result<-list(nn=res[1,],nndiss=res[2,])
  return(result)
}