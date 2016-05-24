#' To get the clusters of variables. 
#' 
#' This function returns the group's membership for the p variables. 
#' The output can be a vector p x 1 of integers between 1 and K, 
#' or a binary matrix of size p x n.
#' 
#' 
#' @param resclv : result of CLV(), CLV_kmeans() or LCLV()
#' @param K : the number of groups chosen (already defined if CLV_kmeans is used)
#' @param type : presented in the form of a "vector" (by default) or a "matrix"

#' @return \item{partition}{the group's membership for the variables) }
#'   
#' @examples data(apples_sh)
#' resclvX <- CLV(X = apples_sh$senso, method = "directional", sX = TRUE)
#' parti4G<-get_partition(resclvX, K = 4) 
#' 
#' @export
#' 
get_partition <-
  function(resclv, K=NULL,type="vector")
  {
    
    if (!inherits(resclv, c("clv","lclv"))) 
      stop("non convenient objects")
    
    X<-resclv$param$X
    libel<-colnames(X)
    
    if(is.null(resclv$param$K)) { 
      if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
      clusters<-as.vector(resclv[[K]]$clusters[2,])
    } else {
      clusters<-as.vector(resclv$clusters[2,])
      K<-resclv$param$K
    }
    names(clusters)<-libel
    
    if (type=="vector") return(partition=clusters)
    if (type=="matrix") { 
      clusters<-as.data.frame(clusters)
      names(clusters)<-"G"
      tab<-tabdisj(clusters)
      return(partition=tab)
    }
      
  }