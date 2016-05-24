#' To get the loadings of the external variables regarding the latent variable in each cluster 
#' 
#' Applies only when external variables (Xr, Xu or both) are involved.
#' 
#' @param resclv : result of CLV(), CLV_kmeans() or LCLV()
#' @param K : the number of groups chosen (already defined if CLV_kmeans is used
#' 
#' @return \item{loading}{the loadings of the external variables \cr
#' For results of LCLV, two types of ladings are defined : \cr
#'       loading_v : loadings of the external Xr variables, \cr
#'       loading_u : loadings of the external Xu variables.
#'      }
#'       
#' @export
#' 
get_load <-
  function(resclv, K=NULL)
  {
    
   if (!inherits(resclv, c("clv","lclv"))) 
      stop("non convenient objects")
    
   X<-resclv$param$X
      
   if(inherits(resclv,"clv")) {
    if(is.null(resclv$param$K)) { 
      if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
      if (is.null(resclv[[K]]$loading)) stop("applies only if external variables are taken into account")
      loading<-resclv[[K]]$loading
    } else {
      if (is.null(resclv$loading)) stop("applies only if external variables are taken into account")
      loading<-resclv$loading
      K<-resclv$param$K
    }
    return(loading=loading)
   }
   
   if(inherits(resclv,"lclv")) {
     if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
     loading_u<-resclv[[K]]$loading_u
     loading_v<-resclv[[K]]$loading_v
     return(list(loading_u=loading_u,loading_v=loading_v))
   }
      
  }