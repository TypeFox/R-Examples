#' @title To get the latent variables associated with each cluster
#' 
#' @param resclv : result of CLV(), CLV_kmeans() or LCLV()
#' @param K : the number of groups chosen (already defined if CLV_kmeans is used)
#' 
#' @return \item{comp}{the group latent variables (centered, but not standardized) \cr
#'      For results of LCLV, two types of latent variables are available : \cr
#'      compt : The latent variables of the clusters defined according to the Xr variables, \cr
#'      compc : The latent variables of the clusters defined according to the Xu variables
#'      }
#' 
#'         
#' @examples data(apples_sh)
#' resclvX <- CLV(X = apples_sh$senso, method = "directional", sX = TRUE)
#' comp4G<-get_comp(resclvX, K = 4) 
#' 
#' @export
#' 
get_comp <-
  function(resclv, K=NULL)
  {
    
   if (!inherits(resclv, c("clv","lclv"))) 
      stop("non convenient objects")
    
   X<-resclv$param$X
      
   if(inherits(resclv,"clv")) {
    if(is.null(resclv$param$K)) { 
      if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
      comp<-resclv[[K]]$comp
    } else {
      comp<-resclv$comp
      K<-resclv$param$K
    }
    return(comp=comp)
   }
   
   if(inherits(resclv,"lclv")) {
     if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
     compt<-resclv[[K]]$compt
     compc<-resclv[[K]]$compc
     return(list(compt=compt,compc=compc))
   }
      
  }