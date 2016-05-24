#' @export
#' @name heterog.Ck
#' @title heterogeneity of a class C_k obtained with the method \code{hclust.geo.H}.
#' @description This function calculate the heterogeneity of a class C_k obtained with the method \code{hclust.geo.H}.
#' heterogeneity= alpha I_intra_var + (1-alpha) I_intra_geo.
#' @param dist.var a square symetric matrix with the distances between individuals measured on the variables.
#' @param dist.geo a square symetric matrix with the geographical distances between individuals.
#' @param partition a vector of size \code{n} giving the belonging of each individuals to each class of individuals.
#' @param classe the number of the class for which we want the heterogeneity.
#' @param wt  a vector of size \code{n} whose values indicate weights of individuals on quantitative variables
#' @param alpha a real between 0 and 1 to fit the global criterion to optimize
#' @return {res} {a list containing the following results : the Inertia of the class based on \code{dist.var},
#'  the Inertia of the class based on \code{dist.geo}, the heterogeneity of the class as defined in description.}
#' @details Note that in the calcul of Inertia, the two matrices of distances are normalized (divided by their own maximum).
#' @keywords internal



heterog.Ck<-function(dist.var, dist.geo, partition, classe, wt, alpha){
  
  dist.var<-(dist.var/max(dist.var))^2
  dist.geo<-(dist.geo/max(dist.geo))^2
  
  ind_k<-which(partition==classe)
  wt_k<-wt[ind_k]
  mu_k<-sum(wt_k)
  
  dist.var_k<-dist.var[ind_k, ind_k]
  dist.geo_k<-dist.geo[ind_k, ind_k]
  
  
  In.var<-sweep(dist.var_k, 1, FUN="*", STATS=wt_k)
  In.var<-sweep(dist.var_k, 2, FUN="*", STATS=wt_k)
  In.var<-In.var/(2*mu_k)
  In.var<-sum(In.var)
  
  In.geo<-sweep(dist.geo_k, 1, FUN="*", STATS=wt_k)
  In.geo<-sweep(dist.geo_k, 2, FUN="*", STATS=wt_k)
  In.geo<-In.geo/(2*mu_k)
  In.geo<-sum(In.geo)
  
  heterog<-alpha*In.var + (1-alpha)*In.geo
  res<-list(In.var.Ck=In.var, In.geo.Ck=In.geo, heterog.Ck=heterog, alpha=alpha)
  return(res)
}
