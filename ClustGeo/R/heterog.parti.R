#' @export
#' @name heterog.parti
#' @title heterogeneity of a partition obtained with the method \code{hclust.geo.H}.
#' @description This function calculate the heterogeneity of a partition as the sum of the heterogeneity of all classes
#' obtained with\code{heterog.Ck}.
#' @param dist.var a square symetric matrix with the distances between individuals measured on the variables.
#' @param dist.geo a square symetric matrix with the geographical distances between individuals.
#' @param partition a vector of size \code{n} giving the belonging of each individuals to each class of individuals.
#' @param wt  a vector of size \code{n} whose values indicate weights of individuals on quantitative variables
#' @param alpha a real between 0 and 1 to fit the global criterion to optimize
#' @return {res} {a list containing the following results : the Inertia of the partition based on \code{dist.var},
#'  the Inertia of the partition based on \code{dist.geo}, the heterogeneity of the partition as defined in description.}
#' @details Note that in the calcul of Inertia, the two matrices of distances are normalized (divided by their own maximum).
#' @keywords  internal


heterog.parti<-function(dist.var, dist.geo, partition, wt, alpha){
  nb.k<-length(unique(partition))
  heterog<-0
  In.var<-0
  In.geo<-0
  for (i in 1:nb.k){
    res.heterog.Ck<-heterog.Ck(dist.var, dist.geo, partition, classe=i, wt, alpha)
    heterog<-heterog+res.heterog.Ck$"heterog.Ck"
    In.var<-In.var+res.heterog.Ck$"In.var.Ck"
    In.geo<-In.geo+res.heterog.Ck$"In.geo.Ck"
    
  }
  res<-list(In.var.parti=In.var, In.geo.parti=In.geo, heterog.parti=heterog)
  return(res)
}
