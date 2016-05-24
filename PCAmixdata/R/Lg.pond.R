#' @export
#' @title Coefficient Lg with ponderation
#' @name Lg.pond
#' @description Compute a data fame with all the Lg coefficients between each matrix of a list ponderated
#' @param liste.mat a list of \code{G} matrix
#' @param ponde a vector of size \code{G} with the ponderation associated to each matrix
#' @return \item{Lg.pond}{a data frame with \code{G} rows and \code{G} colums with all the Lg coefficients ponderated}
#' 
#' @examples
#' V0<-c("a","b","a","a","b")
#' V01<-c("c","d","e","c","e")
#' V1<-c(5,4,2,3,6)
#' V2<-c(8,15,4,6,5)
#' V3<-c(4,12,5,8,7)
#' V4<-c("vert","vert","jaune","rouge","jaune")
#' V5<-c("grand","moyen","moyen","petit","grand")
#' G1<-data.frame(V0,V01,V1)
#' G2<-data.frame(V2,V3)
#' G3<-data.frame(V4,V5)
#' liste.mat<-list(G1,G2,G3)
#' Lg(liste.mat)

Lg.pond<-function(liste.mat,ponde){
  Lg.pond<-Lg(liste.mat)
  Lg.pond<-sweep(Lg.pond,1,STATS=ponde,FUN="*")
  Lg.pond<-sweep(Lg.pond,2,STATS=ponde,FUN="*")
}