## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2014-10-09 10:16 emilio on emilio-despacho>
## ============================================================



##' It computes a principal component analysis with sumplementary quantitative and qualitative variables. It is wrapper of \code{\link[FactoMineR]{PCA}}.
##'
##' This function calls \code{\link[FactoMineR]{PCA}} with the the frequency weights as \code{row.w}. Any variable present in \code{freq} are removed from the data. 
##' @title Principal Component Analysis 
##' @param data a data frame
##' @param tfq a table of frequencies
##' @param freq a name of the variable specifying frequency weights 
##' @param scale.unit a boolean, if TRUE (value set by default) then data are scaled to unit variance
##' @param ncp number of dimensions kept in the results
##' @param quantisup a vector indicating the names of the quantitative supplementary variables
##' @param qualisup a vector indicating the names of the categorical supplementary variables
##' @param colw an optional column weights (by default, uniform column weights)
##' @param graph boolean, if TRUE a graph is displayed
##' @param axes a length 2 vector specifying the components to plot
##' @return
##' It returns a list described in \code{\link[FactoMineR]{PCA}}.
##' @seealso \code{\link[FactoMineR]{PCA}}, \code{link{tablefreq}}
##' @importFrom FactoMineR PCA
##' @export
##' @rdname pcafreq
##' @examples
##' pcafreq(iris,  qualisup="Species", graph=TRUE)
##' 
##' tfq <- tablefreq(iris)
##' .pcafreq(tfq,  qualisup="Species", graph=TRUE)
pcafreq <- function(data, freq=NULL, scale.unit = TRUE, ncp = 5, quantisup = NULL, 
                     qualisup = NULL,  colw=NULL, graph = TRUE, axes = c(1, 2)) {
  tablefreq(data, freq=freq) %>% .pcafreq(scale.unit = scale.unit,
                    ncp = ncp,
                    quantisup = quantisup, 
                    qualisup = qualisup,
                    colw=colw,
                    graph = graph,
                    axes = axes)
 } 

##' @rdname pcafreq
##' @export
.pcafreq <- function(tfq, scale.unit = TRUE, ncp = 5, quantisup = NULL, 
                     qualisup = NULL,  colw=NULL, graph = TRUE, axes = c(1, 2)) {
  quanti.sup <- NULL
  quali.sup <- NULL
  col.w <- NULL
  if(!is.null(quantisup))   quanti.sup <- which( colnames(tfq) %in% quantisup)
  if(!is.null(qualisup))   quali.sup <- which( colnames(tfq) %in% qualisup)
  if(!is.null(colw))   col.w <- which( colnames(tfq) %in% col.w)
  PCA(tfq[, -ncol(tfq),drop=FALSE], scale.unit=scale.unit, ncp = ncp, ind.sup = NULL,
      quanti.sup = quanti.sup, 
      quali.sup = quali.sup,
      row.w = unlist(tfq[,ncol(tfq)]),
      col.w = col.w, graph = graph, axes = axes)
} 


