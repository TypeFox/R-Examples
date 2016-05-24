#'@title Graphical representation of the CLV clustering stages
#'
#'@description
#'This function plots either the CLV dendrogram or the variations of the consolidated CLV criterion.
#'
#' @param x : an object of class \code{clv}
#' @param type : What to plot. \cr 
#'      "dendrogram" : the dendrogram of the hierchical clustering algorithm, \cr    
#'      "delta"      : a barplot showing the variation of the clustering criterium after consolidation.
#' @param cex : Character expansion for labels.
#' @param ... further arguments passed to or from other methods
#' 
#' @seealso CLV
#' 
#' @export

plot.clv =  function (x, type="dendrogram",cex=0.8,...) 
{
  if (!inherits(x, "clv"))    stop("non convenient object")
  
  resclv<-x
  if(is.null(resclv$param$nmax)) stop("plot_clv() shall apply only for hierarchical analysis")
   
 if (type=="dendrogram") {
  par(cex=cex)
  plot(resclv$mydendC, type ="rectangle",  main="CLV Dendrogram", axes=F, cex=cex)
  par(cex=1)
 }
  
 if (type=="delta") {
  p<-resclv$param$p
  nmax<-resclv$param$nmax
  sbegin<-resclv$param$sbegin
  results<-resclv$tabres
  if (p>nmax) gpmax<-nmax
  if (p<=nmax) gpmax<-p

  tempo<-(results[(p-2):(p-gpmax+1),7]-results[(p-1):(p-gpmax+2),7])
  if (results[1,7]>0) tempo<-c(tempo,sbegin-results[1,7])
  if (results[1,7]==0) tempo<-c(tempo,results[p-gpmax,7]-results[p-gpmax+1,7])
  tempo[which(tempo<0)]<-0
  barplot(tempo[(gpmax-1):1],col=4,xlab="Nb clusters", ylab="delta", 
          main="Variation of criterion (after consolidation)",
          axisnames=TRUE,names.arg=paste((gpmax):2,"->",(gpmax-1):1),
          las=2,cex.names=cex,cex.main = 0.8)
 }
 
}