#' @title Draw a coefficient traceplot
#' 
#' @description This function draws a 'traceplot' of coefficient values vs. number of iterations.
#' 
#' @export
#' 
#' @param coef.mat The matrix of coefficients (one coefficient vector per row).
#' @param varnames (Optional) list of variable name labels to make the traceplot more readable.
coef_traceplot <- function(coef.mat,varnames=NULL) {
  numit <- nrow(coef.mat)
  linecolors <- rainbow(ncol(coef.mat))
  plot(NA,NA,main="Coefficient traceplot",type="n",xlim=c(0,numit*1.1),ylim=range(coef.mat),xlab="Iterations",ylab=expression(beta))
  abline(h=0)
  sapply(1:ncol(coef.mat),function(j) { lines(1:numit,coef.mat[,j],col=linecolors[j]) })
  if(!is.null(varnames)) {
    nz.b <- which(coef.mat[nrow(coef.mat),]!=0)
    nz.nms <- varnames[nz.b]
    text(rep(numit,length(nz.nms)),coef.mat[nrow(coef.mat),nz.b],nz.nms,cex=0.7,adj=0,col=linecolors[nz.b])
  }
}