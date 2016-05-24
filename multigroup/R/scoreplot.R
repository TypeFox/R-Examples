#' @title Score plot for multigroup data
#' @description plots of individuals 
#' @param x results of the proposed multigroup methods in the package
#' @param axes a vector of two selected components 
#' @param cex character expansion for text by default .85
#' @param font.lab type of font by default 3
#' @return score plot
#' @export
#'  
#' @examples
#' Data = iris[,-5]
#' Group = iris[,5]
#' res.mgPCA = mgPCA (Data, Group, graph=TRUE)
#' scoreplot(res.mgPCA, axes=c(1,2))
scoreplot <- function(x, axes=c(1,2), cex=NULL, font.lab= NULL){
  #============================================================================
  #                            Preparing inputs
  #============================================================================
  AA=x$loadings.common
  if (max(axes)>ncol(AA))
    stop("\nOops one of the axes value is greater than number of existing dimensions")
  
  if(is.null(cex)) {cex = .85}
  if(is.null(font.lab)) {font.lab = 3}

  Group=x$Group
  xax=axes[1]
  yax=axes[2]
  lab.x <- paste("Dim ", axes[1],sep = "")
  lab.y <- paste("Dim ", axes[2],sep = "")
  #ggplot colors
   cooll <- function(n, alfa=1) {
    hues = seq(15,375,length=n+1)
    hcl(h=hues, l=65, c=100, alpha=alfa)[1:n]
  }
  
  #============================================================================
  #                            Score plot
  #============================================================================
  TT = x$Con.Data %*% x$loadings.common
  Ts = TT [,c(xax,yax)]
  rownames(Ts) = Group
  pervar = x$Lambda[[1]]
  cp1 <- round(var(Ts[,xax])/sum(diag(var(Ts))), digits = 2)
  cp2 <- round(var(Ts[,yax])/sum(diag(var(Ts))), digits = 2)
  lab.x <- paste("Dim ", xax, sep = "")
  lab.y <- paste("Dim ", yax, sep = "")
  labs_col = cooll(nlevels(Group), alfa=1)
  rep.labs_col=rep(labs_col, as.vector(table(Group)))
  plot(Ts[,xax], Ts[,yax], xlab = lab.x, ylab = lab.y, type="n", main = "Individual plot") 
  abline(h = 0, v = 0, col= "gray60")
  text(Ts, labels=rownames(Ts), cex=cex, font.lab=font.lab, col=rep.labs_col)

    
}
