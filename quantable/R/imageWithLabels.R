#' image plot with labels
#'
#' @param x - matrix
#' @param labels - colnames(x)
#' @param cex.axis - size of axis lables
#' @param cex - size of labels
#' @param main - main title
#' @param col - color map for matrix
#' @export
#' @examples
#' x = matrix(rnorm(20*20),ncol=20)
#' imageWithLabels(x)
imageWithLabels = function(x, labels=colnames(x),cex=1,cex.axis=0.5,main=NULL,col = heat.colors(12))
{
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(3,1), heights=c(1,1))

  image(x, axes = F, main =main, col=col)
  axis( 1, at=seq(0,1,length=length((labels))) , labels=labels,cex.axis=cex.axis, las=2, cex=cex )
  axis( 2, at=seq(0,1,length=length((labels))) , labels=labels,cex.axis=cex.axis, las=1, cex=cex )

  colorlevels = seq(min(x,na.rm = TRUE),max(x,na.rm = TRUE),length=length(col))
  image(1, seq(0,1,length=length(colorlevels)),
        matrix(data=colorlevels, nrow=1),
        col=col,xlab="",ylab="",
        axes=FALSE)
  axis( 2, at=seq(0,1,length=length((colorlevels))) , labels=round(colorlevels,digits=2),cex.axis=cex.axis, las=1, cex=cex )
  layout(1)
}
#' if you need an colorscale to you imagelables use this
#' @param data the data matrix
#' @param colors used
#' @export
colorscale = function(data,colors=heat.colors(12)){
  nrc = length(colors)
  z  = seq( min(data) , max(data) , length=nrc)
  image(1, seq(0,1,length=nrc), matrix(z,1,nrc) ,axes=F,ylab="",xlab="")
  axis( 2, at=seq(0,1,length=nrc) , labels=round(z,digits=2), las=2 )
}
