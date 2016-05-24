plotVector <-
function( x1, y1, x2=NULL, y2=NULL, add.vectors=TRUE, col=c("black","red","darkgreen","purple"), 
  lwd=8, font=2, font.lab=2, las=1, cex.lab=1.3, cex.axis=2, usr=NULL, ... )   {
  # Plots one or two vectors and the vector sum.
  # 'x1':  \code{x}-coordinate of the first vector.
  # 'y1':  \code{y}-coordinate of the first vector.
  # 'x2':  \code{x}-coordinate of the second vector.
  # 'y2':  \code{y}-coordinate of the second vector.
  # 'add.vectors':  Logical; prints the vector sum if \code{TRUE} (default); otherwise, does not print the vector sum.
  # 'col':  Scalar or vector, specifying the colors of the plot.  Type \code{colors()} for selections.
  # 'lwd': The line width (see \code{\link[graphics]{par}}).
  # 'font': An integer which specifies which font to use for text (see \code{\link[graphics]{par}}). 
  # 'font.lab': The font to be used for \code{x} and \code{y} labels (see \code{\link[graphics]{par}}). 
  # 'las': Numeric in (0,1,2,3); the style of the axis labels (see \code{\link[graphics]{par}}). 
  # 'cex.lab': The magnification to be used for \code{x} and \code{y} labels (see \code{\link[graphics]{par}}). 
  # 'cex.axis': The magnification to be used for axis annotation (see \code{\link[graphics]{par}}). 
  # 'usr': A vector of the form \code{c(x1, x2, y1, y2)} giving the extremes of the usr coordinates of the 
  #        plotting region (see \code{\link[graphics]{par}}).
  # ...: Optional arguments to \code{\link[graphics]{plot}}.
  # Examples:
  #        # Vectors (2,8) and (4,-3) and their vector sum.
  #        plotVector( 2, 8, 4, -3 ) 
  #     
  #        # Colinear vectors (-3,6) and (-1,2).
  #        plotVector( -3, 6, -1, 2, add=FALSE, col=c("red","black") )
  #
  #        # Colinear vectors (-1,2) and (3,-6).
  #        plotVector( -1, 2, 3, -6, add=FALSE )
  #
  #        # Vectors (2,3) and (5,-4)
  #        plotVector( 2, 3, 5, -4, add=FALSE, usr=c( -5, 5, -4, 7) )
  #
  #        # Vectors (2,3) and (-5,4) and their vector sum.
  #        plotVector( 2, 3, -5, 4, usr=c( -5, 5, -4, 7 ) )
  if (length(col)==1) {col <- rep(col, 4)};  if (length(col)==2) {col <- c(col, col)};  if (length(col)==3) {col <- c(col, col[1])}
  if ( is.null(x2) | is.null(y2) )  x2 <- y2 <- 0
  f1=function(x){x*y1/x1} ;   f2=function(x){x*y2/x2} ;   f3=function(x){x*(y1+y2)/(x1+x2)}; vertical=function(x){y1+y2}
  f4=function(x){ (x2*y1-x1*y2+x*y2) / x2 } ;     f5=function(x){ (x1*y2-x2*y1+x*y1) / x1 }
  x.vec = c(x1,x2,(x1+x2)*add.vectors,0) ;   y.vec = c(y1,y2,(y1+y2)*add.vectors,0)
  x10 = min(x.vec) - ( max(x.vec) - min(x.vec) ) / 15
  x20 = max(x.vec) + ( max(x.vec) - min(x.vec) ) / 15
  y10 = min(y.vec) - ( max(y.vec) - min(y.vec) ) / 15
  y20 = max(y.vec) + ( max(y.vec) - min(y.vec) ) / 15
  if ( !is.null(usr) )  {    
     x10 = min( x10, usr[1] ) ;  x20 = max( x20, usr[2] ) ;  y10 = min( y10, usr[3] ) ;  y20 = max( y20, usr[4] ) }
  plot( c(x10,x10,x20,x20), c(y10,y20,y10,y20), col="white", xlab="X", ylab="Y", 
     lwd=lwd, font=font, font.lab=font.lab, las=las, cex.lab=cex.lab, cex.axis=cex.axis, ... )
  curve( f1, 0, x1, add=TRUE, cex=5, col=col[1],  
     lwd=lwd, font=font, font.lab=font.lab, las=las, cex.lab=cex.lab, cex.axis=cex.axis, ... )
  if ( x2 != 0 )   {
     curve( f2, 0, x2, add=TRUE, col=col[2],
        lwd=lwd, font=font, font.lab=font.lab, las=las, cex.lab=cex.lab, cex.axis=cex.axis, ... )
     if (add.vectors)  {
        if ( x1+x2 != 0 )  curve( f3, 0, x1+x2, add=TRUE, col=col[3], 
           lwd=lwd, font=font, font.lab=font.lab, las=las, cex.lab=cex.lab, cex.axis=cex.axis, ... )
        if ( x1+x2 == 0 )  curve( vertical, 0, 0, n=1, type="h", add=TRUE, col=col[3],
           lwd=lwd, font=font, font.lab=font.lab, las=las, cex.lab=cex.lab, cex.axis=cex.axis, ... ) 
        curve( f4, x1, x1+x2, add=TRUE, lty=2, col=col[4],
           lwd=lwd, font=font, font.lab=font.lab, las=las, cex.lab=cex.lab, cex.axis=cex.axis, ... )
        curve( f5, x2, x1+x2, add=TRUE, lty=2, col=col[4],
           lwd=lwd, font=font, font.lab=font.lab, las=las, cex.lab=cex.lab, cex.axis=cex.axis, ... )  }
  }
  abline( h=0, v=0, lty=3, col="black",
     lwd=lwd, font=font, font.lab=font.lab, las=las, cex.lab=cex.lab, cex.axis=cex.axis, ... )
}
