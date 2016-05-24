lik.contour <- function(x, y, z, levels=NULL, nlevels=11, heat = TRUE, col.heat=NULL, ...)
{
  if(is.null(levels))
  {
    zz <- z[z < Inf & z > -Inf & !is.na(z)]
    levels <- quantile( zz , c(0.001,0.007, seq( 0.05, 0.95, length=nlevels-4),0.993, 0.999) );
    d <- max(-log(min(abs(diff(levels))))/log(10),0);
    levels <- round(levels,d+1);
  }
  if(is.null(col.heat)) col.heat <- heat.colors(length(levels)-1)
  if(heat) image(x,y,z, breaks= sort(levels), col=col.heat,...)
  contour(x,y,z, levels=levels,add=heat,...)
  invisible( list(x=x, y=y, z=z) )
}
