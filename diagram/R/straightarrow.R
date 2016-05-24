
##==============================================================================
# straightarrow: Plot straight arrow at certain distance between two points
##==============================================================================

straightarrow <- function(from, to, lwd=2, lty=1, lcol="black",
    arr.col=lcol, arr.pos=0.5, endhead=FALSE, segment = c(0,1), ...)    {


  if (segment [1] != 0)
    From <- segment[1] * to + (1-segment[1]) * from
  else
    From <- from
    
  if (segment [2] != 1)
    To <- segment[2] * to + (1-segment[2]) * from
  else
    To <- to

  meanpi <- arr.pos*To+(1-arr.pos)*From
  if (endhead) To <- meanpi

  segments(From[1], From[2], To[1], To[2], lwd=lwd, lty=lty, col=lcol)

  mid2  <- c((arr.pos-0.01)*to+(1-arr.pos+0.01)*from)
  mid1  <- c((arr.pos     )*to+(1-arr.pos     )*from)

  Arrows(mid2[1], mid2[2], mid1[1], mid1[2], lcol=lcol, 
    arr.col=arr.col, ...)
  straightarrow <- mid1      # coordinates (x,y) where arrowhead is drawn

}
