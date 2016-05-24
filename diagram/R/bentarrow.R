
##==============================================================================
# bentarrorw: 2- segmented arrow between two points
##==============================================================================

bentarrow <- function(from, to, lwd=2, lty=1, lcol="black", arr.col=lcol, 
   arr.side=2, arr.pos=0.5, path = "H", ...)  {

  sarr <- function(p1,p2,drawarr) {
    if (drawarr)
      m1<<-rbind(m1, straightarrow (from=p1, to=p2, arr.pos=arr.pos, lwd=lwd,
                     lty=lty, lcol=lcol, arr.col=arr.col, ...))  else
    segments(p1[1], p1[2], p2[1], p2[2], lwd=lwd, lty=lty, col=lcol)
  }

  m1    <- NULL


  if (path == "H")   { # horizontal
    sarr(  from,                 c(to[1], from[2]   ), 1 %in% arr.side)
    sarr(c(to[1], from[2]   ),    to                 , 2 %in% arr.side)
  } else {
    sarr(  from,                  c(from[1], to[2]   ), 1 %in% arr.side)
    sarr(c(from[1], to[2]   ),    to                  , 2 %in% arr.side)
  }
  bentarrow <- m1             # coordinates (x,y) where arrowhead is drawn
}
