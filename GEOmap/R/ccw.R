`ccw` <-
function(p0, p1, p2)
  {
  dx1 = p1$x - p0$x;
  dy1 = p1$y - p0$y;
  dx2 = p2$x - p0$x;
  dy2 = p2$y - p0$y;

  if(dx1 * dy2 > dy1 * dx2) return(1);
  if(dx1 * dy2 < dy1 * dx2) return(-1);

  if((dx1 * dx2 < 0) || (dy1 * dy2 < 0)) return(-1);
  if((dx1^2 + dy1^2) < (dx2^2 + dy2^2)) return(1);
  return(0);
  }

