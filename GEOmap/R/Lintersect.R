`Lintersect` <-
function( l1,  l2)
  {
    ######   l1 and l2 are lists (lines)  consisting of points p1 and p2
    ######  l1 = list(p1=list(x=0, y=0), p2=list(x=1,y=1))

  return(((ccw(l1$p1, l1$p2, l2$p1) * ccw(l1$p1, l1$p2, l2$p2)) <= 0)
	 && ((ccw(l2$p1, l2$p2, l1$p1) * ccw(l2$p1, l2$p2, l1$p2)) <= 0));
  }

