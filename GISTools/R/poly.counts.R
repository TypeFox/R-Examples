

poly.counts<-function (pts, polys) 
  colSums(gContains(polys, pts, byid=TRUE))

poly.areas <- function(polys) gArea(polys, byid=TRUE)

generalize.polys <- function(sp,tol) gSimplify(sp,tol)


