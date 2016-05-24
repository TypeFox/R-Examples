tri.dellens <- function(voronoi.obj, exceptions = NULL, inverse = FALSE) {
  ## Return a list of delaunay segment lengths.
  ## If exceptions is a list of site numbers (normally that produced by
  ## voronoi.findrejectsites), exclude the voronoi
  ## triangles associated with that site.  If the inverse flag is
  ## TRUE, just return the segment lengths associated with the
  ## sites in the exceptions list.

  ## TODO - maybe this should be under voronoi.dellens, rather than
  ## tri.dellens?

  check.exceptions <- (length(exceptions) >0)
  nsites <- voronoi.obj$tri$n
  ntri <- length(voronoi.obj$p1)
  dists <- matrix(0, nrow=nsites, ncol=nsites)

  rejtri <- logical(length = ntri) # all FALSE

  if (check.exceptions) {
    anyrej <- exceptions[cbind(voronoi.obj$p1, voronoi.obj$p2, voronoi.obj$p3)];
    dim(anyrej) <- c(ntri,3);
    ## rejtri[i] is true if the ith triangle should be rejected.
    rejtri <- apply(anyrej, 1, any)
    if (inverse) rejtri <- !rejtri
  } else {
    ## accept all delauanay triangles.
  }

  ## ps is an Nx3 array - each row is one Delaunay triangle, giving
  ## the three sites in that triangle.
  ps <- cbind(voronoi.obj$p1, voronoi.obj$p2, voronoi.obj$p3)
  ps <- ps[which(!rejtri),] #throw away those triangles (row) not required.

  ## Find the distances between all sites 1,2 in each valid triangle.
  ## Then store those distances in the dists matrix.
  d <- tri.vordist(voronoi.obj, ps[,1], ps[,2])
  dists[ t(apply(ps[,1:2],1,sort))] <- d;

  d <- tri.vordist(voronoi.obj, ps[,2], ps[,3])
  dists[ t(apply(ps[,2:3],1,sort))] <- d;

  d <- tri.vordist(voronoi.obj, ps[,3], ps[,1])
  dists[ t(apply(cbind(ps[,3],ps[,1]),1,sort))] <- d;

  dists[which(dists>0)]
  
    
}

tri.vordist <- function (vor, p1, p2) {
  ## Return the Euclidean distance between site p1 and p2.
  ## Helper function for tri.dellens.
  ##
  ## Testing tri.vordist
  ## can be called fine with multiple arguments
  ## e.g. return distance between pts 1,3 and between pts 2,5
  ## tri.vordist(vor, c(1,2), c(3,5)) 
  ## Can also calculate output by hand, e.g. for pts (2,5)
  ## sqrt((vor$tri$x[2] - vor$tri$x[5])^2 + (vor$tri$y[2] - vor$tri$y[5])^2)

  dx <- vor$tri$x[p1] - vor$tri$x[p2]
  dy <- vor$tri$y[p1] - vor$tri$y[p2]
  sqrt( (dx*dx) + (dy*dy))
}
