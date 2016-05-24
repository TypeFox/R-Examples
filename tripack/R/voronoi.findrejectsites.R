voronoi.findrejectsites <- function(voronoi.obj, xmin, xmax, ymin, ymax)
{
  ## Given a voronoi object, find the reject sites, i.e. those sites
  ## with one of their vertices outside the bounded rectangle given by
  ## (xmin,ymin) and (xm ax,ymax).
  ## Return a vector `rejects': site N is a reject iff rejects[i] is T.
  nsites <- length(voronoi.obj$tri$x)
  rejects <- logical(nsites)
  outsiders <- ((voronoi.obj$x > xmax) | (voronoi.obj$x < xmin) |
                (voronoi.obj$y > ymax) | (voronoi.obj$y < ymin))


  ## In the list below, each site could be rejected more than once.
  rejects[c(voronoi.obj$p1[outsiders], voronoi.obj$p2[outsiders], voronoi.obj$p3[outsiders])] <- TRUE;

  rejects
}
