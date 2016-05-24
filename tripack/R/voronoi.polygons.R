# from Denis White <white.denis@epamail.epa.gov> :

voronoi.polygons <- function (voronoi.obj)
{
  nsites <- length(voronoi.obj$tri$x)
  polys <- list()
  j <- 0
  for (i in 1:nsites) {
    vs <- voronoi.findvertices(i, voronoi.obj)
    if (length(vs) > 0) {
      j <- j + 1
      polys[[j]] <- cbind (x=voronoi.obj$x[vs],
                           y=voronoi.obj$y[vs])
    }
  }
  class(polys)<-"voronoi.polygons"
  polys
}



