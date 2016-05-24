data(GPS, envir=environment())
data(CH1903, envir=environment())


ReadPoly <-
function(shapefile, proj = GPS, fix.holes = TRUE, ID = NULL)
{
  polys <-  readShapePoly(shapefile, proj4string = proj, IDvar = ID, repair=TRUE, delete_null_obj=TRUE)
  if (fix.holes) {
    polys <- fix.holes(polys)
  }
  return(polys)
}
