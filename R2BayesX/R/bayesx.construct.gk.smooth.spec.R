bayesx.construct.gk.smooth.spec <- bayesx.construct.geokriging.smooth.spec <- function(object, dir, prg, data)
{
  return(geo.smooth.spec(object, dir, prg, data, "geokriging"))
}

