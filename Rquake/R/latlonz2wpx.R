latlonz2wpx<-function(twpx, stas)
{
  ####  populate a pick list with the lat-lon-z from a station list
  msta = match(twpx$name, stas$name)

  latz = stas$lat[msta]
  lonz = stas$lon[msta]
  zees =stas$z[msta]

  twpx$lat = latz
  twpx$lon = lonz
  twpx$z = zees
  return(twpx)

}
