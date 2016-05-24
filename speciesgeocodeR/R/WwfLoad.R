WwfLoad <-
function(x){
  if ( x == ""){x <- getwd()}
  download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
  destfile = file.path(x, "wwf_ecoregions.zip"))#?1349272619")
  unzip(file.path(x, "wwf_ecoregions.zip"), exdir = file.path(x, "WWF_ecoregions"))
  file.remove(file.path(x, "wwf_ecoregions.zip"))
  wwf <- maptools::readShapeSpatial(file.path(x, "WWF_ecoregions", "official", "wwf_terr_ecos.shp"))
  return(wwf)
}
