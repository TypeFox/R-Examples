missingvalsexclude <-
function (rst, dat) 
{
  #NAvalue(rst) <- -9999
  rst@file@nodatavalue<- -9999
  fieldsmissing(dat, fields = c("ID", "x", "y", "Species", 
                                "Exclude","Reason"))
  spp <- unique(dat$Species)
  xy <- SpatialPoints(data.frame(x = dat$x, y = dat$y))
  cid <- cellFromXY(rst, xy)
  e <- extract(rst, xy)
  f1 <- which(is.na(e))
  dat$Exclude[f1] <- 1
  dat$Reason <- as.character(dat$Reason)
  dat$Reason[f1] <- 'No raster values'
  return(dat)
}
