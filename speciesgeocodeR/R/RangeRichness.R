RangeRichness <- function(ra, limits = c(-180, 180, -90, 90), reso = 60){
  limits <- extent(limits)
  xmin <- limits[1]
  xmax <- limits[2]
  ymin <- limits[3]
  ymax <- limits[4]
  if (xmin * xmax > 0) {
    cols <- abs(abs(slot(limits, "xmax")) - abs(slot(limits, "xmin")))
  } else {
    cols <- abs(abs(slot(limits, "xmax")) + abs(slot(limits, "xmin")))
  }
  if (ymin * ymax > 0) {
    rows <- abs(abs(slot(limits, "ymax")) - abs(slot(limits, "ymin")))
  } else {
    rows <- abs(abs(slot(limits, "ymax")) + abs(slot(limits, "ymin")))
  }
  reso <- 60/reso
  rasto <- raster(limits, ncol = cols* reso, nrow = rows * reso, vals = 1)
  
#   rasto <- raster(xmn=limits[1], xmx=limits[2], ymn=limits[3], ymx=limits[4], ncol = ncol, nrow = nrow, vals = 1)
  
  rang.ras <- lapply(ra, function(x) rasterize(x, rasto))
  
  for(i in 1:length(rang.ras)){
    rang.ras[[i]][is.na(rang.ras[[i]])] <-0
  }
  
  rang.ras <- stack(rang.ras)
  div.ras <- sum(rang.ras)
  div.ras[div.ras == 0] <- NA
  return(div.ras)
}