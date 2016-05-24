plotProbcast <-
function(forecast, longitude, latitude, nGrid = 65, 
         type = c("image", "contour", "persp"), ..., 
         interpolate = FALSE, span = 0.75, maps = NULL)
{
## col=rev(rainbow(100,start=0,end=0.85))

# FIELDS <- exists("US", mode = "function") 
# FIELDS <- FIELDS && exists("world", mode = "function")
# FIELDS <- FIELDS && exists("image.plot", mode = "function")

# MAPS <- exists("map", mode = "function") 
# MAPS <- MAPS & exists("area.map", mode = "function") 
# MAPS <- MAPS & exists("match.map", mode = "function") 

 type <- type[1]

 lonRange <- range(longitude)
 lon <- seq(from=lonRange[1],to=lonRange[2],length=nGrid)

 latRange <- range(latitude)
 lat <- seq(from=latRange[1],to=latRange[2],length=nGrid)

 if (interpolate) {
# use prediction on a grid from loess fit

   fit <- loess( forecast ~ longitude*latitude, span = span)

   pred <- predict( fit, expand.grid(list( longitude = lon, latitude = lat)))
  }
  else {
# use prediction on a grid via binning

    pred <- binGrid( forecast, longitude, latitude, nGrid)
  }

   if (type == "image") {
##     fields:::image.plot(lon, lat, pred, xlim=lonRange, ylim=latRange,
##     if (!exists("image.plot")) library(fields)
       fields::image.plot(lon, lat, pred, xlim=lonRange, ylim=latRange,
                  xlab = "", ylab = "", horizontal = TRUE, ...)
  }
  else  {
      do.call( type, c(list(x = lon, y = lat, 
               z = pred, xlab = "", ylab = "", ...)))
  }
    

if (is.null(maps)) maps <- type == "image"| type == "contour"

if (maps) {

if(min(lon) <= -124.7 & max(lon) >= -124.7){
  US.map <- 1
}
if(min(lon) >= -124.7 & max(lon) <= -67.1){
  US.map <- 1
}
if(min(lon) <= -67.1 & max(lon) >= -67.1){
  US.map <- 1
}
if(min(lat) <= 25.2 & max(lat) >= 25.2){
  US.map <- 1 + US.map
}
if(min(lat) >= 25.2 & max(lat) <= 49.4){
  US.map <- 1 + US.map
}
if(min(lat) <= 49.4 & max(lat) >= 49.4){
  US.map <- 1 + US.map
}

lonRange <- range(lon)
latRange <- range(lat)

if(US.map==2){
##  fields:::US( xlim = lonRange, ylim = latRange,  add=TRUE, col=1, lwd=2)
##  if (!exists("US")) library("fields")
    fields::US( xlim = lonRange, ylim = latRange,  add=TRUE, col=1, lwd=2)
##  maps:::map('world', 'Canada', interior = TRUE, 
##   if (!exists("map")) library(maps)
    maps::map('world', 'Canada', interior = TRUE, 
           xlim = lonRange, ylim = latRange, add=TRUE, col=1, lwd=2 )
}
 else {
##   fields:::world( xlim=lonRange, ylim=latRange, add=TRUE, col=1, lwd=2)     
##   if (!exists("world")) library("fields")
     fields::world( xlim=lonRange, ylim=latRange, add=TRUE, col=1, lwd=2)     
}
}  

invisible()
}

