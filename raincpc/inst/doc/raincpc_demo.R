## ----, echo = FALSE, message = FALSE-------------------------------------
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE)

## ----, message = FALSE---------------------------------------------------
library(raincpc)
library(SDMTools)
library(raster)
library(ggplot2)

## ------------------------------------------------------------------------
cpc_get_rawdata(2014, 7, 1, 2014, 7, 7) 

## ------------------------------------------------------------------------
rain1 <- cpc_read_rawdata(2014, 7, 1)
print(rain1)

## ------------------------------------------------------------------------
rain2 <- cpc_read_rawdata(2014, 7, 2)
rain3 <- cpc_read_rawdata(2014, 7, 3)
rain4 <- cpc_read_rawdata(2014, 7, 4)
rain5 <- cpc_read_rawdata(2014, 7, 5)
rain6 <- cpc_read_rawdata(2014, 7, 6)
rain7 <- cpc_read_rawdata(2014, 7, 7)

rain_tot <- rain1 + rain2 + rain3 + rain4 + rain5 + rain6 + rain7
print(rain_tot)

## ------------------------------------------------------------------------
plot(rain_tot, 
     breaks = c(0, 1, 90, 180, 270, 360),
     col = c("grey", "red", "green", "blue"), 
     main = "Rainfall (mm) July 1 - July 7, 2014")

## ------------------------------------------------------------------------
raster_ggplot <- function(rastx) {
  require(SDMTools)

  stopifnot(class(rastx) == "RasterLayer")
  
  gfx_data <- getXYcoords(rastx)
  # lats need to be flipped
  gfx_data <- expand.grid(lons = gfx_data$x, lats = rev(gfx_data$y), 
                            stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  gfx_data$rain <- rastx@data@values
  
  return (gfx_data)
}

## ------------------------------------------------------------------------
rain_gg <- raster_ggplot(rain_tot)

rain_gg$rain_chunks <- cut(rain_gg$rain, breaks = c(0, 1, 90, 180, 270, 360), 
                         include.lowest = TRUE)

gfx_gg <- ggplot(data = rain_gg)
gfx_gg <- gfx_gg + geom_raster(aes(lons, lats, fill = rain_chunks))
gfx_gg <- gfx_gg + scale_fill_manual(values = c("lightgrey", "grey", "red", "green", "blue"))
gfx_gg <- gfx_gg + theme(axis.text = element_blank(), axis.ticks = element_blank())
gfx_gg <- gfx_gg + labs(x = NULL, y = NULL, fill = "Rain (mm)")
gfx_gg <- gfx_gg + ggtitle("Global Rainfall July 1 - July 7, 2014")
  
plot(gfx_gg)


## ------------------------------------------------------------------------
lon_vals <- seq(100.75, 130.75, 0.5)
lat_vals <- seq(-10.75, 20.75, 0.5)
reg_box <- expand.grid(lons = lon_vals, lats = lat_vals, 
                            stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
reg_box$rain <- extract.data(reg_box, rain_tot)
reg_box$rain_chunks <- cut(reg_box$rain, breaks = c(0, 1, 90, 180, 270, 360), 
                         include.lowest = TRUE)

## ------------------------------------------------------------------------
gfx_gg <- ggplot(data = reg_box)
gfx_gg <- gfx_gg + geom_raster(aes(lons, lats, fill = rain_chunks))
gfx_gg <- gfx_gg + scale_fill_manual(values = c("lightgrey", "grey", "red", "green", "blue"))
gfx_gg <- gfx_gg + theme(axis.text = element_blank(), axis.ticks = element_blank())
gfx_gg <- gfx_gg + labs(x = NULL, y = NULL, fill = "Rain (mm)")
gfx_gg <- gfx_gg + ggtitle("Rainfall over Southeast Asia July 1 - July 7, 2014")
  
plot(gfx_gg)

## ------------------------------------------------------------------------
cpc_get_rawdata(2014, 7, 1, 2014, 7, 7, usa = TRUE) 
rain1 <- cpc_read_rawdata(2014, 7, 1, usa = TRUE)
print(rain1)

rain2 <- cpc_read_rawdata(2014, 7, 2, usa = TRUE)
rain3 <- cpc_read_rawdata(2014, 7, 3, usa = TRUE)
rain4 <- cpc_read_rawdata(2014, 7, 4, usa = TRUE)
rain5 <- cpc_read_rawdata(2014, 7, 5, usa = TRUE)
rain6 <- cpc_read_rawdata(2014, 7, 6, usa = TRUE)
rain7 <- cpc_read_rawdata(2014, 7, 7, usa = TRUE)

rain_tot <- rain1 + rain2 + rain3 + rain4 + rain5 + rain6 + rain7
print(rain_tot)

## ------------------------------------------------------------------------
rain_gg <- raster_ggplot(rain_tot)

rain_gg$rain_chunks <- cut(rain_gg$rain, breaks = c(0, 1, 45, 90, 135, 180), 
                         include.lowest = TRUE)

gfx_gg <- ggplot(data = rain_gg)
gfx_gg <- gfx_gg + geom_raster(aes(lons, lats, fill = rain_chunks))
gfx_gg <- gfx_gg + scale_fill_manual(values = c("lightgrey", "grey", "red", "green", "blue"))
gfx_gg <- gfx_gg + theme(axis.text = element_blank(), axis.ticks = element_blank())
gfx_gg <- gfx_gg + labs(x = NULL, y = NULL, fill = "Rain (mm)")
gfx_gg <- gfx_gg + ggtitle("USA Rainfall July 1 - July 7, 2014")
  
plot(gfx_gg)


## ------------------------------------------------------------------------
rain1 <- cpc_read_rawdata(2014, 7, 1, usa = TRUE, write_output = TRUE)

