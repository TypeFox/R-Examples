library(geotopbricks)


tz <-  "Etc/GMT+1"
start <- as.POSIXct("2002-03-22",tz=tz)
end <- as.POSIXct("2002-03-25",tz=tz)
day <- 3600*24
when <- seq(from=start,to=end,by=day)
#' # The data containing in the link are only for educational use
wpath <- "http://www.rendena100.eu/public/geotopbricks/simulations/idroclim_test1"
x <- "SoilLiqContentTensorFile"
when <- as.POSIXct("2002-03-22 UTC",tz="A")

# Not Run because it elapses too long time!!!
# Please Uncomment the following lines to run by yourself!!!


# wpath <- '/Users/ecor/attivita/2013/fem-idroclima/Trentino_500_dstr_GEOtop_1_225_9_002'

#kpsi <- "SoilLiqWaterPressTensorFile" ## soil water pressure head

out <-listFromOutputSoil3DTensor(x,when=when,wpath=wpath,tz=tz,use.read.raster.from.url=FALSE)