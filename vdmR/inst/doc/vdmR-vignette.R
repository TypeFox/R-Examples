## ----setup, results="hide"-----------------------------------------------
library(vdmR)

## ------------------------------------------------------------------------
data(vsfuk2012)
head(vsfuk2012[,1:5])

## ----warning=FALSE, message=FALSE----------------------------------------
vscat(MortalityRate, FertilityRate, vsfuk2012, "scat01", "vsfuk2012")
vhist(MarriageRate, vsfuk2012, "hist01", "vsfuk2012")

## ----eval=FALSE----------------------------------------------------------
#  vlaunch(vsfuk2012, "main", "vsfuk2012")

## ------------------------------------------------------------------------
vscat(MortalityRate, FertilityRate, vsfuk2012, "scat01", "vsfuk2012",
      color=Type, size=pop_male)

## ------------------------------------------------------------------------
vscat(MortalityRate, FertilityRate, vsfuk2012, "scat01", "vsfuk2012",
      color=I("darkgreen"), size=pop_male)

## ----message=FALSE-------------------------------------------------------
vhist(MarriageRate, vsfuk2012, "hist01", "vsfuk2012",
      fill=I("darkgreen"), color=I("black"))

## ------------------------------------------------------------------------
vpcp(vsfuk2012, 4:17, "pcp1", "vsfuk2012",
     groupColumn="Type", scale="uniminmax", missing="min10")

## ------------------------------------------------------------------------
library(maptools)
shp.path <- file.path(system.file(package="vdmR"), "etc/shapes/fukuoka2012.shp")
vsfuk2012.spdf <- readShapeSpatial(shp.path, IDvar="CityCode")
head(vsfuk2012.spdf@data)

## ------------------------------------------------------------------------
frcol <- ggplot2::scale_fill_gradient2(low="blue", mid="white", high="red",
                              midpoint=median(vsfuk2012$FertilityRate))
vcmap(shp.path, vsfuk2012, "CityCode", "CityCode", "map1", "vsfuk2012",
      fill=FertilityRate, ggscale=frcol)

