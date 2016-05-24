## ----setup, echo=FALSE---------------------------------------------------
library(knitr)
opts_chunk$set(dev="png", fig.path='figs/')
library(mptools)

## ----mpfile--------------------------------------------------------------
mp <- system.file('example.mp', package='mptools')

## ----res, message=FALSE--------------------------------------------------
res <- results(mp=mp)
res

## ----res_head------------------------------------------------------------
head(res$results[,, 'ALL'])

## ----res_head40----------------------------------------------------------
head(res$results[,, 'Pop 40'])

## ----meta, message=FALSE-------------------------------------------------
met <- meta(mp=mp)
head(met)

## ----xy, fig.width=6.6, fig.height=4.7-----------------------------------
xy <- mp2xy(mp=mp, r=habitat, cell.length=9.975)
head(xy)

## ----mp2sp, message=FALSE------------------------------------------------
spdf <- mp2sp(mp=mp, coords=xy, start=2000, 
              s_p4s='+init=epsg:3577', t_p4s='+init=epsg:4326') 

## ----plot_spplot, fig.width=6.6, fig.height=3.6--------------------------
library(sp)
library(viridis)
spplot(subset(spdf, time==2000), zcol='N',  
       cuts=c(-Inf, 0, 10^(0:6)), key.space='right', 
       col.regions=c('gray80', viridis(100)))

## ----write_shp, message=FALSE--------------------------------------------
library(rgdal)
writeOGR(spdf, dsn=tempdir(), layer='mp', driver='ESRI Shapefile', 
         overwrite_layer=TRUE)

## ----write_kml-----------------------------------------------------------
library(rgdal)
writeOGR(subset(spdf, time==2000), dsn=file.path(tempdir(), 'mp.kml'), 
         layer='', driver='KML')

## ----show_kml, eval=FALSE------------------------------------------------
#  file.show(file.path(tempdir(), 'mp.kml'))

## ----kch-----------------------------------------------------------------
k <- kch(meta=met, path=dirname(mp))
str(k)

## ----knt, message=FALSE, fig.width=6.6, fig.height=6.6-------------------
knt(meta=met, kch=k, pops=c('Pop 169', 'Pop 170', 'Pop 174', 'Pop 175'), 
    show_N=TRUE, results=res, samelims=TRUE, layout=c(2, 2))

## ----animate, message=FALSE----------------------------------------------
library(raster)
spdf <- mp2sp(mp=mp, coords=xy, start=2000)
mp_animate(spdf, habitat=habitat, outfile='dynamics.gif', zlim=c(0, 800), 
           width=630, height=615, overwrite=TRUE)

