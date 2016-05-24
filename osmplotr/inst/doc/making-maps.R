## ----load, message=FALSE-------------------------------------------------
library (osmplotr)
library (maptools) # Needed for this vignette

## ---- echo=FALSE, message=FALSE------------------------------------------
setwd ('../..')
devtools::load_all ('osmplotr', export_all=FALSE)
setwd ('./osmplotr/vignettes')

## ---- echo=FALSE, message=FALSE------------------------------------------
# Combining (dat_B, dat_BC) and (dat_H, dat_HP) requires removing the repeated
# objects
indx <- which (!london$dat_BR$id %in% london$dat_BNR$id)
dat_B <- maptools::spRbind (london$dat_BR [indx,], london$dat_BNR)
indx <- which (!london$dat_H$id %in% london$dat_HP$id)
dat_H <- maptools::spRbind (london$dat_H [indx,], london$dat_HP)
dat_T <- london$dat_T

## ------------------------------------------------------------------------
bbox <- get_bbox (c(-0.13,51.50,-0.11,51.52))

## ---- eval=FALSE---------------------------------------------------------
#  dat_B <- extract_osm_objects (key='building', bbox=bbox)

## ----map1----------------------------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')

## ------------------------------------------------------------------------
map <- add_osm_objects (map, dat_B, col='gray40')

## ---- eval=FALSE---------------------------------------------------------
#  print (map)

## ----map1-print, echo=FALSE----------------------------------------------
mapwd <- 500
mapht <- mapwd * diff (bbox [2,]) / diff (bbox [1,])
png (height=mapht, width=mapwd, file='map_a1.png')
print (map)
graphics.off ()

## ---- eval=FALSE---------------------------------------------------------
#  bbox <- get_bbox (c(-0.13,51.51,-0.11,51.52))
#  dat_B <- extract_osm_objects (key='building', bbox=bbox)
#  dat_H <- extract_osm_objects (key='highway', bbox=bbox)

## ------------------------------------------------------------------------
class (dat_B); class (dat_H); class (dat_T)

## ------------------------------------------------------------------------
length (dat_B); length (dat_H); length (dat_T)

## ---- eval=FALSE---------------------------------------------------------
#  dat_T <- extract_osm_objects (key='natural', value='tree', bbox=bbox)

## ---- eval=FALSE---------------------------------------------------------
#  dat_NT <- extract_osm_objects (key='natural', value='!tree', bbox=bbox)

## ---- eval=FALSE---------------------------------------------------------
#  dat_H <- extract_osm_objects (key='highway', value='!primary', bbox=bbox)

## ---- eval=FALSE---------------------------------------------------------
#  extra_pairs <- c ('name', 'Royal.Festival.Hall')
#  dat <- extract_osm_objects (key='building', extra_pairs=extra_pairs,
#                                         bbox=bbox)

## ---- eval=FALSE---------------------------------------------------------
#  extra_pairs <- list (c ('addr:street', 'Stamford.St'),
#                       c ('addr:housenumber', '150'))
#  dat <- extract_osm_objects (key='building', extra_pairs=extra_pairs,
#                                        bbox=bbox)

## ------------------------------------------------------------------------
osm_structures ()

## ------------------------------------------------------------------------
osm_structures()$value [1:4]

## ---- eval=FALSE---------------------------------------------------------
#  dat <- make_osm_map (structures=osm_structures (), bbox=bbox)

## ---- echo=FALSE---------------------------------------------------------
dat1 <- list (dat_BU=NULL, dat_A=NULL, dat_W=NULL, dat_G=NULL,
              dat_N=NULL, dat_P=NULL, dat_H=NULL, dat_BO=NULL, dat_T=NULL)
dat <- list (osm_data=dat1, map=ggplot2::ggplot ())

## ------------------------------------------------------------------------
names (dat); sapply (dat, class); names (dat$osm_data)

## ------------------------------------------------------------------------
osm_structures (structures=c('building', 'highway'))

## ------------------------------------------------------------------------
osm_structures (structures='grass')

## ------------------------------------------------------------------------
structures <- c ('highway', 'highway', 'building', 'building', 'building',
                 'amenity', 'park', 'natural', 'tree')   
structs <- osm_structures (structures=structures, col_scheme='dark')   
structs$value [1] <- '!primary'   
structs$value [2] <- 'primary'
structs$suffix [2] <- 'HP'
structs$value [3] <- '!residential'
structs$value [4] <- 'residential'
structs$value [5] <- 'commercial'
structs$suffix [3] <- 'BNR'
structs$suffix [4] <- 'BR'
structs$suffix [5] <- 'BC'

## ---- eval=FALSE---------------------------------------------------------
#  london <- make_osm_map (structures=structs, bbox=bbox)
#  london <- london$osm_data

## ---- echo=FALSE---------------------------------------------------------
highways1 <- london$highways1

## ---- eval=FALSE---------------------------------------------------------
#  highways <- c ('Monmouth.St', 'Short.?s.Gardens', 'Endell.St', 'Long.Acre',
#                 'Upper.Saint.Martin')
#  highways1 <- connect_highways (highways=highways, bbox=bbox)

## ------------------------------------------------------------------------
class (highways1); length (highways1); head (coordinates (highways1))

## ---- eval=FALSE---------------------------------------------------------
#  highways <- c ('Endell.St', 'High.Holborn', 'Drury.Lane', 'Long.Acre')
#  highways2 <- connect_highways (highways=highways, bbox=bbox)
#  highways <- c ('Drury.Lane', 'High.Holborn', 'Kingsway', 'Great.Queen.St')
#  highways3 <- connect_highways (highways=highways, bbox=bbox)

## ----connect_highways, fig.width=4, message=FALSE, eval=TRUE-------------
# TODO: Set eval=FALSE before cran re-sub!!
bbox_big <- get_bbox (c(-0.15,51.5,-0.10,51.52))
highways <- c ('Kingsway', 'Holborn', 'Farringdon.St', 'Strand',
               'Fleet.St', 'Aldwych')
highway_list <- connect_highways (highways=highways, bbox=bbox_big, plot=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  dat_B <- extract_osm_objects (key='building', bbox=bbox)
#  dat_H <- extract_osm_objects (key='highway', bbox=bbox)
#  dat_T <- extract_osm_objects (key='natural', value='tree', bbox=bbox)

## ----map2----------------------------------------------------------------
bbox_small <- get_bbox (c(-0.13,51.51,-0.11,51.52))
map <- plot_osm_basemap (bbox=bbox_small, bg='gray20')
map <- add_osm_objects (map, dat_H, col='gray70')
map <- add_osm_objects (map, dat_B, col='gray40')

## ---- fig.width=4, fig.height=2, eval=FALSE------------------------------
#  dev.new (width=8, height=8 * diff (bbox2 [2,]) / diff (bbox2 [1,]))
#  # or png, pdf, or whatever device is used for printing
#  print (map)

## ----map2-print, echo=FALSE----------------------------------------------
mapht <- mapwd * diff (bbox_small [2,]) / diff (bbox_small [1,])
png (height=mapht, width=mapwd, file='map_a2.png')
print (map)
graphics.off ()

## ----map3, eval=FALSE----------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox_small, bg='gray20')
#  map <- add_osm_objects (map, dat_B, col='gray40', border='orange', size=0.2)
#  print (map)

## ----map3-print, echo=FALSE----------------------------------------------
map <- plot_osm_basemap (bbox=bbox_small, bg='gray20')
map <- add_osm_objects (map, dat_B, col='gray40', border='orange', size=0.2)
png (height=mapht, width=mapwd, file='map_a3.png')
print (map)
graphics.off ()

## ----map4, eval=FALSE----------------------------------------------------
#  map <- add_osm_objects (map, dat_H, col='gray70', size=0.7)
#  map <- add_osm_objects (map, dat_T, col='green', size=2, shape=1)
#  print (map)

## ----map4-print, echo=FALSE----------------------------------------------
map <- add_osm_objects (map, dat_H, col='gray70', size=0.7)
map <- add_osm_objects (map, dat_T, col='green', size=2, shape=1)
png (height=mapht, width=mapwd, file='map_a4.png')
print (map)
graphics.off ()

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  png (height=mapht, width=mapwd, file='map.png')
#  print (map)
#  graphics.off ()

## ---- eval=FALSE---------------------------------------------------------
#  dat_HP <- extract_osm_objects (key='highway', value='primary', bbox=bbox)
#  dat_H <- extract_osm_objects (key='highway', value='!primary', bbox=bbox)

## ---- echo=FALSE---------------------------------------------------------
dat_HP <- london$dat_HP
dat_H <- london$dat_H

## ----map5, eval=FALSE----------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox_small, bg='gray20')
#  map <- add_osm_objects (map, dat_H, col='gray50')
#  map <- add_osm_objects (map, dat_HP, col='gray80', size=2)
#  print (map)

## ----map5-print, echo=FALSE----------------------------------------------
map <- plot_osm_basemap (bbox=bbox_small, bg='gray20')
map <- add_osm_objects (map, dat_H, col='gray50')
map <- add_osm_objects (map, dat_HP, col='gray80', size=2)
png (height=mapht, width=mapwd, file='map_a5.png')
print (map)
graphics.off ()

## ---- echo=FALSE---------------------------------------------------------
dat_RFH <- london$dat_RFH
dat_ST <- london$dat_ST

## ----map7, eval=FALSE----------------------------------------------------
#  bbox_small2 <- get_bbox (c (-0.118, 51.504, -0.110, 51.507))
#  map <- plot_osm_basemap (bbox=bbox_small2, bg='gray95')
#  map <- add_osm_objects (map, dat_H, col='gray80')
#  map <- add_osm_objects (map, dat_HP, col='gray60', size=2)
#  map <- add_osm_objects (map, dat_RFH, col='orange', border='red', size=2)
#  map <- add_osm_objects (map, dat_ST, col='skyblue', border='blue', size=2)
#  dev.new (width=8, height=8 * diff (bbox_small2 [2,]) / diff (bbox_small2 [1,]))
#  print (map)

## ----map7-print, echo=FALSE----------------------------------------------
bbox_small2 <- get_bbox (c (-0.118, 51.504, -0.110, 51.507))
map <- plot_osm_basemap (bbox=bbox_small2, bg='gray95')
map <- add_osm_objects (map, dat_H, col='gray80')
map <- add_osm_objects (map, dat_HP, col='gray60', size=2)
map <- add_osm_objects (map, dat_RFH, col='orange', border='red', size=2)
map <- add_osm_objects (map, dat_ST, col='skyblue', border='blue', size=2)
mapht <- mapwd * diff (bbox_small2 [2,]) / diff (bbox_small2 [1,])
png (height=mapht, width=mapwd, file='map_a7.png')
print (map)
graphics.off ()

## ----map8, eval=FALSE----------------------------------------------------
#  structs <- c ('highway', 'building', 'park', 'grass', 'tree')
#  structures <- osm_structures (structures=structs, col_scheme='light')
#  dat <- make_osm_map (structures=structures, bbox=bbox)
#  map <- dat$map
#  dev.new (width=8, height=8 * diff (bbox [2,]) / diff (bbox [1,]))
#  print (map)

## ----map8-print, echo=FALSE----------------------------------------------
structs <- c ('highway', 'building', 'amenity', 'park', 'grass', 'tree')   
structures <- osm_structures (structures=structs, col_scheme='light')   
osm_dat <- list (dat_B=dat_B, dat_H=dat_H, dat_P=london$dat_P,
                 dat_A=london$dat_A, dat_G=london$dat_G, dat_T=london$dat_T)
dat <- make_osm_map (structures=structures, osm_data=osm_dat)
map <- dat$map
mapht <- mapwd * diff (bbox [2,]) / diff (bbox [1,])
png (height=mapht, width=mapwd, file='map_a8.png')
print (map)
graphics.off ()

## ----map9, eval=FALSE----------------------------------------------------
#  dat <- make_osm_map (osm_data=dat$osm_data, structures=structures,
#                       bbox=bbox_small)
#  print (dat$map)

## ----map9-print, echo=FALSE----------------------------------------------
dat <- make_osm_map (structures=structures, osm_data=osm_dat, bbox=bbox_small)
mapht <- mapwd * diff (bbox_small [2,]) / diff (bbox_small [1,])
png (height=mapht, width=mapwd, file='map_a9.png')
print (dat$map)
graphics.off ()

## ------------------------------------------------------------------------
structs <- c ('amenity', 'building', 'grass', 'highway', 'park')
osm_structures (structs, col_scheme='light')

## ----map10---------------------------------------------------------------
structures <- osm_structures (structures=structs, col_scheme='dark')   
dat <- make_osm_map (structures=structures, osm_data=dat$osm_dat, bbox=bbox)
map <- add_axes (dat$map, colour='black')

## ---- eval=FALSE---------------------------------------------------------
#  print (map)

## ----map10-print, echo=FALSE---------------------------------------------
png (height=mapht, width=mapwd, file='map_a10.png')
print (map)
graphics.off ()

## ----map11, eval=FALSE---------------------------------------------------
#  map <- add_axes (map, colour='blue', pos=c(0.1,0.2))
#  print (map)

## ----map11-print, echo=FALSE---------------------------------------------
map <- add_axes (map, colour='blue', pos=c(0.1,0.2))
png (height=mapht, width=mapwd, file='map_a11.png')
print (map)
graphics.off ()

