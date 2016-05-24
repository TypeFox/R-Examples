## ----load, message=FALSE-------------------------------------------------
library (maptools) # Needed for this vignette

## ---- echo=FALSE, message=FALSE------------------------------------------
setwd ('../..')
devtools::load_all ("osmplotr", export_all=FALSE)
setwd ('./osmplotr/vignettes')

## ------------------------------------------------------------------------
bbox <- get_bbox (c(-0.13,51.51,-0.11,51.52))

## ---- echo=FALSE---------------------------------------------------------
mapwd <- 500
mapht <- mapwd * diff (bbox [2,]) / diff (bbox [1,])

## ---- echo=FALSE, message=FALSE------------------------------------------
indx <- which (!london$dat_BR$id %in% london$dat_BNR$id)
dat_B <- maptools::spRbind (london$dat_BR [indx,], london$dat_BNR)
indx <- which (!london$dat_H$id %in% london$dat_HP$id)
dat_H <- maptools::spRbind (london$dat_H [indx,], london$dat_HP)
dat_HP <- london$dat_HP

## ---- eval=FALSE---------------------------------------------------------
#  dat_B <- extract_osm_objects (key='building', bbox=bbox)

## ----map1, eval=FALSE----------------------------------------------------
#  pts <- sp::SpatialPoints (cbind (c (-0.115, -0.125, -0.125, -0.115),
#                               c (51.513, 51.513, 51.517, 51.517)))
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_groups (map, dat_B, groups=pts, cols='orange', bg='gray40')
#  print (map)

## ----map1-print, echo=FALSE----------------------------------------------
pts <- sp::SpatialPoints (cbind (c (-0.115, -0.125, -0.125, -0.115),
                             c (51.513, 51.513, 51.517, 51.517)))
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_groups (map, dat_B, groups=pts, cols='orange', bg='gray40')
png (height=mapht, width=mapwd, file='map_b1.png')
print (map)
graphics.off ()

## ----map2, eval=FALSE----------------------------------------------------
#  pts2 <- sp::SpatialPoints (cbind (c (-0.111, -0.1145, -0.1145, -0.111),
#                               c (51.517, 51.517, 51.519, 51.519)))
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_groups (map, dat_B, groups=list (pts, pts2),
#                         cols=c ('orange', 'tomato'), bg='gray40')
#  print (map)

## ----map2-print, echo=FALSE----------------------------------------------
pts2 <- sp::SpatialPoints (cbind (c (-0.111, -0.1145, -0.1145, -0.111),
                             c (51.517, 51.517, 51.519, 51.519)))
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_groups (map, dat_B, groups=list (pts, pts2), 
                       cols=c ('orange', 'tomato'), bg='gray40')

png (height=mapht, width=mapwd, file='map_b2.png')
print (map)
graphics.off ()

## ----map3, eval=FALSE----------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_groups (map, dat_B, groups=list (pts, pts2),
#                         cols=c ('orange', 'tomato'))
#  print (map)

## ----map3-print, echo=FALSE----------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_groups (map, dat_B, groups=list (pts, pts2), 
                       cols=c ('orange', 'tomato'))
png (height=mapht, width=mapwd, file='map_b3.png')
print (map)
graphics.off ()

## ----map5, eval=FALSE----------------------------------------------------
#  pts <- sp::SpatialPoints (rbind (coordinates (pts), c (-0.12, 51.515)))
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_groups (map, dat_B, groups=pts, cols='orange', bg='gray40')
#  print (map)

## ----map5-print, echo=FALSE----------------------------------------------
pts <- sp::SpatialPoints (rbind (coordinates (pts), c (-0.12, 51.515)))
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_groups (map, dat_B, groups=pts, cols='orange', bg='gray40')
png (height=mapht, width=mapwd, file='map_b5.png')
print (map)
graphics.off ()

## ----map6, eval=FALSE----------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_groups (map, dat_B, groups=list (pts, pts2), make_hull=TRUE,
#                         cols=c ('orange', 'tomato'), bg='gray40', boundary=1)
#  print (map)

## ----map6-print, echo=FALSE----------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_groups (map, dat_B, groups=list (pts, pts2), make_hull=TRUE,
                       cols=c ('orange', 'tomato'), bg='gray40', boundary=1)
png (height=mapht, width=mapwd, file='map_b6.png')
print (map)
graphics.off ()

## ----map7, eval=FALSE----------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_groups (map, dat_B, groups=list (pts, pts2), make_hull=TRUE,
#                         cols=c ('orange', 'tomato'), bg='gray40', boundary=0)
#  print (map)

## ----map7-print, echo=FALSE----------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_groups (map, dat_B, groups=list (pts, pts2), make_hull=TRUE,
                       cols=c ('orange', 'tomato'), bg='gray40', boundary=0)
png (height=mapht, width=mapwd, file='map_b7.png')
print (map)
graphics.off ()

## ---- eval=FALSE---------------------------------------------------------
#  dat_P <- extract_osm_objects (key='park', bbox=bbox)

## ------------------------------------------------------------------------
osm_structures (structure='park')

## ---- echo=FALSE---------------------------------------------------------
dat_P <- london$dat_P

## ----map8, eval=FALSE----------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_groups (map, dat_B, groups=list (pts, pts2), make_hull=TRUE,
#                         cols=c ('orange', 'tomato'), bg='gray40', boundary=0)
#  col_park_in <- rgb (50, 255, 50, maxColorValue=255)
#  col_park_out <- rgb (50, 155, 50, maxColorValue=255)
#  map <- add_osm_groups (map, dat_P, groups=list (pts, pts2),
#                         cols=rep (col_park_in, 2), bg=col_park_out, boundary=0)
#  print (map)

## ----map8-print, echo=FALSE----------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_groups (map, dat_B, groups=list (pts, pts2), make_hull=TRUE,
                       cols=c ('orange', 'tomato'), bg='gray40', boundary=0)
col_park_in <- rgb (50, 255, 50, maxColorValue=255)
col_park_out <- rgb (50, 155, 50, maxColorValue=255)
map <- add_osm_groups (map, dat_P, groups=list (pts, pts2), 
                       cols=rep (col_park_in, 2), bg=col_park_out, boundary=0)
png (height=mapht, width=mapwd, file='map_b8.png')
print (map)
graphics.off ()

## ----map9, eval=FALSE----------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_objects (map, dat_P, col=col_park_out)
#  map <- add_osm_groups (map, dat_P, groups=list (pts, pts2),
#                         cols=rep (col_park_in, 2), bg=col_park_out, boundary=0)
#  map <- add_osm_groups (map, dat_B, groups=list (pts, pts2), make_hull=TRUE,
#                         cols=c ('orange', 'tomato'), bg='gray40', boundary=0)
#  print (map)

## ----map9-print, echo=FALSE----------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_objects (map, dat_P, col=col_park_out)
map <- add_osm_groups (map, dat_P, groups=list (pts, pts2), 
                       cols=rep (col_park_in, 2), bg=col_park_out, boundary=0)
map <- add_osm_groups (map, dat_B, groups=list (pts, pts2), make_hull=TRUE,
                       cols=c ('orange', 'tomato'), bg='gray40', boundary=0)
png (height=mapht, width=mapwd, file='map_b9.png')
print (map)
graphics.off ()

## ----map10, eval=FALSE---------------------------------------------------
#  dat_HP <- extract_osm_objects (key="highway", value="primary", bbox=bbox)
#  dat_H <- extract_osm_objects (key="highway", value="!primary", bbox=bbox)
#  # darken colours by aboud 20%
#  cols_adj <- adjust_colours (c ('orange', 'tomato'), adj=-0.2)
#  map <- add_osm_groups (map, dat_HP, groups=list (pts, pts2), make_hull=TRUE,
#                         cols=cols_adj, bg=adjust_colours ('gray40', adj=-0.4),
#                         boundary=1, size=2)
#  map <- add_osm_groups (map, dat_H, groups=list (pts, pts2), make_hull=TRUE,
#                         cols=cols_adj, bg=adjust_colours ('gray40', adj=-0.2),
#                         boundary=1, size=1)
#  print (map)

## ----map10-print, echo=FALSE---------------------------------------------
dat_H <- london$dat_H
dat_HP <- london$dat_HP
cols_adj <- adjust_colours (c ('orange', 'tomato'), adj=-0.2)
map <- add_osm_groups (map, dat_HP, groups=list (pts, pts2), make_hull=TRUE,
                       cols=cols_adj, bg=adjust_colours ('gray40', adj=-0.4), 
                       boundary=1, size=2)
map <- add_osm_groups (map, dat_H, groups=list (pts, pts2), make_hull=TRUE,
                       cols=cols_adj, bg=adjust_colours ('gray40', adj=-0.2), 
                       boundary=1, size=1)
png (height=mapht, width=mapwd, file='map_b10.png')
print (map)
graphics.off ()

## ----map11, eval=FALSE---------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray95')
#  map <- add_osm_groups (map, dat_B, groups=pts, cols='gray40', bg='gray85',
#                     boundary=1)
#  map <- add_osm_groups (map, dat_H, groups=pts, cols='gray20', bg='gray70',
#                     boundary=0)
#  map <- add_osm_groups (map, dat_HP, groups=pts, cols='gray10', bg='white',
#                     boundary=0, size=1)
#  print (map)

## ----map11-print, echo=FALSE---------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray95')
map <- add_osm_groups (map, dat_B, groups=pts, cols='gray40', bg='gray85',
                   boundary=1)
map <- add_osm_groups (map, dat_H, groups=pts, cols='gray20', bg='gray70',
                   boundary=0)
map <- add_osm_groups (map, dat_HP, groups=pts, cols='gray10', bg='white',
                   boundary=0, size=1)
png (height=mapht, width=mapwd, file='map_b11.png')
print (map)
graphics.off ()

## ------------------------------------------------------------------------
set.seed (2)
ngroups <- 12
x <- bbox [1,1] + runif (ngroups) * diff (bbox [1,])
y <- bbox [2,1] + runif (ngroups) * diff (bbox [2,])
groups <- cbind (x, y)
groups <- apply (groups, 1, function (i) 
              sp::SpatialPoints (matrix (i, nrow=1, ncol=2)))

## ----map12, eval=FALSE---------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray95')
#  map <- add_osm_groups (map, dat_B, groups=groups,
#                         cols=rainbow (length (groups)))
#  print (map)

## ----map12-print, echo=FALSE---------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray95')
map <- add_osm_groups (map, dat_B, groups=groups, 
                       cols=rainbow (length (groups)))
png (height=mapht, width=mapwd, file='map_b12.png')
print (map)
graphics.off ()

## ----map13, eval=FALSE---------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray95')
#  map <- add_osm_groups (map, dat_B, groups=groups, borderWidth=2,
#                         cols=heat.colors (length (groups)))
#  print (map)

## ----map13-print, echo=FALSE---------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray95')
map <- add_osm_groups (map, dat_B, groups=groups, borderWidth=2,
                       cols=heat.colors (length (groups)))
png (height=mapht, width=mapwd, file='map_b13.png')
print (map)
graphics.off ()

## ---- fig.width=4--------------------------------------------------------
plot.new ()
cmat <- colour_mat (plot=TRUE)

## ---- fig.width=4--------------------------------------------------------
plot.new ()
cmat <- colour_mat (n=c(4,8), rotate=90, plot=TRUE)

## ----map14, eval=FALSE---------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray95')
#  map <- add_osm_groups (map, dat_B, groups=groups, borderWidth=2, colmat=TRUE,
#                         cols=c('red','green','yellow','blue'), rotate=180)
#  print (map)

## ----map14-print, echo=FALSE---------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray95')
map <- add_osm_groups (map, dat_B, groups=groups, borderWidth=2, colmat=TRUE,
                       cols=c('red','green','yellow','blue'), rotate=180)
png (height=mapht, width=mapwd, file='map_b14.png')
print (map)
graphics.off ()

## ----connect-highways, eval=FALSE----------------------------------------
#  highways <- c ('Monmouth.St', 'Short.?s.Gardens', 'Endell.St', 'Long.Acre',
#                 'Upper.Saint.Martin')
#  highways1 <- connect_highways (highways=highways, bbox=bbox)
#  highways <- c ('Endell.St', 'High.Holborn', 'Drury.Lane', 'Long.Acre')
#  highways2 <- connect_highways (highways=highways, bbox=bbox, plot=T)
#  highways <- c ('Drury.Lane', 'High.Holborn', 'Kingsway', 'Great.Queen.St')
#  highways3 <- connect_highways (highways=highways, bbox=bbox, plot=T)

## ---- echo=FALSE---------------------------------------------------------
highways1 <- london$highways1
highways2 <- london$highways2
highways3 <- london$highways3

## ------------------------------------------------------------------------
class (highways1); length (highways1); length (highways2); length (highways3)

## ---- echo=TRUE----------------------------------------------------------
groups <- list (london$highways1, london$highways2, london$highways3)
cols_B <- c ('red', 'orange', 'tomato') # for the 3 groups
cols_H <- adjust_colours (cols_B, -0.2)
bg_B <- 'gray40'
bg_H <- 'gray60'

## ----map15, eval=FALSE---------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_objects (map, dat_P, col=col_park_out)
#  map <- add_osm_groups (map, dat_B, groups=groups, boundary=1,
#                     bg=bg_B, col=cols_B)
#  map <- add_osm_groups (map, dat_H, groups=groups, boundary=1,
#                     bg=bg_H, col=cols_H)
#  map <- add_osm_groups (map, dat_HP, groups=groups, boundary=0,
#                     cols=cols_H, bg=bg_H, size=1)
#  print (map)

## ----map15-print, echo=FALSE---------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_objects (map, dat_P, col=col_park_out)
map <- add_osm_groups (map, dat_B, groups=groups, boundary=1,
                   cols=cols_B, bg=bg_B)
map <- add_osm_groups (map, dat_H, groups=groups, boundary=0,
                   cols=cols_H, bg=bg_H)
map <- add_osm_groups (map, dat_HP, groups=groups, boundary=0,
                   cols=cols_H, bg=bg_H, size=1)
png (height=mapht, width=mapwd, file='map_b15.png')
print (map)
graphics.off ()

## ------------------------------------------------------------------------
n <- 5
x <- seq (bbox [1,1], bbox [1,2], length.out=n)
y <- seq (bbox [2,1], bbox [2,2], length.out=n)
dat <- data.frame (
    x=as.vector (array (x, dim=c(n, n))),
    y=as.vector (t (array (y, dim=c(n, n)))),
    z=x * y
    )
head (dat)

## ----map16, eval=FALSE---------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_surface (map, dat_B, dat=dat, cols=heat.colors (30))
#  print (map)

## ----map16-print, echo=FALSE---------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_surface (map, dat_B, dat=dat, cols=heat.colors (30))
png (height=mapht, width=mapwd, file='map_b16.png')
print (map)
graphics.off ()

## ---- eval=TRUE----------------------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_surface (map, dat_HP, dat=dat, cols=heat.colors (30))
map <- add_osm_surface (map, dat_H, dat=dat, cols=heat.colors (30))

## ----map17, eval=FALSE---------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_surface (map, dat_B, dat=dat, cols=heat.colors (30))
#  cols_adj <- adjust_colours (heat.colors (30), -0.2)
#  map <- add_osm_surface (map, dat_HP, dat=dat, cols=cols_adj, size=1.5)
#  map <- add_osm_objects (map, dat_P, col=rgb (0.1,0.3,0.1))
#  map <- add_osm_objects (map, dat_H, col='gray60')
#  print (map)

## ----map17-print, echo=FALSE---------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_surface (map, dat_B, dat=dat, cols=heat.colors (30))
cols_adj <- adjust_colours (heat.colors (30), -0.2)
map <- add_osm_surface (map, dat_HP, dat=dat, cols=cols_adj, size=1.5)
map <- add_osm_objects (map, dat_P, col=rgb (0.1,0.3,0.1))
map <- add_osm_objects (map, dat_H, col='gray60')
png (height=mapht, width=mapwd, file='map_b17.png')
print (map)
graphics.off ()

## ----map18, eval=FALSE---------------------------------------------------
#  map <- add_colourbar (map, cols=terrain.colors (100), zlims=range (dat$z))
#  print (map)

## ----map18-print, echo=FALSE---------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_surface (map, dat_B, dat=dat, cols=heat.colors (30))
cols_adj <- adjust_colours (heat.colors (30), -0.2)
map <- add_osm_surface (map, dat_HP, dat=dat, cols=cols_adj, size=1.5)
map <- add_osm_objects (map, dat_P, col=rgb (0.1,0.3,0.1))
map <- add_osm_objects (map, dat_H, col='gray60')
map <- add_colourbar (map, cols=terrain.colors (100), zlims=range (dat$z))
png (height=mapht, width=mapwd, file='map_b18.png')
print (map)
graphics.off ()

## ----map19, eval=FALSE---------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_surface (map, dat_B, dat=dat, cols=heat.colors (30))
#  cols_adj <- adjust_colours (heat.colors (30), -0.2)
#  map <- add_osm_surface (map, dat_HP, dat=dat, cols=cols_adj, size=1.5)
#  map <- add_colourbar (map, cols=heat.colors (100), zlims=range (dat$z),
#                        alpha=0.9, vertical=FALSE,
#                        barwidth=c(0.1,0.12), barlength=c(0.5,0.9),
#                        text_col="blue", fontsize=5, fontface=3,
#                        fontfamily="Times")
#  map <- add_axes (map, colour="blue", fontsize=5, fontface=3,
#                   fontfamily="Times")
#  print (map)

## ----map19-print, echo=FALSE---------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_surface (map, dat_B, dat=dat, cols=heat.colors (30))
cols_adj <- adjust_colours (heat.colors (30), -0.2)
map <- add_osm_surface (map, dat_HP, dat=dat, cols=cols_adj, size=1.5)
map <- add_colourbar (map, cols=heat.colors (100), zlims=range (dat$z),
                      alpha=0.9, vertical=FALSE, 
                      barwidth=c(0.1,0.12), barlength=c(0.5,0.9),
                      text_col="blue", fontsize=5, fontface=3,
                      fontfamily="Times")
map <- add_axes (map, colour="blue", fontsize=5, fontface=3,
                 fontfamily="Times")
png (height=mapht, width=mapwd, file='map_b19.png')
print (map)
graphics.off ()

## ------------------------------------------------------------------------
d <- sqrt ((dat$x - mean (dat$x)) ^ 2 + (dat$y - mean (dat$y)) ^ 2)
range (d)

## ----map20, eval=FALSE---------------------------------------------------
#  dat <- dat [which (d < 0.01),]
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_surface (map, dat_B, dat=dat, cols=heat.colors (30),
#                          bg="gray40")
#  cols_adj <- adjust_colours (heat.colors (30), -0.2)
#  map <- add_osm_surface (map, dat_HP, dat=dat, cols=cols_adj, size=c (1.5, 0.5),
#                          bg="gray70")
#  print (map)

## ----map20-print, echo=FALSE---------------------------------------------
dat <- dat [which (d < 0.01),]
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_surface (map, dat_B, dat=dat, cols=heat.colors (30),
                        bg="gray40")
cols_adj <- adjust_colours (heat.colors (30), -0.2)
map <- add_osm_surface (map, dat_HP, dat=dat, cols=cols_adj, size=c (1.5, 0.5),
                        bg="gray70")
png (height=mapht, width=mapwd, file='map_b20.png')
print (map)
graphics.off ()

## ---- eval=FALSE---------------------------------------------------------
#  dat_T <- extract_osm_objects (key='tree', bbox=bbox)

## ---- echo=FALSE---------------------------------------------------------
dat_T <- london$dat_T

## ----map21, eval=FALSE---------------------------------------------------
#  map <- plot_osm_basemap (bbox=bbox, bg='gray20')
#  map <- add_osm_surface (map, dat_HP, dat=dat, cols=terrain.colors (30),
#                          size=c (1.5, 0.5), bg="gray70")
#  map <- add_osm_surface (map, dat_H, dat=dat, cols=terrain.colors (30),
#                          size=c (1, 0.5), bg="gray70")
#  map <- add_osm_surface (map, dat_T, dat=dat, cols=heat.colors (30),
#                                bg="lawngreen", size=c(3, 2), shape=c(8, 1))
#  print (map)

## ----map21-print, echo=FALSE---------------------------------------------
map <- plot_osm_basemap (bbox=bbox, bg='gray20')
map <- add_osm_surface (map, dat_HP, dat=dat, cols=terrain.colors (30), 
                        size=c (1.5, 0.5), bg="gray70")
map <- add_osm_surface (map, dat_H, dat=dat, cols=terrain.colors (30), 
                        size=c (1, 0.5), bg="gray70")
map <- add_osm_surface (map, dat_T, dat=dat, cols=heat.colors (30),
                              bg="lawngreen", size=c(3, 2), shape=c(8, 1))
png (height=mapht, width=mapwd, file='map_b21.png')
print (map)
graphics.off ()

