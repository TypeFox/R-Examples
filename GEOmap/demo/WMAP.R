
data(worldmap)
data(ETOPO5)

data(USAmap)
data(usacity)
data(worldcity)


######### plot world map:
plotworldmap(worldmap, add = FALSE, col = 2)

#######  set longitude rotation to 180 to get Greenwich in center of map
shiftlon = 180

readline("Hit Return>")
######  plot world map:
 plotworldmap(worldmap, shiftlon = shiftlon, add = FALSE, col = 2)

readline("Hit Return>")
###########  get USA in center of Map
shiftlon = 60
plotworldmap(worldmap, shiftlon = shiftlon, add = FALSE, col = 2)


#######  choose a "target region" to focus
A = locworld(shiftlon)

#####  for New York area this is:
A = list(lon=c(281.9203, 291.9571), lat=c(35.58868, 44.55468),
  LON=c(281.9203, 291.9571), LAT=c(35.58868, 44.55468),
  utmbox=list(x=17, y="S", lon=-84, lat=32),
  UTM0 =list(lam=-81 , phi=36 ),
  shiftlon=60)
  

PROJ = setPROJ(type=2, LAT0=mean(A$LAT) , LON0=mean(A$LON) )


NYLIM = c(A$LON[1], A$LAT[1],A$LON[2], A$LAT[2] )


#######  very crude map:
readline("Hit Return>")
plotGEOmapXY(worldmap, PROJ=PROJ, LIM=NYLIM ,  add=FALSE)

#######  better USA map
readline("Hit Return>")
plotGEOmapXY(USAmap, PROJ=PROJ, LIM=NYLIM ,  add=FALSE)

##########  set up map polygons for continents
 LMAP = SETPOLIMAP()

P = list(lat=mean(A$LAT), lon=mean(A$LON))

J = LOCPOLIMAP(P, LMAP)
                
pmap2 = selectPOLImap(which(J==1), 2)

####pmap2$STROKES$code

pmap2$STROKES$col[pmap2$STROKES$code=="2"]='green'
pmap2$STROKES$col[pmap2$STROKES$code=="3"]='purple'
pmap2$STROKES$col[pmap2$STROKES$code=="4"]='red'

####  zz = map('state', region = c('new york', 'new jersey', 'penn'))
###  how to extract info from the
###    mapdata database


library(maps)

zz = map('county', 'new york', plot = FALSE)

cntyNY.xy = GLOB.XY(zz$y, zz$x, PROJ)
lines(cntyNY.xy$x,  cntyNY.xy$y, col='red' )
plotGEOmapXY(pmap2, PROJ=PROJ  , LIM=NYLIM, add=TRUE)

###  the lakes and Islands are not defined in the cia data bases, not good.

readline("Hit Return>")

#########  run the interactive map program to see topography
DOTOPOMAPI(TOPO=ETOPO5, themap=worldmap, shiftlon=180)

readline("Hit Return>")


DOTOPOMAPI(TOPO=ETOPO5, themap=worldmap, shiftlon=180,  DOCONT=TRUE, DOIMG=TRUE,
                    PNTS=NULL, PCOL='red', PCH=1, PCEX=1, PS=FALSE, MAPPATH=NULL,
                     polybase = NULL,
                     USAmap= USAmap,
                     usacity=usacity,
                     worldcity=worldcity)


######### not run: #########
#########
#########  library(ETOPO); data(ETOPO2)
#########  DOTOPOMAPI(TOPO=ETOPO2, themap=worldmap, shiftlon=180)
#########

