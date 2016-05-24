library('mapmisc')


coords = rbind(Alert = c(-62.338889, 82.501389),
    Qaanaaq = c(-69.238685,77.466335), 	
    'Alex Fjord' = c(-75.999722, 78.9),
    'Hans island' = c(-66.459722, 80.828056)
)

x = SpatialPointsDataFrame(
    coords, 
    data=data.frame(name=rownames(coords)),
    proj4string=crsLL
    )
        
if(require('rgdal', quietly=TRUE)) {

map = openmap(x, path='osm', verbose=TRUE, maxTiles=20, buffer=c(30,3))    

if(!interactive()) pdf(tempfile("osmplot", tmpdir=".", fileext=".pdf"))
map.new(map)
plot(map,add=TRUE)
points(x)
text(x, label=x$name, pos=4)
scaleBar(x, 'bottom')
scaleBar(x, 'left', seg.len=0, bty='n')
if(!interactive()) dev.off()


mapSat = openmap(x, path='mapquest', verbose=TRUE, maxTiles=20, buffer=c(30,3))    

if(!interactive()) pdf(tempfile("osmplot", tmpdir=".", fileext=".pdf"))
map.new(mapSat)
if(nlayers(mapSat) > 2) plotRGB(mapSat,add=TRUE)
points(x)
text(x, label=x$name, pos=4)
scaleBar(x, 'bottom')
scaleBar(x, 'left', seg.len=0, bty='n')
if(!interactive()) dev.off()

mapSat = openmap(x[x$name=='Hans island',], path='mapquest', 
    verbose=TRUE, buffer=c(4,0.01), zoom=6)    

if(!interactive()) pdf(tempfile("osmplot", tmpdir=".", fileext=".pdf"))
map.new(x[x$name=='Hans island',], buffer=0.3)
if(nlayers(mapSat) > 2) plotRGB(mapSat,add=TRUE)
points(x, pch=4, col='#FF000040', cex=5)
text(x, label=x$name, pos=4)
scaleBar(x, 'bottomright')
if(!interactive()) dev.off()


mapSat = openmap(x[x$name=='Qaanaaq',], path='mapquest', 
    verbose=TRUE, buffer=c(2,0.1), zoom=4)    

if(!interactive()) pdf(tempfile("osmplot", tmpdir=".", fileext=".pdf"))
map.new(x[x$name=='Qaanaaq',], buffer=0.5)
if(nlayers(mapSat) > 2) plotRGB(mapSat,add=TRUE)
points(x, pch=4, col='red', cex=5)
text(x, label=x$name, pos=4, col='red')
scaleBar(x, 'bottomleft')
if(!interactive()) dev.off()



xMerc = spTransform(x, omerc(x))
mapMerc = openmap(xMerc, path='osm', verbose=TRUE, 
    maxTiles=20, buffer=c(50,200)*1000)    

if(!interactive()) pdf(tempfile("osmplot", tmpdir=".", fileext=".pdf"))
map.new(xMerc, buffer=50000)
plot(mapMerc,add=TRUE)
points(xMerc)
text(xMerc, label=xMerc$name, pos=4)
scaleBar(xMerc, 'bottomleft')
scaleBar(xMerc, 'left', seg.len=0, bty='n')
if(!interactive()) dev.off()

}

if( FALSE) {
  # this takes too long
  xRot = spTransform(x, omerc(x, angle=-80))
  mapRot = openmap(xRot, path='osm', verbose=TRUE, 
    maxTiles=16, buffer=c(100,20,20,100)*1000)    

map.new(xRot, buffer=50000)
plot(mapRot,add=TRUE)
points(xRot)
text(xRot, label=xRot$name, pos=2+2*(xRot$name=='Qaanaaq'))
scaleBar(xRot, 'bottom')
scaleBar(xRot, 'left', seg.len=0, bty='n')
scaleBar(xRot, 'top', seg.len=0, bty='n')
scaleBar(xRot, 'topright', seg.len=0, bty='n')
scaleBar(xRot, 'bottomright', seg.len=0, bty='n')

}