`latlon2idx` <-
function(lat,lon, layer){
nlon= length(layer$x)
xres=(max(layer$x)-min(layer$x))/nlon
ylon= length(layer$y)
yres=(max(layer$y)-min(layer$y))/ylon
#if(lon>=179.9999) lon=179.9999
x = floor((lon-min(layer$x))/xres)
y = floor((max(layer$y)-lat)/yres)
#idx = nlon*y+x+1
idx=list(x=(x+1),y=(ylon-y))
return(idx)
}

