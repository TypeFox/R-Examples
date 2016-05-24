`LatLon2XY` <- structure(function#computes the coordinate transformation from lat/lon to map tile coordinates 
### The function LatLon2XY(lat,lon,zoom) computes the coordinate transformation from lat/lon to map tile coordinates given a zoom level. 
### It returns the tile coordinates as well as the pixel coordinates within the Tile itself. 
### Thanks to Neil Young (see \url{http://groups.google.com/group/Google-Maps-API/browse_thread/thread/d2103ac29e95696f?hl=en}) for providing the formulae used.
(
  lat, ##<< latitude values to transform
  lon, ##<< longitude values to transform
  zoom ##<< zoom level.lat,lon,zoom
){
  #The int part may be used directly to obtain the Tiles (Attention: Illegal) using a URL formed like this:

  #String.Format("http://mt.google.com/mt?x={0}&y={1}&zoom={2}", c.x, c.y, 17-zoom)

  ##note<< The fractional part times 256 is the pixel coordinate within the Tile itself. 
  
  #double sin_phi = System.Math.Sin(this.lat * System.Math.PI /180);
  #double norm_x = this.lon / 180;
  #double norm_y = (0.5 * System.Math.Log((1 + sin_phi) / (1 -sin_phi))) / System.Math.PI;
  #this.y = System.Math.Pow(2, this.zoom) * ((1 - norm_y) / 2);
  #this.x = System.Math.Pow(2, this.zoom) * ((norm_x + 1) / 2); 
  SinPhi = sin(lat * pi /180);
  normX = lon / 180;
  normY = (0.5 * log((1 + SinPhi) / (1 -SinPhi))) / pi;
  Y = (2^zoom) * ((1 - normY) / 2);
  X = (2^zoom) * ((normX + 1) / 2); 	

  x = 256 *(X- floor(X));
  y = 256 *(Y- floor(Y));
  
  return(list(Tile = cbind(X=floor(X),Y=floor(Y)), Coords = cbind(x=x,y=y)))
###  A list with values
###  \item{Tile}{integer numbers specifying the tile}
###  \item{Coords}{pixel coordinate within the Tile}
}, ex = function(){
  LatLon2XY(38.45, -122.375, 11)
})


