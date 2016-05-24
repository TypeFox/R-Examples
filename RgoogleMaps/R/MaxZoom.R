`MaxZoom` <-structure(function#computes the maximum zoom level which will contain the given lat/lon range
### computes the maximum zoom level which will contain the given lat/lon range
(
  latrange, ##<< range of latitude values
  lonrange, ##<< range of longitude values
  size = c(640,640) ##<< desired size of the map tile image. defaults to maximum size returned by the Gogle server, which is 640x640 pixels
){
   SinPhi = sin(latrange * pi /180);
   normX = lonrange / 180;
   normY = (0.5 * log(abs((1 + SinPhi) / (1 -SinPhi) )) ) / pi;
   
   MaxZoom.lon <- floor(1 + log2(abs(size[1]/256/diff(normX))));
   MaxZoom.lat <- floor(1 + log2(abs(size[2]/256/diff(normY))));
   
   return(min(c(MaxZoom.lat=MaxZoom.lat,MaxZoom.lon=MaxZoom.lon)))
### zoom level
 })

