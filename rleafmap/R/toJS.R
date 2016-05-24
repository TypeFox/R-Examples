#' Generate Leaflet JS code for a given layer
#'
#' This function is used internally by \code{\link{writeMap}} to generate Leaflet JavaScript code for a given layer.
#' @param x a \code{spl} or \code{basemap} object.
#' @param url a character string giving the path for the raster files.
#' 
#' @return A character string of JavaScript Code.
toJS <- function(x, url=""){     # x is an spl or bm object

  if(is(x, "basemap")){
    if(!is.null(x$URL.cr)){
      credits <- paste("attribution: '", x$URL.cr, "',\n", sep="")
    } else {
      credits <- ""
    }
    res <- paste("var ", safeVar(x$name), "BaseMap = L.tileLayer('", x$URL, "', {\n",
                  "opacity: ", x$alpha, ",\n",
                  "minZoom: ", x$minZoom, ",\n",
                  "maxZoom: ", x$maxZoom, ",\n",
                  "tileSize: ", x$tileSize, ",\n",
                  "tms: ", x$tms,  ",\n",
                  credits,
                  "}).addTo(map);", sep="")
  }

  if(is(x, "splpoints")){
    res <- paste("var ", safeVar(x$name), "Points = L.geoJson(", safeVar(x$name), ", {
                 pointToLayer: function (feature, latlng) {
                 return L.circleMarker(latlng);
                 },
                 style: function (feature){
                 return {
                 radius: feature.properties.size,  
                 stroke: feature.properties.stroke,
                 color: feature.properties.strokeCol,
                 weight: feature.properties.strokeLwd,
                 dashArray: feature.properties.strokeLty,
                 opacity: feature.properties.strokeAlpha,
                 fill: feature.properties.fill,
                 fillColor: feature.properties.fillCol,
                 fillOpacity: feature.properties.fillAlpha
                 };
                 },
                 onEachFeature: function(feature, layer){
                 if (feature.properties.popup) {
                 layer.bindPopup(feature.properties.popup);
                 }
                 }
  }).addTo(map);", sep="")
  }

  if(is(x, "splicons")){
    res <- paste("var ", safeVar(x$name), "Icons = L.geoJson(", safeVar(x$name), ", {
                 pointToLayer: function (feature, latlng) {
                 return L.marker(latlng, {icon: L.icon({iconUrl: feature.properties.png, iconSize: feature.properties.size})});
                 },
                 onEachFeature: function(feature, layer){
                 if (feature.properties.popup) {
                 layer.bindPopup(feature.properties.popup);
                 }
                 }
  }).addTo(map);", sep="")
  }

  if(is(x, "spllines")){
    res <- paste("var ", safeVar(x$name), "Lines = L.geoJson(", safeVar(x$name), ", {
                 style: function (feature){
                 return {
                 stroke: feature.properties.stroke,
                 color: feature.properties.strokeCol,
                 weight: feature.properties.strokeLwd,
                 dashArray: feature.properties.strokeLty,
                 opacity: feature.properties.strokeAlpha,
                 };
                 },
                 onEachFeature: function(feature, layer){
                 if (feature.properties.popup) {
                 layer.bindPopup(feature.properties.popup);
                 }
                 }
  }).addTo(map);", sep="")
  }

  if(is(x, "splpolygons")){
    res <- paste("var ", safeVar(x$name), "Polygons = L.geoJson(", safeVar(x$name), ", {
                 style: function (feature){
                 return {
                 stroke: feature.properties.stroke,
                 color: feature.properties.strokeCol,
                 weight: feature.properties.strokeLwd,
                 dashArray: feature.properties.strokeLty,
                 opacity: feature.properties.strokeAlpha,
                 fill: feature.properties.fill,
                 fillColor: feature.properties.fillCol,
                 fillOpacity: feature.properties.fillAlpha
                 };
                 },
                 onEachFeature: function(feature, layer){
                 if (feature.properties.popup) {
                 layer.bindPopup(feature.properties.popup);
                 }
                 }
  }).addTo(map);", sep="")
  }
  
  if(is(x, "splgrid")){
    bb <- x$x.bbox
    res <- paste("var ", safeVar(x$name), "Raster = L.imageOverlay(\"", url, "\", [[",
                 bb[2,1], ",", bb[1,1], "],[",bb[2,2], ",",bb[1,2], "]]",
                 ").addTo(map);", sep="")
    }
  
  return(res)
}
