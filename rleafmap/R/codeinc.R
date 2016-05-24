
incLeaflet <- function(loc){
  if(loc == "online"){
    res <- paste("<link rel=\"stylesheet\" href=\"http://cdn.leafletjs.com/leaflet-0.7.2/leaflet.css\" />
<script src=\"http://cdn.leafletjs.com/leaflet-0.7.2/leaflet.js\"></script>")
  } else {
    res  <- paste("<link rel=\"stylesheet\" href=\"", loc, "/leaflet.css\" />
<script src=\"", loc, "/leaflet.js\"></script>", sep="")
  }
return(res)
}

incEncoding <- function(code){
  res <- paste("<?xml version=\"1.0\" encoding=\"utf-8\"?>
<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\" 
\"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\">
<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"fr\">
<meta http-equiv=\"Content-Type\" content=\"application/xhtml+xml; charset=utf-8\" />", sep="")
  return(res)
}

incPopupCSS <- function(height, width){
  res <- paste0(
".leaflet-popup-content {
  padding-right:20px !important;
  width:auto !important;
  max-width:", width*0.75, "px !important;
  max-height:", height*0.75, "px !important;
  overflow:auto !important;
}")
  return(res)
}

incData <- function(prefix){
  paste("<script src=\"", prefix,"_data/", prefix,"_datapoints.js\"></script>
<script src=\"", prefix,"_data/", prefix,"_dataicons.js\"></script>
<script src=\"", prefix,"_data/", prefix,"_datalines.js\"></script>
<script src=\"", prefix,"_data/", prefix,"_datapolygons.js\"></script>",  # Fix le 2eme prefix + faire une fonction
  sep="")
}

initMap0 <- function(height, width){
  paste("<div id=\"map\" class=\"map\" style=\"height: ", height,"px; width: ",width, "px\"></div>", sep="")
}

initMap1 <- function(setView, setZoom){
  if(is.null(setView)) setView <- c(0, 0)
  if(is.null(setZoom)) setZoom <- 1
  paste("var map = L.map('map', {zoomControl:false, attributionControl:false}).setView([", setView[1], ", ", setView[2],"], ", setZoom, ");", sep="")
}

incInfoPanelCSS <- function(){
  ".info {
    padding: 6px 8px;
    background: rgba(255,255,255,0.8);
    box-shadow: 0 0 15px rgba(0,0,0,0.2);
    border-radius: 5px;
  }"
}

incLegendCSS <- function(){
  ".legend i {
    clear: left;
    float: left;
    margin-right: 8px;
    margin-top: 5px;
  }
  .legend p {
    float: left;;
    line-height: 5px;
  }
  .legend h1 {
    font-size: 10px;
    margin-top: 0px;
    margin-bottom: 2px;
  }
  .legend hr {
    margin-top: 2px;
    margin-bottom: 2px;
    display: block;
	  height: 1px;
    border: 0;
    border-top: 1px solid #AAAAAA;
  }"
}