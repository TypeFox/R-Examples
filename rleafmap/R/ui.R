#' Options settings for map interface
#'
#' Allow the user to choose which interface elements are displayed on the map and their positions.
#' @param zoom a character string indicating if and how should the zoom control be displayed.
#' This must be one of "\code{topleft}", "\code{topright}", "\code{bottomleft}", "\code{bottomright}", "\code{none}"
#' @param layers a character string indicating if and how should the layers control be displayed.
#' This must be one of "\code{topleft}", "\code{topright}", "\code{bottomleft}", "\code{bottomright}", "\code{none}"
#' @param attrib a character string indicating if and how should the attribution control be displayed.
#' This must be one of "\code{topleft}", "\code{topright}", "\code{bottomleft}", "\code{bottomright}", "\code{none}"
#' @param attrib.text a character string for additionnal credits. HTML tags are accepted.
#' 
#' @export
#' @return An object of class \code{ui} which can be directly given as \code{interface} argument of \code{\link{writeMap}}.
ui <- function(zoom=c("topleft", "topright", "bottomleft", "bottomright", "none"),
               layers=c("none", "topright", "topleft", "bottomleft", "bottomright"),
               attrib=c("bottomright", "topleft", "topright", "bottomleft", "none"),
               attrib.text=""){
  res <- list(zoom=match.arg(zoom),
              layers=match.arg(layers),
              attrib=match.arg(attrib),
              attrib.text=attrib.text)
  class(res) <- "ui"
  return(res)
}

#' Generate user interface JS code
#'
#' This function is used internally by \code{\link{writeMap}} to generate JavaScript code
#' related to the user interface.
#' @param interface an \code{ui} object created with \code{\link{ui}}.
#' @param ar a list of \code{basemap} and \code{spl} objects.
uiJS <- function(interface, ar){
  if(is.null(interface)){
    interface <- ui()
  }
  if(!is(interface, "ui")){
    stop("interface must be of class 'ui'")
  }
  if(interface$zoom != "none"){
    zoomInterface <- paste("L.control.zoom({position:'", interface$zoom, "'}).addTo(map);", sep="")
  } else {
    zoomInterface <- ""
  }
  if(interface$attrib != "none"){
    attribInterface <- paste("var attrib = L.control.attribution({prefix:false, position:'", interface$attrib, "'}).addTo(map);\n
                               attrib.addAttribution(\"", interface$attrib.text,"\");", sep="")
  } else {
    attribInterface <- ""
  }
  
  if(interface$layers != "none"){
    arNames <- xvarnames(ar)
    arNames.bl <- arNames[arNames$xclass == "basemap",]
    arNames.ol <- arNames[arNames$xclass != "basemap",]
    if(dim(arNames.bl)[1] == 0){
      bl.js <- "var baseMaps = 1"
    }else{
      bl.list <- paste("\"", arNames.bl$xname, "\" : ", arNames.bl$xvarname, sep="", collapse=",\n")
      bl.js <- paste("var baseMaps = {\n", bl.list, "\n};", sep="")
    }
    if(dim(arNames.ol)[1] == 0){
      ol.js <- "var overlayMaps = 1"
    }else{
      ol.list <- paste("\"", arNames.ol$xname, "\" : ", arNames.ol$xvarname, sep="", collapse=",\n")
      ol.js <- paste("var overlayMaps = {\n", ol.list, "\n};", sep="")
    }
    ctrl <- paste("L.control.layers(baseMaps, overlayMaps, {position:'", interface$layers, "'}).addTo(map);", sep="")
    layersInterface <- paste(bl.js, ol.js, ctrl, sep="\n")
  } else {
    layersInterface <- ""
  }
  
  res <- list(ui.1 = paste(zoomInterface, attribInterface, sep = "\n\n"), ui.2 = layersInterface)
  return(res)
}