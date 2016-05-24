#' Export and display the map
#'
#' This function combines all the elements specified by the user and write
#' the corresponding HTML and Javascript code in a local directory.
#'
#' @param ... \code{basemap} and \code{spl} objects to embed into the map
#' @param dir a character string giving the directory path to export the map.
#' Default is the working directory.
#' @param prefix a character string to add a prefix to file names.
#' This allows multiple exportations in the same directory.
#' @param width,height the width and height of the map, in pixels.
#' @param setView a numeric vector of the form \code{c(x, y)}
#' setting the initial geographical center of the map.
#' @param setZoom a numeric value setting the initial map zoom.
#' @param interface an \code{ui} object created with \code{\link{ui}}
#' to customize the interface controls.
#' @param lightjson logical. Should GeoJSON files size be reducedby supressing
#' extra whitespace characters and rounding numeric values? Default is \code{FALSE}.
#' This is currently not compatible with RMarkdown popups.
#' @param directView a character string indicating if and how should the map be displayed.
#' Default option "\code{viewer}" uses (if available) the RStudio HTML viewer to display the map,
#' "\code{browser}" opens the map into the web browser and "\code{disabled}" disables direct display.
#' @param leaflet.loc a character string specifying the location (directory) of the leaflet library.
#' If set to "\code{online}" (default), the library is loaded from the leaflet
#' official CDN and requires an internet connection.
#' 
#' @export
writeMap <- function(..., dir=getwd(), prefix="", width=700, height=400,
                     setView=NULL, setZoom=NULL,
                     interface=NULL, lightjson=FALSE,
                     directView=c("viewer", "browser", "disabled"),
                     leaflet.loc="online"){
  ar <- list(...)
  depsub <- deparse(substitute(list(...)))
  depsub <- unlist(cleanDepsub(depsub))
  for(i in 1:length(ar)){
    if(is.null(ar[[i]]$name)){
      ar[[i]]$name <- depsub[i]
    }
    if(!is.null(ar[[i]]$legend)){
      ar[[i]]$legend$layer <- depsub[i]
      ar[[i]]$legend$layer.name <- ar[[i]]$name
      if(is.null(ar[[i]]$legend$title)){
        ar[[i]]$legend$title <- ar[[i]]$name
      }
    }
  }
  
  viewer <- getOption("viewer")
  user.view <- match.arg(directView)
  map.file <- paste(dir, "/", prefix,"_map.html", sep="")
  map.file.tmp <- paste(tempdir(), "/", prefix,"_map.html", sep="")
  
  writeMapInternal(ar=ar, dir=dir, prefix=prefix, width=width, height=height, setView=setView,
                   setZoom=setZoom, interface=interface,
                   lightjson=lightjson, leaflet.loc=leaflet.loc)
  
  if(user.view == "browser"){
    browseURL(map.file)
  }
  if(user.view == "viewer"){
    if(is.null(viewer)){
      warning("Cannot render in HTML viewer pane (require RStudio v.>0.98.5).
              Use 'directView' argument to render into your browser.")
    } else {

      writeMapInternal(ar=ar, dir=tempdir(), prefix=prefix, width=width,
                       height=height, setView=setView,
                       setZoom=setZoom, interface=interface,
                       lightjson=lightjson, leaflet.loc=leaflet.loc)
      if(leaflet.loc != "online"){
        #       file.copy(from=paste(leaflet.loc, "/leaflet.css", sep=""), to=tempdir(), overwrite=T)
        #       file.copy(from=paste(leaflet.loc, "/leaflet.js", sep=""), to=tempdir(), overwrite=T)
        #       file.copy(from=paste(leaflet.loc, "/images", sep=""), to=tempdir(), overwrite=T)
        warning("Cannot render in HTML viewer pane with a local copy of Leaflet.")
      }
      viewer(map.file.tmp)
    }
  }
}




writeMapInternal <- function(ar, dir, prefix, width, height, setView, setZoom,
                             interface, lightjson, leaflet.loc){
  
  ar.valid.class <- sapply(ar, function(x) is(x, "basemap") || is(x, "splpoints") || is(x, "splicons") || is(x, "spllines") || is(x, "splpolygons") || is(x, "splgrid"))
  if (any(ar.valid.class==FALSE)){
    stop("Invalid data format")
  }
  ar.names <- sapply(ar, function(x) return(safeVar(x$name)))
  if(any(duplicated(ar.names))){
    stop("Elements with duplicated names: ",
         paste(ar.names[duplicated(ar.names) | duplicated(ar.names, fromLast = TRUE)], collapse=", ")
    )
  }
    
  data.dir <- paste(dir, "/", prefix,"_data", sep="")
  if(!file.exists(data.dir))
    dir.create(data.dir)
  icons.dir <- paste0(data.dir, "/", prefix, "_icons")
  if(!file.exists(icons.dir))
    dir.create(icons.dir)
  icons.legend.dir <- paste0(icons.dir, "/", prefix, "_legend_icons")
  if(!file.exists(icons.legend.dir))
    dir.create(icons.legend.dir)
  
  #Include html+css code
  inc.encoding <- incEncoding(code = "UTF-8")
  inc.leaflet <- incLeaflet(loc = leaflet.loc)
  inc.data <- incData(prefix = prefix)
  init.map0 <- initMap0(height = height, width = width)
  init.map1 <- initMap1(setView, setZoom)
  inc.extra.css <- paste("<style type=\"text/css\">",
                         incPopupCSS(height = height, width = width),
                         incLegendCSS(),
                         incInfoPanelCSS(),
                         "</style>", sep = "\n\n")
  
  #Interface Controls
  ui.js <- uiJS(interface = interface, ar = ar)
  ui.js.1 <- ui.js$ui.1
  ui.js.2 <- ui.js$ui.2
  legend.js <- lapply(ar[sapply(ar, function(x) !is.null(x$legend))], function(x) processLegend(x$legend, icons.legend.dir = icons.legend.dir, prefix = prefix))
  legend.js <- paste0(legend.js, collapse = "\n\n")
  
  # Base Map
  bm <- ar[sapply(ar, function(x) is(x, "basemap"))]
  bm.js <- lapply(bm, toJS)
  bm.js <- do.call("paste", c(bm.js, sep="\n\n"))
  
  # Points  
  sppts <- ar[sapply(ar, function(x) is(x, "splpoints"))]
  sppts.json <- lapply(sppts, toGeoJSON, lightjson=lightjson)
  sppts.json <- do.call("paste", c(sppts.json, sep="\n\n\n\n"))
  write(sppts.json, paste(data.dir, "/", prefix,"_datapoints.js", sep=""))
  sppts.js <- lapply(sppts, toJS)
  sppts.js <- do.call("paste", c(sppts.js, sep="\n\n"))
  
  # Icons
  spico <- ar[sapply(ar, function(x) is(x, "splicons"))]
  icons.list <- lapply(spico, function(x) x$png)
  icons.list <- levels(as.factor(do.call("c", icons.list)))
  file.copy(from = icons.list, to = icons.dir)
  spico <- lapply(spico, function(x) {
    x$png <- paste("\"", prefix, "_data", "/", prefix, "_icons", "/",
                   gsub("(.*\\/)([^.]+\\.[[:alnum:]]+$)","\\2", x$png), "\"",
                   sep="")
    return(x)
  })
  spico.json <- lapply(spico, toGeoJSON, lightjson=lightjson)
  spico.json <- do.call("paste", c(spico.json, sep="\n\n\n\n"))
  write(spico.json, paste(data.dir, "/", prefix,"_dataicons.js", sep=""))
  spico.js <- lapply(spico, toJS)
  spico.js <- do.call("paste", c(spico.js, sep="\n\n"))
  
  # Lines  
  splns <- ar[sapply(ar, function(x) is(x, "spllines"))]
  splns.json <- lapply(splns, toGeoJSON, lightjson=lightjson)
  splns.json <- do.call("paste", c(splns.json, sep="\n\n\n\n"))
  write(splns.json, paste(data.dir, "/", prefix,"_datalines.js", sep=""))
  splns.js <- lapply(splns, toJS)
  splns.js <- do.call("paste", c(splns.js, sep="\n\n"))
  
  # Polygons  
  sppol <- ar[sapply(ar, function(x) is(x, "splpolygons"))]
  sppol.json <- lapply(sppol, toGeoJSON, lightjson=lightjson)
  sppol.json <- do.call("paste", c(sppol.json, sep="\n\n\n\n"))
  write(sppol.json, paste(data.dir, "/", prefix,"_datapolygons.js", sep=""))
  sppol.js <- lapply(sppol, toJS)
  sppol.js <- do.call("paste", c(sppol.js, sep="\n\n"))
  
  #Rasters
  spgrid <- ar[sapply(ar, function(x) is(x, "splgrid"))]
  raster.dir <- paste(data.dir, "/", prefix, "_rasters", sep="")
  if(!file.exists(raster.dir))
    dir.create(raster.dir)
  spgrid.js <- lapply(spgrid, function(x){
    url <- paste(raster.dir, "/", prefix, "_", safeVar(x$name), ".png", sep="")
    spgrid.js <- toJS(x, paste(prefix, "_data","/", prefix, "_rasters", "/", prefix, "_", safeVar(x$name), ".png", sep=""))
    png(url, bg = "transparent", width = 1200, height = round(1200 * pngasp(x$x.bbox)), type = "cairo")
    par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")    
    image(x$x, col = x$cells.col, asp = 1/cos((mean(x$x.bbox[1, ]) * pi)/180))
    dev.off()
    return(spgrid.js)
  })
  
  if(is.null(setView)){
    extbox <- .getExtBox(sppts, spico, splns, sppol, spgrid)
    init.map1 <- paste0(init.map1, "\n",
                       "map.fitBounds([[", extbox[3], ",", extbox[1],"], [", extbox[4], ",", extbox[2], "]]);")
  }
  
  # Compile
  map.file <- paste(dir, "/", prefix, "_map.html", sep="")
  write(paste(inc.encoding, inc.leaflet, inc.extra.css, inc.data, init.map0,
              "<script>", init.map1, ui.js.1, bm.js,
              sppol.js, spgrid.js, splns.js, sppts.js, spico.js,
              ui.js.2, legend.js, "</script>",
              sep="\n\n"), map.file)

}
