#' Renders a map for each of the species provided in names
#' 
#' This function map the species occurrences in the Iberian Peninsula.
#' 
#' For constructing these maps we used free online map repositories. 
#' We downloaded the Spanish administrative map from  http://www.diva-gis.org/ 
#' and the Spanish physical map of http://www.openstreetmap.org/ 
#' using the R- library OpenStreetMap.
#' 
#' @usage avisMapSpecies(names, maptype = "admin", ...)
#' @param names scientific name of the species 
#' (it could be a list of scientific names). E.g. "Pica pica"
#' @param maptype Available types of map are 'admin', 
#' administrative provinces of Spain (by default) 
#' or 'phys', physical map of Spain.
#' @param ... other filters passed to the observations query with avisQuerySpecies
#' @return a plot with the occurrences of the species in the Iberian Peninsula. Maps have high resolution, so they could be printed.
#' @export 
#' @examples \dontrun{
#' 
#' avisMapSpecies("Bubo bubo", "phys")
#' 
#'# if interested in several species, you can explore the database using avisMapSpecies
#'avisMapSpecies (list("Tyto alba", "Athene noctua", "Bubo bubo", "Strix aluco"), 
#'                maptype="phys")
#'                
#'# and you can save those maps individually using the tiff function
#'
#'directory<- "C:/your_directory"
#'species<- list("Tyto alba", "Athene noctua", "Bubo bubo", "Strix aluco")
#'for (x in species){
#'  tiff (file.path (directory, paste ("/", x, ".tiff", sep=""))) 
#'  avisMapSpecies (x)
#'  dev.off() 
#' }
#' 
#' }
#' 
avisMapSpecies<- function (names, maptype = 'admin', ...)
{
  if(!(is.list(names) | is.character (names))) {
    stop ("species names should be a list: e.g. avisMapSpecies (list ('Pica pica','Bucanetes githagineus')); or a character: avisMapSpecies ('Pica pica')")
  }
  
  for (name in names) 
  {
    obs<- avisQuerySpecies (name, args = list(...))
    avisMap(obs, name, maptype)
  }
}


#' Renders a map for the observations provided in 'obs'
#' 
#' This function should be used with avisQuerySpecies, to set a particular
#' query (with or without filters) and get the observations that we want to map. 
#' It just allow to map one species. See avisMapSpecies for multiple maps. 
#'  
#' @usage avisMap(obs, label = "", maptype = "admin")
#' @param obs set of observations returned by any of the avisQueryXXX functions
#' @param label label for the map. E.g. "Occurrences of Pica pica in Proyecto AVIS"
#' @param maptype Available types of map are 'admin', 
#' administrative provinces of Spain (by default) 
#' or 'phys, physical map of Spain.
#' @return a plot with the occurrences of the species in the Iberian Peninsula. Maps have high resolution, so they could be printed.
#' @export 
#' @examples \dontrun{
#' obs<- avisQuerySpecies ("Pica pica", args = list(habitat = "bosque"))
#' avisMap(obs, label = "Pica pica")
#' avisMap(obs, label = "Pica pica", maptype = "phys")
#'}
 
avisMap<-function(obs, label = '', maptype = 'admin')
{
  if(is.null(obs$x) || is.null(obs$y)){
    stop("missing x or y columns in provided obs parameter")
  }

  par(mar = c(0, 0, 0, 0))
  layout(matrix(c(1,1,1,1,1,1,1,1,2), 3, 3, byrow = TRUE))
  if (maptype=='phys'){
    .avisRenderMapPhysical(obs, label)
  } else if(maptype=='admin'){
    .avisRenderMapAdmin(obs, label)
  } else {
    stop(paste("Map type '", maptype, "' not available. Available map types are: 'admin' and 'phys'"));
  }
}

.avisRenderMapPhysical<-function(obs, label)
{
  peninsulaImg<-.avisReadPeninsulaImg()
  canariasImg<-.avisReadCanariasImg()

  plotRGB (peninsulaImg)
  points(obs$x, obs$y, col=alpha ("red", 0.5), pch=19, cex=1.2)
  text(-9.5, 34.2, label,  font=3, cex=2, adj=c(0,0))
  plotRGB (canariasImg)
  points(obs$x, obs$y, col=alpha ("red", 0.5), pch=19, cex=1.2)
}

.avisRenderMapAdmin<-function(obs, label)
{
  # ravis_shape_spain: shape in package data folder
  
  # hack for avoiding NOTE on check: 'no visible binding for global variable'
  # see: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
  ravis_shape_spain <- NULL
  rm(ravis_shape_spain)

  plot (ravis_shape_spain, border="grey75", ylim=c(34,44), xlim=c(-10,5))
  points(obs$x, obs$y, col=alpha ("red", 0.5), pch=19, cex=1.2)
  text(-9.5, 34.2, label,  font=3, cex=2, adj=c(0,0))
  plot (ravis_shape_spain, border="grey75", ylim=c(27.5, 29.5), xlim=c(-18.5,-13.5))
  rect(-18.5, 21, -11, 30, density = NULL, angle = 45,
   col = NA, border = "grey40", lwd=2)
  points(obs$x, obs$y, col=alpha ("red", 0.5), pch=19, cex=1.2)
}

.avisReadPeninsulaImg<-function()
{
  .avisCacheReturnOrSetup(".ravis_img_ipeninsula", function(){
  brick ( system.file('tif/peninsula.tif', package="rAvis"))
  })
}

.avisReadCanariasImg<-function()
{
  .avisCacheReturnOrSetup(".ravis_img_canarias", function(){

      brick ( system.file('tif/canarias.tif', package="rAvis"))

  })
}