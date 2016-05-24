#' Provide national boundaries of a country.
#' 
#' 
#' This function uses the data set \code{\link{TMWorldBorders}}. One or several countries
#' can be selected and will be returned as an object of class 
#' \code{\link[sp:SpatialPolygons-class]{SpatialPolygons}}.
#' 
#' @param countries a vector of character string corresponding to the name of the countries 
#' you want to extract from the dataset. If \code{NULL}, a dialog box will be appear in order to select
#' the desired country.
#' @param multiple should the dialog box allow multiple selection (unused if \code{countries} is specified)?
#' @param proj map projection to use for the result (longitude and latitude in decimal degrees by default).
#' 
#' @details \code{proj} could be a character string corresponding to a \emph{PROJ.4} 
#' projection (see \url{http://trac.osgeo.org/proj/} for more details) 
#' or an object of class \code{\link[sp:CRS-class]{CRS}}\{\pkg{sp}\}.
#' 
#' @return Object of class \code{\link[sp:SpatialPolygons-class]{SpatialPolygons}}\{\pkg{sp}\}.
#' 
#' @note The result will be automatically plotted.
#' 
#' @seealso \code{\link{TMWorldBorders}}.
#' 
#' @examples 
#' \dontrun{
#' boundary <- create.boundary()
#' }
#' \dontshow{par(ask = TRUE)}
#' boundary <- create.boundary("Burkina Faso")
#' boundary <- create.boundary("Burkina Faso",
#'  proj="+proj=utm +zone=30 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
#' boundary <- create.boundary(countries = c("Burkina Faso", "Ghana", "Benin"))
#' \dontshow{par(ask = FALSE)}
#' 
#' @export
#' @keywords manip spatial

create.boundary = function(countries = NULL, multiple = F, proj = "+proj=longlat +ellps=WGS84"){
  ##################################################################################################
  # Cette fonction extrait de la base de donnees TM_WORLD_BORDERS-0.3 le contour des pays
  # Elle renvoie un object (boundary) de classe spatialPolygons contenant les frontieres des pays selectionnees
  # Les arguments sont
  # countries : Un vecteur character contenant le nom des pays que vous desirez extraire de la base de donnee
  #             Si countries est NULL vous pourrez selectionner  ala souris dans une liste  un ou plusieurs pays
  # multiple :  Un logical contenant T =>  on peut selectionner plusieurs pays
  #                                  F => on ne peut selectionner qu'un pays
  # proj :  La projection dans laquelle vous desirez boundary en sortie (character ou objet CRS)
  #
  #  Exemple d'utilisation
  # boundary      = create.boundary()
  # boundary      = create.boundary(proj="+proj=ortho +lat_0=0 +lon_0=0 +x0=0 +y0=0 +units=km") 
  # boundary      = create.boundary(countries = c("France","Italy","Belgium","Switzerland","Luxembourg","Germany","Austria","Slovakia","Czech Republic"))
  #
  ##################################################################################################
  data     = slot(prevR::TMWorldBorders,"data") 
  polygons = slot(prevR::TMWorldBorders,"polygons")
  if(is.null(countries)) {
      title = c(gettext("Select one country:",domain="R-prevR"),gettext("Select countries (use CTRL):",domain="R-prevR"))[c(!multiple,multiple)]
      countries = select.list(sort(as.character(data[["NAME"]])), multiple = multiple, title =  title)
  }
  if(length(countries)==0) return(NULL)
  ind       = match(countries,as.character(data[["NAME"]]),nomatch=0)
  if(any(ind==0)){
    invalid.names = paste(countries[ind==0],collapse=", ")
    n.invalid = length(countries[ind==0])
    sprintf(ngettext(n.invalid,"%s is not a valid country name.","%s are not valid country names.",domain="R-prevR"),invalid.names) -> stop.mess
    stop(stop.mess, call.=F)
  } 

  boundary = SpatialPolygons(polygons[ind],proj4string =CRS("+proj=longlat +ellps=WGS84"))
  if(!is.null(proj)){
    if(!is(proj,"CRS")){
      isOk = try(CRS(proj),silent=T)
      if(attr(isOk,"class") == "try-error"){
        stop(gettextf("the projection %s, defined in the 'proj' argument, is incorect.",proj,domain="R-prevR"), call.=F)
      }
      proj = CRS(proj)
    }
    boundary = spTransform(boundary,proj)
  }
  get("plot",pos="package:sp")(boundary,col = 1 + (1:length(boundary)),asp =1,pbg = "white",axes=T)
  if(length(countries)==1) title(main = countries)
    message("Source: World Borders Dataset 0.3 (2008)\nProvided by Bjorn Sandvik, http://thematicmapping.org/downloads/world_borders.php\nThe dataset was derived by Schuyler Erle from public domain sources.\nSean Gilles did some clean up and made some enhancements.\nThe dataset is available under a Creative Commons Attribution-Share Alike License.\nhttp://creativecommons.org/licenses/by-sa/3.0/\nThe boundaries, names designations used do not imply official endorsement or acceptance by the authors.",domain="R-prevR")
  boundary
}

