#' Open a postscript device with the color scheme given by geoplotpalette.
#' 
#' The function starts a postscript device and the arguments are the same as
#' the arguments to postscript (height, width, file, bg etc).
#' 
#' <!--explain details here-->
#' 
#' @param ... The arguments to postscript are optional. See help postscript.
#' @return No value returned.
#' @section Side Effects: A graphics device is opened (often a file).  It must
#' be closed again by dev.off()
#' @seealso postscript, geoplotpalette,litir
#' @keywords <!--Put one or more s-keyword tags here-->
#' @examples
#' 
#' \dontrun{
#' colps(file="map1.ps",height=6,width=5) 
#' geoplot(xlim=c(-28,-10),ylim=c(64,69))
#' geosymbols(data,z=data$value,circles=0.2,sqrt=T)
#' geopolygon(island,col="white")# paint white over the symbols 
#' geolines(island) # that are inside the country.  (island) 
#' dev.off() 
#' 
#' # same example in a different way.  
#' colps(file="map1.ps",height=6,width=5,bg="white") 
#' geoplot(xlim=c(-28,-10),ylim=c(64,69))
#' geosymbols(data,z=data$value,circles=0.2,sqrt=T)
#' geopolygon(island,col=0)#col 0 is now white
#' geolines(island) # was transparent earlier
#' dev.off() 
#' }
#' @export colps
colps <-
function(...){
  postscript(...) 
  geoplotpalette()
}

