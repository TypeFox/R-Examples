# Purpose        : Open and close KML file;
# Maintainer     : Pierre Roudier (pierre.roudier@landcare.nz);
# Contributions  : Dylan Beaudette (debeaudette@ucdavis.edu); Tomislav Hengl (tom.hengl@wur.nl);   
# Status         : tested (but not on the MacOS / Linux)
# Note           : See [http://code.google.com/apis/kml/documentation/kmlreference.html] for more info.

kml_open <- function(
  file.name,
  folder.name = file.name,
  kml_open = TRUE,
  kml_visibility = TRUE,
  overwrite = TRUE,
  use.Google_gx = FALSE,
  kml_xsd = get("kml_xsd", envir = plotKML.opts),
  xmlns = get("kml_url", envir = plotKML.opts),
  xmlns_gx = get("kml_gx", envir = plotKML.opts)
  ){

  if (file.exists(file.name) & overwrite==FALSE) {
    stop(paste("File", file.name, "already exists. Set the overwrite option to TRUE or choose a different name."))
  }

  # header
  if(use.Google_gx){
    kml.out <- newXMLNode("kml", attrs=c(version="1.0"), namespaceDefinitions = c("xsd"=kml_xsd, "xmlns"=xmlns, "xmlns:gx"=xmlns_gx))
  }
  else {
    kml.out <- newXMLNode("kml", attrs=c(version="1.0"), namespaceDefinitions = c("xsd"=kml_xsd, "xmlns"=xmlns))
  }
  
  h2 <- newXMLNode("Document", parent = kml.out)
  h3 <- newXMLNode("name", folder.name, parent = h2)
  if(kml_visibility==FALSE){
    h4 <- newXMLNode("visibility", as.numeric(kml_visibility), parent = h2)
  }
  h5 <- newXMLNode("open", as.numeric(kml_open), parent = h2)
  
  # init connection to an XML object: 
  assign("kml.out", kml.out, envir=plotKML.fileIO)
  message("KML file opened for writing...")
  
}

## Closes the current KML canvas
kml_close <- function(file.name, overwrite = FALSE) {
  
  # get our invisible file connection from custom evnrionment
  kml.out <- get("kml.out", envir=plotKML.fileIO)
  saveXML(kml.out, file.name)
  message(paste("Closing ", file.name))
  
}

## Open the KML file using the default OS application:
kml_View <- function(file.name){
  if(.Platform$OS.type == "windows") {
      ext <- raster::extension(file.name)
      x <- NULL # set default value for error checking
      if(!inherits(try({ x <- utils::readRegistry(ext, hive="HCR") }, silent = TRUE), "try-error")){
        if(! is.null(x[which(names(x) %in% c('Content Type', '(Default)'))])){
          system(paste("open ", utils::shortPathName(normalizePath(paste(getwd(), "/", file.name, sep=""))), sep=""))
        }
      }
      else{
        warning(paste("No MIME type detected for", ext, "file extension."))
      }
  } 

  ## DB: this is untested on Mac OS X / UNIX!
  else {
      if(length(grep("mac.binary", .Platform$pkgType))>0){
        if(inherits( try(system(paste("open ", normalizePath(paste(getwd(), "/", file.name, sep="")), sep="")), silent = TRUE), "try-error")){
          warning(paste("'open' resulted in a try-error. Consider re-installing the default", ext, "browser."))
         } 
        }
      else{
        if(inherits( try(system(paste("gnome-open ", normalizePath(paste(getwd(), "/", file.name, sep="")), sep="")), silent = TRUE), "try-error")){
          warning(paste("'gnome-open' resulted in a try-error. Consider re-installing the default", ext, "browser."))
         } 
      }
  }
}

# end of script;
