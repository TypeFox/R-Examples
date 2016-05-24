#' Convert a spl object to GeoJSON format
#'
#' This function is used internally by \code{\link{writeMap}} to convert a given \code{spl} object to GeoJSON format.
#' @param x a \code{spl} object.
#' @param lightjson logical. Should GeoJSON code size be reduced by supressing extra whitespace characters and rounding numeric values?
#' 
#' @return A character string of GeoJSON formatted code.
#' @seealso \code{\link[rgdal]{writeOGR}} and its driver "\code{GeoJSON}" for GeoJSON export.
toGeoJSON <- function(x, lightjson=F){
  x.class <- class(x)
  x <- as.list(x)
  name.x <- safeVar(x["name"])
  coords.x <- x$coords
  if(lightjson){
    coords.x <- lapply(coords.x, function(x) lapply(x, round, digits=5))
  }
  holes.x <- x$holes
  order.x <- x$order
  x <- x[!names(x) %in% c("name", "coords", "holes", "order", "legend")]
  props <- paste("\"", names(x), "\"", sep="")
  lix <- list()
  
  for(i in 1:length(x[[1]])){
    lix[[i]] <- vector()
    for(j in 1:length(props)){
      lix[[i]][j] <- paste("\t\t\t", props[j], ": ", x[[j]][i], sep="")
    }
  }
  
  lix <- lapply(lix, paste, collapse=",\n")
  
  if(x.class == "splpoints"||x.class == "splicons"){
    coor <- coords2json(coords.x, x.class)
    res <- mapply(function(a, b) paste("\n\n{\"type\":\"Feature\",\n\"properties\": {\n", a, "},\n\"geometry\": {\n\"type\": \"Point\",\n\"coordinates\": ", b, "}}", collapse=",\n"),
                  a=lix, b=coor)
  }
  
  if(x.class == "spllines"){
    coor <- coords2json(coords.x, x.class)
    res <- mapply(function(a, b) paste("\n\n{\"type\":\"Feature\",\n\"properties\": {\n", a, "},\n\"geometry\": {\n\"type\": \"MultiLineString\",\n\"coordinates\": ", b, "}}", collapse=",\n"),
                  a=lix, b=coor)
  }
  
  if(x.class == "splpolygons"){
    coor <- coords2json(coords.x, x.class, holes.x, order.x)
    res <- mapply(function(a, b) paste("\n\n{\"type\":\"Feature\",\n\"properties\": {\n", a, "},\n\"geometry\": {\n\"type\": \"MultiPolygon\",\n\"coordinates\": ", b, "}}", collapse=",\n"),
                  a=lix, b=coor)
  }  
  
  res <- paste(res, collapse=",")
  if(lightjson){
    res <- gsub(" ", "", res)
    res <- gsub("\n", "", res)
    res <- gsub("\t", "", res)
  }
  res <- paste("var ", name.x, " = [", res, "\n]", sep="")
  return(res)
}
