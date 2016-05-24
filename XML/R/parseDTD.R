parseDTD <- 
function(extId, asText = FALSE, name = "", isURL = FALSE, error = xmlErrorCumulator())
{
  extId <- as.character(extId)
  if(!asText && missing(isURL)) {
    isURL <- length(grep("(http|ftp)://", extId, useBytes = TRUE))  > 0 
  }

  if(missing(name))
     name <- extId

  .oldErrorHandler = setXMLErrorHandler(error)
  on.exit(.Call("RS_XML_setStructuredErrorHandler", .oldErrorHandler, PACKAGE = "XML"), add = TRUE)

  if(asText) {
    f <- gsub("\\", "/", tempfile(), fixed=TRUE)    
    cat(extId, "\n", file = f)
    extId = f
    asText = FALSE
  }
  
 .Call("RS_XML_getDTD", as.character(name), as.character(extId),  
                          as.logical(asText), as.logical(isURL), error, PACKAGE = "XML")
}
