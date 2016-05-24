##
## Function: CreateXmlString
##
CreateXmlString <- function(neosxml, cdatalist){
  if(!(class(neosxml) == "NeosXml")){
    stop("\nPlease provide an object of class 'NeosXml' for argument neosxml.\n")
  }
  xmlstr <- neosxml@xml
  lnames <- names(cdatalist)
  if(!(all(lnames %in% names(xmlstr)))){
    stop("\nNamed list object 'cdatalist' does contain entries that are not node names of object 'xmlstr'.\n")
  }
  idx <- 1:length(cdatalist)
  tmp <- removeChildren(xmlstr, kids = lnames)
  for(i in idx){
    tmp <- append.XMLNode(tmp, newXMLNode(lnames[i], cdatalist[[i]], cdata = TRUE))
  }
  ans <- toString(tmp)
  return(ans)
}             
