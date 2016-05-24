xmlInitializeCatalog = 
function()
   .C("R_xmlInitializeCatalog", PACKAGE = "XML")


catalogResolve =
function(id, type = "uri", asIs = FALSE, debug = FALSE)
{
   xmlInitializeCatalog()
   type = rep(type, length = length(id))
  
   types = c("uri", "public", "system")
   i = pmatch(tolower(type), types, duplicates.ok = TRUE)
   if(any(is.na(i)))
     stop("don't recognize type. Must be one of ", paste(types, collapse = ", "))

   ans = .Call("R_xmlCatalogResolve", as.character(id), i, as.logical(debug), PACKAGE = "XML")

   if(asIs)
     ans[is.na(ans)] = id[is.na(ans)]

   ans
}  


catalogLoad =
function(fileNames)
{
  .Call("RS_XML_loadCatalog", path.expand(fileNames), PACKAGE = "XML")
}  

catalogClearTable =
function()
{
  .Call("RS_XML_clearCatalog", PACKAGE = "XML")
}


XMLCatalogTypes = c("public", "system", "rewriteSystem", "rewriteURI", "uri", "delegateSystem", "delegatePublic", "delegateURI", "nextCatalog", "catalog")


catalogAdd =
function(orig, replace, type = "rewriteURI")
{
  if(missing(replace)) {
     replace = orig
     orig = names(replace)
  }
  else
     length(replace) = length(orig)
  
  idx = pmatch(type, XMLCatalogTypes)
  if(any(is.na(idx))) {
    stop("unrecognized XML catalog type(s) ", type[is.na(idx)], ". Must be one of ",
           paste("'", XMLCatalogTypes, "'", sep = "", collapse = ", "))
  }
  type = XMLCatalogTypes[idx]
  type = rep(as.character(type), length = length(orig))

  xmlInitializeCatalog()
  .Call("RS_XML_catalogAdd", as.character(orig), as.character(replace), as.character(type), PACKAGE = "XML")
}  


catalogDump =
  #
  # Get a snapshot of the current contents of the global catalog table
  # parsing it or writing it to a file for further use 
  # If asText = TRUE and you don't specify a value for fileName,
  # it returns the XML content as a string for easier viewing.
  
function(fileName = tempfile(), asText = TRUE)
{
  xmlInitializeCatalog()
  ans = .Call("RS_XML_catalogDump", as.character(fileName), PACKAGE = "XML")
  if(missing(fileName)) {
    ans = xmlParse(fileName)
    if(asText)
      ans = saveXML(ans)
    unlink(fileName)
  }

  ans
}  

