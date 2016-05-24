setGeneric("readSolrDoc", 
           function(doc, ...)
             standardGeneric("readSolrDoc"))

setMethod("readSolrDoc", "character",
           function(doc, ...)
              readSolrDoc(xmlParse(doc), ...))

setMethod("readSolrDoc", "AsIs",
           function(doc, ...)
              readSolrDoc(xmlParse(doc), ...))

setMethod("readSolrDoc", "XMLInternalDocument",
           function(doc, ...)
              readSolrDoc(xmlRoot(doc), ...))


setMethod("readSolrDoc", "XMLInternalNode",
           function(doc, ...) {
             kids = xmlChildren(doc)
             kids = kids[!sapply(kids, inherits, "XMLInternalCommentNode")]
             if(length(kids) == 0)
                return(list())
             
             keys = sapply(kids, xmlGetAttr, "name")
             structure(lapply(kids, readSolrNodeValue),
                       names = keys)
           })


readSolrNodeValue =
function(node)
{
  id = xmlName(node)
  switch(id, int = if(abs(tmp <- as.numeric(xmlValue(node))) > .Machine$integer.max) tmp else as.integer(xmlValue(node)),
             long = as.numeric(xmlValue(node)),
             str = xmlValue(node),
             lst = readSolrDoc(node),
             bool = as.logical(xmlValue(node)),
             date = as.POSIXct(strptime(xmlValue(node), "%Y-%m-%dT%H:%M:%SZ")),
           )
}




  
