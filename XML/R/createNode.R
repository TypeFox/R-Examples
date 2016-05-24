xmlNode <-
function(name, ..., attrs = NULL, namespace = "",
          namespaceDefinitions = NULL, .children = list(...))
{
  kids <- lapply(.children, asXMLNode)
  kids = addNames(kids)

    # Check the names paces
  node <- list(name = name, attributes = attrs, children = kids, namespace=namespace,
                namespaceDefinitions = as(namespaceDefinitions, "XMLNamespaceDefinitions"))
  class(node) <- oldClass("XMLNode") # , "XMLAbstractNode")

  node
}

setAs("NULL", "XMLNamespaceDefinitions", function(from) structure(list(), class = "XMLNamespaceDefinitions"))

addNames =
function(kids, fromTag = TRUE)
{
  if(fromTag)
     names(kids) = sapply(kids, xmlName)
  else if(length(names(kids)) == 0)
      names(kids) <- sapply(kids,xmlName)
  else if(any( idx <- names(kids) == "")) 
      names(kids)[idx] <- sapply(kids[idx], xmlName)

  kids
}

setGeneric("xmlChildren<-",
function(x, ..., value) {
  standardGeneric("xmlChildren<-")
})

setMethod("xmlChildren<-", "ANY",
function(x, append = FALSE, ..., value) {
  #value = addNames(value)
  #x$children <- value
  addChildren(x, append = append, ..., kids = value)
})

# Don't need this.
# xmlName.character =
# function(node, full = FALSE)
#    "text"

setMethod("xmlChildren<-", "XMLInternalNode",
function(x, append = TRUE, ..., value) {
  if(!append)
     removeNodes(xmlChildren(x))
  if(!is.list(value))
    value = list(value)
  addChildren(x, kids = value)
})




addChildren =
function(node, ..., kids = list(...), at = NA, cdata = FALSE, append = TRUE)
  UseMethod("addChildren")


addChildren.XMLNode =  
function(node, ..., kids = list(...), at = NA, cdata = FALSE, append = TRUE)
{
  kids = lapply(kids,
                function(i) {
                  if(!inherits(i, "XMLNode"))
                    xmlTextNode(as.character(i), cdata = cdata)
                  else
                    i
                })

  node$children = if(append) c(node$children, kids) else kids
  node$children = addNames(node$children)
  
  node
}


# It would be better tokenize this, but ...
XMLEntities =
  c("&" = "amp",  # order is important as if we insert an entity, then the next time we will see the &.
    ">" = "gt",
    "<" = "lt",
    "'" = "apos",
    '"' = "quot")


insertEntities =
function(value, entities = XMLEntities)
{
    pat = names(entities)
    subs = paste("&", entities, ";", sep = "")
    for(i in seq(along = entities)) 
      value = gsub(pat[i], subs[i], value)

    value
}

xmlTextNode <- 
function(value, namespace = "", entities = XMLEntities, cdata = FALSE)
{
  node <- xmlNode("text", namespace = namespace)

  if(length(entities) && !inherits(value, "AsIs")) 
   value = insertEntities(value, XMLEntities)

  if(cdata)
    value = xmlCDataNode(value)
  
  node$value <- value
  if(!cdata)
     class(node) <- oldClass("XMLTextNode") # , class(node))
  if(length(entities))
    class(node) <- c(class(node), "EntitiesEscaped") #"XMLEntitiesEscapedTextNode"
  
  node
}


xmlPINode <-
function(sys, value, namespace="")
{
  x <- xmlNode(name=sys, namespace=namespace)
  x$value <- value
  class(x) <- oldClass("XMLProcessingInstruction") # , class(x))

 x
}

xmlCommentNode <-
function(text)
{
  node <- xmlTextNode(text)
  class(node) <- oldClass("XMLCommentNode") # , class(node))

  node
}

xmlCDataNode <-
function(...)
{
  txt <- paste(..., collapse="")  
 
  node <- xmlNode("text")
  node$value <- txt
  class(node) <- oldClass("XMLCDataNode") # , class(node))

 node
}

asXMLNode <-
function(x)
{
   #XXX
  if(!inherits(x, "XMLNode")) {
    xmlTextNode(x)
  } else {
    x
  }
}

