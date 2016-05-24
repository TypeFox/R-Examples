xmlTree <-
  #
  # Create an XML document using internal nodes and help to manage
  # the state for the user rather than requiring them to manage
  # the individual nodes. For the most part, the two approaches
  # are relatively similar in complexity.
  #
  #
function(tag = NULL, attrs = NULL, dtd = NULL, namespaces = list(),
          doc = newXMLDoc(dtd, namespaces))
  # Allows a DOCTYPE, etc. at the beginning by specifying dtd as 
  # a vector of 1, 2, 3 elements passed to newXMLDTDNode() or
  # as an XMLDTDNode directly.
  
  #
{
 currentNodes <- list(doc)  # the stack of nodes
 
 isXML2 <- libxmlVersion()$major != "1" 

     # if we are given a DTD, add it to the document.
 if(!is.null(dtd)) {
   if(isXML2) {
     node = NULL
     if(inherits(dtd, "XMLDTDNode"))
       node = dtd
     else if(is.character(dtd) && dtd[1] != "")
       node = newXMLDTDNode(dtd, doc = doc)

     if(!is.null(node)) {
       addChildren(doc, node)
       currentNodes[[2]] <- node #???XXX
     }
   } else
     warning("DTDs not supported in R for libxml 1.*. Use libxml2 instead.")
 }
 
 definedNamespaces = list()
 defaultNamespace = NULL
 addNamespaceDefinitions = is.null(tag)
 
 setActiveNamespace = function(ns) {
                         defaultNamespace <<- ns
                      }
 
 asXMLNode <- function(x) {
        if(inherits(x, "XMLInternalNode"))
          return(x)
        
        v = if(is.list(x)) 
               lapply(x, asXMLNode)
            else 
               newXMLTextNode(as.character(x), doc = doc, escapeEntities = is(x, "AsIs"))

        v 
      }


 
 setNamespace <- function(node, namespace = defaultNamespace) {
  
         # if there is no namespace or if we have one and no names on the namespace
      if(length(namespace) == 0 || ! ( length(namespace) == 1 && is.null(names(namespace)) ) )
       return(NULL)

     if(is.list(namespace))
       return(NULL)

      
     if(!is.na(match(namespace, names(namespaces))) && is.na(match(namespace, names(definedNamespaces)))) {
       ns <- .Call("R_xmlNewNs", node, namespaces[[namespace]], namespace, PACKAGE = "XML")
       definedNamespaces[[namespace]] <<- ns
     }

     setXMLNamespace(node,  definedNamespaces[[namespace]])
#old     setInternalNamespace( node, definedNamespaces[[namespace]])
 }

 # namespace is intended to be the namespace for this node
 # and not any definitions.
 # How do we define new namespaces with this function?
 # Can we add them to attrs. No!
 addTag <- function(name, ..., attrs = NULL,
                    close = TRUE, namespace = defaultNamespace, .children = list(...) )
 {
   if(inherits(name, "XMLInternalNode")) {
     addChildren(currentNodes[[1]], name)
     currentNodes <<- c(node, currentNodes)
     addChildren(node, kids = .children)
     if(close)
        currentNodes <<- currentNodes[-1]
     return(name)
   }
   
   # if the user gives us something like "r" for the namespace as opposed to
   #  c(r = "http:...") then we try to match the prefix in an earlier node
   # ??? Should we use the defined namespaces in the document?
if(FALSE) {   
   if(length(namespace) == 1 && length(names(namespace)) == 0) {
     tmp = namespace
     if(length(currentNodes)) {
       defs = namespaceDeclarations(currentNodes[[1]], TRUE)
       i = match(namespace, names(defs))
       if(!is.na(i))
         namespace = defs[[i]]
     } 
   }
 }

   
   if(!is.null(attrs))
      storage.mode(attrs) <- "character"

   if(inherits(name, "XMLInternalNode"))
      node = name
   else {
      parent = if(length(currentNodes) > 1) 
                    currentNodes[[1]] 
               else 
                    xmlRoot(currentNodes[[1]])
      node <- newXMLNode(name, attrs = attrs, namespace = namespace, 
                         doc = doc,  parent = parent,
                          namespaceDefinitions = if(addNamespaceDefinitions) namespaces else NULL)

      if(addNamespaceDefinitions) {
#       lapply(seq(along = namespaces),
#               function(i)
#                   setXMLNamespace(node, namespaces[[i]], names(namespaces)[i]))
        addNamespaceDefinitions <<- FALSE
      }
   }

#   if(length(currentNodes) > 1) 
#      addChildren(currentNodes[[1]], node)

   currentNodes <<- c(node, currentNodes)

#   if(!inherits(name, "XMLInternalNode"))
#      setNamespace(node, namespace)      

   for(i in .children) 
      addChildren(node, asXMLNode(i))  # vectorize XXX

   if(close == TRUE)
     closeTag()
   
   invisible(node)
 }


 closeTag <- function(name="") {

   if(nargs() == 0) {
     tmp <- currentNodes[[1]]
     currentNodes <<- currentNodes[-1]
   } else if( is.character(name) ) {

     w = sapply(currentNodes, inherits, "XMLInternalElementNode")
     useNamespace = length(grep(":", name)) > 0
     ids = sapply(currentNodes[ w ], xmlName, useNamespace)
     tmp = list()
     for(id in name) {
        i = which(id == ids)
        if(length(i) == 0)
          stop("Cannot close tag for node with name ", id, " - no such node open")
        tmp = c(tmp, currentNodes[1:i])
        currentNodes <<- currentNodes[-c(1:i)]
        ids = ids[-(1:i)]
     }
     
   } else if(inherits(name, "numeric")) {
       num = name
       if(is.na(num) || num == -1) 
              # close all of the nodes, except the document node.
           w = seq(along = currentNodes[- length(currentNodes)])
       else if(length(num) == 1) 
           w = 1:num
       else
           w = num
       tmp = currentNodes[ w ]
       currentNodes <<- currentNodes[ - w ]
   }


  invisible(tmp)
 }


 add = function(node, parent = currentNodes[[1]], close = TRUE) {
        if(!is.null(parent)) {
            addChildren(parent, node)
            if(!close)
              currentNodes <<- c(node, currentNodes)
        }
        invisible(node)
       }
 
 addComment <- function(...) {
   add(newXMLCommentNode(paste(as.character(list(...)), sep=""), doc = doc))
 }


 addCData <- function(text) {
   add(newXMLCDataNode(text, doc = doc))
 }

 addPI <- function(name, text) {
   add(newXMLPINode(name, text, doc = doc), NULL)
 }


   # deal with the top-level node the user may have supplied.
 if(!is.null(tag)) {
   if(is.character(tag)) {
     node = addTag(tag, attrs = attrs, namespace = namespaces, close = FALSE)
   } else if(inherits(tag, "XMLInternalNode")) {
     if(is.null(xmlParent(node))) # if we have a DTD node, need to add it to that or parallel to that?
       addChildren(doc, node)
   }
   

 }

 v <- list(
         addTag = addTag,
         addNode = addTag,           
         addCData = addCData,
         addPI = addPI,
         closeTag = closeTag,
         closeNode = closeTag,
         addComment = addComment,
         setNamespace = setActiveNamespace,
         value = function() doc,
         doc = function() doc,
         add = function(...){}
       )

 #class(v) <- c("XMLInternalDOM", "XMLOutputStream")
 # v
 ans = new("XMLInternalDOM", v)
 names(ans) = names(v)
 ans
}

setAs("XMLInternalNode", "XMLNode",
        function(from) 
           asRXMLNode(from)
        )


xmlRoot.XMLInternalDOM =
function(x, skip = TRUE, ...)
{
  xmlRoot(x$doc(), skip = skip)
}  


 #??? This was XMLInternalElement and not ...Node
xmlRoot.XMLInternalElement = xmlRoot.XMLInternalNode =
function(x, skip = TRUE, ...)
{
  doc = as(x, "XMLInternalDocument")
  if(is.null(doc))
    getRootNode(x) # skip = skip - getRootNode doesn't have a skip argument
  else
    xmlRoot(doc, skip = skip)
}


     # Get the name of the file/URI for the document.
setGeneric("docName", function(doc, ...) standardGeneric("docName"))

setMethod("docName", "NULL",
           function(doc, ...)
             as.character(NA)
         )

setMethod("docName", "XMLNode",
           function(doc, ...)
             as.character(NA)
         )

setMethod("docName", "XMLHashTreeNode",
           function(doc, ...)
             docName(doc$env, ...)              
         )

docName.XMLInternalDocument =
function(doc, ...)
{
  .Call("RS_XML_getDocumentName", doc, PACKAGE = "XML")
}

setMethod("docName", "XMLInternalDocument", docName.XMLInternalDocument)

docName.XMLInternalNode =
function(doc, ...)
{
  docName(as(doc, "XMLInternalDocument"))
}
setMethod("docName", "XMLInternalNode", docName.XMLInternalNode)

docName.XMLDocument =
function(doc, ...)
{
  doc$doc$file
}
setMethod("docName", "XMLDocument", docName.XMLDocument)

docName.XMLDocumentContent =
function(doc, ...)
{
  doc$file
}

setOldClass("XMLDocumentContent")
setMethod("docName", "XMLDocumentContent", docName.XMLDocumentContent)

setGeneric("docName<-", function(x, value)
                         standardGeneric("docName<-"))



setMethod("docName<-", "XMLInternalDocument",
function(x, value)
{
  .Call("RS_XML_setDocumentName", x, value, PACKAGE = "XML")
  x
})


# See hashTree.R
setMethod("docName<-", "XMLHashTree",
function(x, value)
{
  assign(".doc", value, x)
  x
})


parseXMLAndAdd =
function(txt, parent = NULL, top = "tmp", nsDefs = character())
{
  txt = paste(txt, collapse = "")
  if(!inherits(txt, "AsIs") && length(top) > 0) {
     open = sprintf("%s%s", top,
                            paste(sprintf(' xmlns%s%s="%s"', ifelse(names(nsDefs) != "", ":", ""),
                                                             names(nsDefs),
                                                             nsDefs),
                                       collapse = ""))
     tmp = sprintf("<%s>%s</%s>", open, txt, top)
  } else
     tmp = txt
  
  doc = xmlParse(tmp, asText = TRUE)
  if(!is.null(parent))
     invisible(.Call("R_insertXMLNode", xmlChildren(xmlRoot(doc)), parent, -1L, FALSE,  PACKAGE = "XML"))
  else
     xmlRoot(doc)
}



