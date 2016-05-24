if(!exists("Sys.setenv", baseenv()))
    Sys.setenv <- get("Sys.putenv", "package:base")


xmlRoot <-
function(x, skip = TRUE, ...)
{
 UseMethod("xmlRoot")
}

xmlRoot.XMLDocument <-
function(x, skip = TRUE,...)
{
#  x$children[[1]]
# x$doc

  xmlRoot(x$doc, skip = skip,...)
}

xmlRoot.XMLDocumentContent <-
function(x, skip = TRUE, ...)
{
  args <- list(...)

  a <- x$children[[1]]
  if(skip & inherits(a, "XMLCommentNode")) {
     which <- sapply(x$children, function(x) !inherits(x, "XMLCommentNode"))
     if(any(which)) {
       which <- (1:length(x$children))[which]
       a <- x$children[[which[1]]]
     } 
  }

 a
}


xmlRoot.HTMLDocument <-
function(x, skip = TRUE, ...)
{
   x$children[[1]]
}

xmlApply <-
function(X, FUN, ...)
{
  UseMethod("xmlApply")
}

xmlSApply <-
function(X, FUN, ...)
{
  UseMethod("xmlSApply")
}

xmlApply.XMLNode <- 
function(X, FUN, ...) { 
  lapply(xmlChildren(X), FUN, ...) 
} 


xmlApply.XMLDocument <-
function(X, FUN, ...)
{
  xmlApply(xmlRoot(X), FUN, ...)
}

xmlSApply.XMLDocument <-
function(X, FUN, ...)
{
  xmlSApply(xmlRoot(X), FUN, ...)
}


xmlSApply.XMLNode <- 
function(X, FUN, ...) { 
  sapply(xmlChildren(X), FUN, ...) 
} 

xmlApply.XMLDocumentContent <-
function(X, FUN, ...)
{
  xmlSApply(X$children, FUN, ...)
}

xmlSApply.XMLDocumentContent <-
function(X, FUN, ...)
{
  xmlSApply(X$children, FUN, ...)
}


xmlValue <- 
function(x, ignoreComments = FALSE, recursive = TRUE, encoding = getEncoding(x), trim = FALSE)
{
 UseMethod("xmlValue")
}

if(useS4)
  setGeneric("xmlValue", function(x, ignoreComments = FALSE, recursive = TRUE, encoding = getEncoding(x))
        standardGeneric("xmlValue"))


xmlValue.XMLNode <- 
function(x, ignoreComments = FALSE, recursive = TRUE, encoding = getEncoding(x), trim = FALSE)
{
 if(recursive && xmlSize(x) > 0) {
   kids = xmlChildren(x)
   if(ignoreComments)
     kids = kids[ !sapply(kids, "XMLCommentNode") ]
   return(paste(unlist(lapply(kids, xmlValue, ignoreComments, trim = trim)), collapse = ""))
 } else if(!recursive && xmlSize(x) > 0) {
        #XXX If !recursive but have text nodes e.g. in the second child.
    i = sapply(xmlChildren(x), inherits, "XMLTextNode")
    if(any(i))
      return(paste(unlist(lapply(xmlChildren(x)[i], xmlValue, ignoreComments, trim = trim)), collapse = ""))
 }
 
 # if(xmlSize(x) == 1) # && (inherits(x[[1]], "XMLTextNode"))
 #    return(xmlValue(x[[1]], ignoreComments))

 if(is.null(x$value))
   character()
 else
   if(trim) trim(x$value) else x$value
}

setS3Method("xmlValue", "XMLNode")

xmlValue.XMLTextNode <- 
function(x, ignoreComments = FALSE, recursive = TRUE, encoding = getEncoding(x), trim = FALSE)
{
  if(!is.null(x$value))
     if(trim) trim(x$value) else x$value
  else
     character(0)
}

setS3Method("xmlValue", "XMLTextNode")

xmlValue.XMLComment <-  xmlValue.XMLCommentNode <-
function(x, ignoreComments = FALSE, recursive = TRUE, encoding = getEncoding(x), trim = FALSE)
{
 if(ignoreComments)
   return("")

  if(!is.null(x$value))
     if(trim) trim(x$value) else x$value
  else
     character(0)
}

setS3Method("xmlValue", "XMLCommentNode")

xmlValue.XMLCDataNode <- 
function(x, ignoreComments = FALSE, recursive = TRUE, encoding = getEncoding(x), trim = FALSE)
{
  if(trim) trim(x$value) else x$value
}

setS3Method("xmlValue", "XMLCDataNode")

xmlValue.XMLProcessingInstruction <- 
function(x, ignoreComments = FALSE, recursive = TRUE, encoding = getEncoding(x), trim = FALSE)
{
  if(trim) trim(x$value) else x$value
}

setS3Method("xmlValue", "XMLProcessingInstruction")

"xmlValue.NULL" =
function(x, ignoreComments = FALSE, recursive = TRUE, encoding = getEncoding(x), trim = FALSE)
          as.character(NA)

#setS3Method("xmlValue", "NULL")

getSibling.XMLInternalNode =
  # Access the next field in the xmlNodePtr object.
  # not exported.
function(node, after = TRUE, addFinalizer = NA, ...)
{
  if(!inherits(node, "XMLInternalNode"))
    stop("can only operate on an internal node")

  .Call("RS_XML_getNextSibling", node, as.logical(after), addFinalizer, PACKAGE = "XML")
}

xmlNamespaceDefinitions <-
function(x, addNames = TRUE, recursive = FALSE, simplify = FALSE, ...)
{
  UseMethod("xmlNamespaceDefinitions")
}

xmlNamespaces = xmlNamespaceDefinitions

xmlNamespaceDefinitions.XMLInternalDocument =
function(x, addNames = TRUE, recursive = FALSE, simplify = FALSE, ...)
{
  r = xmlRoot(x, addFinalizer = FALSE)
  while(!is.null(r) && !inherits(r, "XMLInternalElementNode")) 
     r = getSibling(r, addFinalizer = FALSE)

  if(is.null(r))
    return(if(simplify) character() else NULL)
  
  xmlNamespaceDefinitions(r, addNames, recursive, simplify)
}

xmlNamespaceDefinitions.XMLNode =
  function(x, addNames = TRUE, recursive = FALSE, simplify = FALSE, ...) {
    ans = unclass(x)$namespaceDefinitions

    if(recursive == TRUE) {
                   #      warning("recursive facility not yet implemented.")
      f = function(node) {
            if(!inherits(node, "XMLNode") || xmlName(node) == "")
              return(FALSE)
            ans <<- append(ans, unclass(node)$namespaceDefinitions)
            xmlApply(node, f)
          }
      xmlApply(x, f)
    }

    if(addNames && length(ans) && length(names(ans)) == 0)
        names(ans) = sapply(ans, function(x) x$id)

    if(simplify) {
      if(length(ans) == 0)
        return(character())
      
      ans = structure(sapply(ans, function(x) x$uri),
                      class = c("SimplifiedXMLNamespaceDefinitions", "XMLNamespaceDefinitions"))
    } else if(!is.null(ans))
      class(ans) = "XMLNamespaceDefinitions"

    ans
}

xmlNamespaceDefinitions.XMLInternalNode =
  function(x, addNames = TRUE, recursive = FALSE, simplify = FALSE, ...)
{
    ans = .Call("RS_XML_internalNodeNamespaceDefinitions", x, as.logical(recursive), PACKAGE = "XML")
    if(addNames && length(ans) > 0)
      names(ans) = sapply(ans, function(x) x$id)

    if(simplify) {
      if(length(ans) == 0)
        return(character(0))
      ans = sapply(ans, function(x) x$uri)
      ans = structure(removeDuplicateNamespaces(ans), class = c("SimplifiedXMLNamespaceDefinitions", "XMLNamespaceDefinitions"))
    } else if(!is.null(ans))
       class(ans) = "XMLNamespaceDefinitions"

    ans
  }

setGeneric("getEffectiveNamespaces",
function(node, ...)
 standardGeneric("getEffectiveNamespaces"))

tmp =
function(node, ...)
{  
   ans = xmlNamespaceDefinitions(node)
   merge = function(to, what) {
      i = !(names(what) %in% names(to))
      if(any(i))
        ans[names(what)[i]] <<- what[i] 
   }
     
   tmp = xmlParent(node, manageMemory = FALSE)
   while(!is.null(tmp)) {
      merge(ans, xmlNamespaceDefinitions(tmp))
      tmp = xmlParent(tmp, manageMemory = FALSE)
   }
   ans
}

setMethod("getEffectiveNamespaces", "XMLInternalNode", tmp)
setMethod("getEffectiveNamespaces", "XMLHashTreeNode", tmp)

setMethod("getEffectiveNamespaces", "XMLNode",
           function(node)
               xmlNamespaceDefinitions(node))


removeDuplicateNamespaces =
function(ns)
{
  dups = duplicated(names(ns))
  if(!any(dups))
    return(ns)

  tapply(ns, names(ns),
           function(els) {
             if(length(els) == 1)
               return(TRUE)

             if(length(unique(els)) > 1)
               stop("different URIs for the same name space prefix ", names(els)[1])
             TRUE
           })
  
  ns[!dups]
}  

xmlNamespace <-
function(x)
{
 UseMethod("xmlNamespace")
}


xmlNamespace.XMLNode <-
function(x)
{
  x$namespace
}

#setMethod("xmlNamespace", "character",
xmlNamespace.character = 
function(x) {
    a = strsplit(x, ":")[[1]]
    if(length(a) == 1)
      character()
    else
      a[1]
}
#)

verifyNamespace =
  # Check that the namespace prefix in tag (if any)
  # has a definition in def that matches the definition of the same prefix in node.
function(tag, def, node)
{
   # could have prefix: with no name, but that should never be allowed earlier than this.
  ns = strsplit(tag, ":")[[1]]
  if(length(ns) == 1)
    return(TRUE)

  if(! (ns[1] %in% names(def)) )
     return(FALSE)

  defs = xmlNamespaceDefinitions(node)

  if( defs[[ ns[1] ]]$uri != def[ ns[1] ]) 
      stop("name space prefix ", ns, " does not match ", def[ ns[1] ], " but ", defs[[ ns[1] ]] $uri)

  TRUE
}  


xmlGetAttr <-
  #Added support for name spaces.
function(node, name, default = NULL, converter = NULL, namespaceDefinition = character(),
          addNamespace = length(grep(":", name)) > 0)
{
  a <- xmlAttrs(node, addNamespace)
  if(is.null(a) || is.na(match(name, names(a)))) 
    return(default)

  if(length(namespaceDefinition))
     verifyNamespace(name, namespaceDefinition, node)

  if(!is.null(converter))
    converter(a[[name]])
  else
    a[[name]]
}  


getXInclude =
function(node, parse = FALSE, sourceDoc = NULL)
{
  href = xmlGetAttr(node, "href")
  xpointer = xmlGetAttr(node, "xpointer")

  if(parse) {
     #
     # Perhaps just reload the original document
     # and see what the difference is. Not guaranteed
     # to work since people may have already altered
     # the source document.
    
    if(!is.na(href)) {
       fileName = paste(dirname(docName(sourceDoc)), href, sep = .Platform$file.sep)
       doc = xmlParse(fileName)
    } else
      doc = sourceDoc
    
    if(!is.na(xpointer)) {

    }
  } else 
    c(href = href, xpointer = xpointer)
}  

getInclude =
  #
  #XXX  getXIncludeInfo is not defined!
  #
function(doc, parse = FALSE)
{
  xpathApply(doc, "//xi:include", getXIncludeInfo, parse, docName(doc), doc,
                 namespaces = c(xi="http://www.w3.org/2001/XInclude"))
}  

getXIncludeInfo =
function(node, parse = FALSE, baseURL = character(), doc = NULL)
{

}
