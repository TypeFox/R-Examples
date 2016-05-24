#xmlRoot.HTMLInternalDocument =
  xmlRoot.XMLInternalDocument = 
function(x, skip = TRUE, addFinalizer = NA, ...)
{
  .Call("R_xmlRootNode", x, as.logical(skip), addFinalizer, PACKAGE = "XML")
}

setAs("XMLNode", "XMLInternalNode",
       function(from) {
           con = textConnection("tmp", "w", local = TRUE)
           sink(con)
           on.exit({sink(file = NULL); close(con)})
           print(from)
           
           doc = xmlParse(tmp, asText = TRUE)
           node = xmlRoot(doc)
           removeChildren(node)
           node
        }
      )


setAs("XMLInternalDocument", "character", function(from) saveXML(from))
setAs("XMLInternalDOM", "character", function(from) saveXML(from))


setAs("XMLInternalDocument", "XMLInternalNode",
       function(from) xmlRoot(from))


setAs("XMLInternalNode", "XMLInternalDocument",
        function(from) {
           doc = .Call("R_getXMLNodeDocument", from, PACKAGE = "XML")
           addDocFinalizer(doc, TRUE)
           if(is(doc, "HTMLInternalDocument"))
              class(doc) = c(class(doc), "XMLInternalDocument", "XMLAbstractDocument")
	   doc
      })



setGeneric("free", function(obj) standardGeneric("free"))

setMethod("free", "XMLInternalDocument",
           function(obj) {
              invisible(.Call("R_XMLInternalDocument_free", obj, PACKAGE = "XML"))
           })


addFinalizer =
function(obj, fun, ...)
{
  UseMethod("addFinalizer")
}

addCFinalizer.XMLInternalDocument =
function(obj, fun, ...)
{
  if(missing(fun) || fun == NULL)
    fun = getNativeSymbolInfo("RSXML_free_internal_document")$address
  else if(!is.function(obj)) {

  }
    
  .Call("R_addXMLInternalDocument_finalizer", obj, fun, PACKAGE = "XML")
}   


asRXMLNode =
function(node, converters = NULL, trim = TRUE, ignoreBlanks = TRUE)
{  
   .Call("R_createXMLNode", node, converters, as.logical(trim), as.logical(ignoreBlanks), PACKAGE = "XML")[[1]]
}

"[.XMLInternalDocument" =
function(x, i, j, ..., namespaces = xmlNamespaceDefinitions(x, simplify = TRUE), addFinalizer = NA)
{
  if(is.character(i)) {
    getNodeSet(x, i, ..., addFinalizer = addFinalizer)
  } else
     stop("No method for subsetting an XMLInternalDocument with ", class(i))
}  

"[[.XMLInternalDocument" =
function(x, i, j, ..., exact = NA, namespaces = xmlNamespaceDefinitions(x, simplify = TRUE),
           addFinalizer = NA)
{
  ans = x[i, addFinalizer = addFinalizer]
  if(length(ans) > 1)
    warning(length(ans), " elements in node set. Returning just the first one! (Use [])")
  ans[[1]]
}






xmlName.XMLInternalNode =
function(node, full = FALSE)
{
  ans = .Call("RS_XML_xmlNodeName", node, PACKAGE = "XML")
  if((is.logical(full) && full) || (!is.logical(full) && length(full))) {
    tmp = xmlNamespace(node)
    if(length(tmp) && length(names(tmp)) > 0 && names(tmp) != "")
       ans = paste(names(tmp), ans, sep = ":")
    else if(is.character(full) && full != "")
       ans = paste(full, ans, sep = ":")
  }
  ans
}

if(useS4)
 setMethod("xmlName", "XMLInternalNode", xmlName.XMLInternalNode)


xmlNamespace.XMLInternalNode =
function(x)
{
  .Call("RS_XML_xmlNodeNamespace", x, PACKAGE = "XML")
}



xmlAttrs.XMLInternalNode = 
function(node, addNamespacePrefix = FALSE, addNamespaceURLs = TRUE, ...)
{
  ans = .Call("RS_XML_xmlNodeAttributes",  node, as.logical(addNamespacePrefix), as.logical(addNamespaceURLs), PACKAGE = "XML")
  if(length(attr(ans, "namespaces")))
    ans = new("XMLAttributes", ans) # class(ans) = "XMLAttributes"
  
  ans
}

#setOldClass(c("XMLAttributes", "character"))
setClass("XMLAttributes", contains = "character")

setMethod("show", "XMLAttributes",
           function(object)
              print(unclass(object)))

setMethod('[', c('XMLAttributes', "ANY"),
function(x, i, j, ...)
{
  ans = callNextMethod()
  i = match(i, names(x))
  structure(ans, namespaces = attr(x, "namespaces")[i], class = class(x))
})



xmlChildren.XMLInternalNode =
function(x, addNames = TRUE, omitNodeTypes = c("XMLXIncludeStartNode", "XMLXIncludeEndNode"), addFinalizer = NA, ...)
{
  kids = .Call("RS_XML_xmlNodeChildrenReferences", x, as.logical(addNames), addFinalizer, PACKAGE = "XML")
  
if(length(omitNodeTypes))
    kids = kids[! sapply(kids, function(x) any(inherits(x, omitNodeTypes)) )]   

  structure(kids, class = c("XMLInternalNodeList", "XMLNodeList"))
}


xmlChildren.XMLInternalDocument =
function(x, addNames = TRUE, ...)
{
# .Call("RS_XML_xmlDocumentChildren", x, as.logical(addNames), PACKAGE = "XML")
 xmlChildren.XMLInternalNode(x, addNames, ...)
}


if(useS4) {
setMethod("xmlAttrs", "XMLInternalNode", xmlAttrs.XMLInternalNode)
setMethod("xmlChildren", "XMLInternalNode", xmlChildren.XMLInternalNode)
setMethod("xmlChildren", "XMLInternalDocument", xmlChildren.XMLInternalDocument)
}


xmlSize.XMLInternalNode =
function(obj)
  .Call("RS_XML_xmlNodeNumChildren", obj, PACKAGE = "XML")

"[[.XMLInternalNode" <-
#setMethod("[[", "XMLInternalNode",
function(x, i, j, ..., addFinalizer = NA)
{
  if(inherits(i, "formula")) {
    return(getNodeSet(x, i, if(missing(j)) character() else j, addFinalizer = addFinalizer, ...)[[1]])
  }
  
  if(is.na(i))
     return(NULL)
     # Get the individual elements rather than all the children and then subset those
  return(
       if(is(i, "numeric"))
          .Call("R_getChildByIndex", x, as.integer(i), as.logical(addFinalizer), PACKAGE = "XML")
       else
          .Call("R_getChildByName", x, as.character(i), as.logical(addFinalizer), PACKAGE = "XML")
       )
  
  kids = xmlChildren(x, addFinalizer = addFinalizer)
  if(length(kids) == 0)
    return(NULL)
  
  if(is.numeric(i)) 
     kids[[i]]
  else {
     id = as.character(i)
     which = match(id, sapply(kids, xmlName))
     kids[[which]]
  }
}  



"[.XMLInternalNode" <-
function(x, i, j, ..., addFinalizer = NA)
{
  kids = xmlChildren(x, addFinalizer = addFinalizer)
  if(is.logical(i))
    i = which(i)

  if(is(i, "numeric"))
     structure(kids[i], class = c("XMLInternalNodeList", "XMLNodeList"))
  else {
     id = as.character(i)
     which = match(sapply(kids, xmlName), id)
     structure(kids[!is.na(which)], class = c("XMLInternalNodeList", "XMLNodeList"))     
  }
}  


xmlValue.XMLInternalNode =
function(x, ignoreComments = FALSE, recursive = TRUE, encoding = getEncoding(x), trim = FALSE) #CE_NATIVE)
{

  encoding = if(is.integer(encoding)) 
               encoding 
             else
               getEncodingREnum(encoding)

  if(!recursive) {
     if(xmlSize(x) == 0)
       return(character())

    kids = xmlChildren(x, addFinaliizer = FALSE)
    i = sapply(kids, inherits, "XMLInternalTextNode")
    if(any(i))
      return(paste(unlist(lapply(kids[i], xmlValue, ignoreComments, recursive = TRUE, encoding = encoding)), collapse = ""))
    else
      return(character())
   }
  
  ans = .Call("R_xmlNodeValue", x, NULL, encoding, PACKAGE = "XML") # 2nd argument ignored.
  if(trim)
     trim(ans)
  else
     ans
}

setS3Method("xmlValue", "XMLInternalNode")

setGeneric("xmlValue<-", function(x, ..., value) standardGeneric("xmlValue<-"))

setMethod("xmlValue<-", "XMLInternalTextNode",
           function(x, ..., value) {
             .Call("R_setXMLInternalTextNode_value", x, as.character(value), PACKAGE = "XML")
             x
           })

setMethod("xmlValue<-", "XMLTextNode",
           function(x, ..., value) {
              x$value = as.character(value)
              x
           })

setMethod("xmlValue<-", "XMLAbstractNode",
           function(x, ..., value) {
             if(xmlSize(x) == 0) {
               x = addChildren(x, as.character(value))
             } else if(xmlSize(x) == 1 && any(inherits(x[[1]], c("XMLTextNode", "XMLInternalTextNode")))) {
               #XXX Fix the assignment to children.
               #   should be xmlValue(x[[1]]) = value
               tmp = x[[1]]
               xmlValue(tmp) = as.character(value)
               if(inherits(x[[1]], "XMLTextNode"))
                  x$children[[1]] = tmp
             } else
                 stop("Cannot set the content of a node that is not an XMLInternalTextNode or a node containing a text node")
             x
           })

          

names.XMLInternalNode =
function(x)
  xmlSApply(x, xmlName, addFinalizer = FALSE)

xmlApply.XMLInternalNode =
function(X, FUN, ..., omitNodeTypes = c("XMLXIncludeStartNode", "XMLXIncludeEndNode"), addFinalizer = NA)
{
   kids = xmlChildren(X, addFinalizer = addFinalizer)
   if(length(omitNodeTypes))
     kids = kids[! sapply(kids, function(x) any(inherits(x, omitNodeTypes)) )]   
   lapply(kids, FUN, ...)
}  

xmlSApply.XMLInternalNode =
function(X, FUN, ..., omitNodeTypes = c("XMLXIncludeStartNode", "XMLXIncludeEndNode"), addFinalizer = NA)
{
   kids = xmlChildren(X, addFinalizer = addFinalizer)
   if(length(omitNodeTypes))
     kids = kids[! sapply(kids, function(x) any(inherits(x, omitNodeTypes)) )]
   sapply(kids, FUN, ...)
}  


xmlSApply.XMLNodeSet =
function(X, FUN, ..., omitNodeTypes = c("XMLXIncludeStartNode", "XMLXIncludeEndNode"), addFinalizer = NA)
{
  sapply(X, FUN, ...)
}

xmlApply.XMLNodeSet =
function(X, FUN, ..., omitNodeTypes = c("XMLXIncludeStartNode", "XMLXIncludeEndNode"), addFinalizer = NA)
{
  lapply(X, FUN, ...)
}

getChildrenStrings =
function(node, encoding = getEncoding(node), asVector = TRUE, len = xmlSize(node),
          addNames = TRUE)
{
   encoding = getEncodingREnum(encoding)
   .Call("R_childStringValues", node, as.integer(len), as.logical(asVector), as.integer(encoding),
               as.logical(addNames), PACKAGE = "XML")
}


setMethod("xmlParent", "XMLInternalNode",
function(x, addFinalizer = NA, ...)
{
  .Call("RS_XML_xmlNodeParent", x, addFinalizer, PACKAGE = "XML")
})


newXMLDTDNode <-
function(nodeName, externalID = character(), systemID = character(), doc = NULL, addFinalizer = NA)  
{
  if(length(nodeName) > 1 && missing(externalID))
    externalID = nodeName[2]
  if(length(nodeName) > 2 && missing(systemID))
    systemID = nodeName[3]  

  .Call("R_newXMLDtd", doc, as.character(nodeName), as.character(externalID), as.character(systemID),
          addFinalizer, PACKAGE = "XML")
}

setInternalNamespace =
function(node, ns)
{
  .Call("R_xmlSetNs", node, ns, FALSE, PACKAGE = "XML") # as.logical(append))
}


addDocFinalizer =
function(doc, finalizer)
{
  fun = NULL
  if(is.logical(finalizer)) {
    if(is.na(finalizer) || !finalizer)
      return()
    else
      fun = NULL
  } else {
    fun = finalizer
    if(inherits(fun, "NativeSymbolInfo"))
      fun = fun$address
  }

  if(!is.null(fun) && !is.function(fun) && typeof(fun) != "externalptr")
    stop("need an R function, address of a routine or NULL for finalizer")

  .Call("R_addXMLInternalDocument_finalizer", doc, fun, PACKAGE = "XML")
}

HTML_DTDs =
  c("http://www.w3.org/TR/html4/frameset.dtd",
    "http://www.w3.org/TR/html4/loose.dtd",
    "http://www.w3.org/TR/html4/strict.dtd",
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd",
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd",
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-frameset.dtd",
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"    
    )

newHTMLDoc =
function(dtd = "loose", addFinalizer = TRUE, name = character(),
          node = newXMLNode("html", newXMLNode("head", addFinalizer = FALSE), newXMLNode("body", addFinalizer = FALSE),
                             addFinalizer = FALSE))
{
  if(is.na(dtd) || dtd == "")
     dtd = ""
  else if(tolower(dtd) %in% c("html5", "5"))
     dtd = "5"
  else {
     i = grep(dtd, HTML_DTDs)
     if(length(i)) {
        if(length(i) > 1)
           warning("matched multiple DTDs. Using the first")
        dtd = HTML_DTDs[i[1]]
     } else
        dtd = ""
   } 
     
  
  doc = newXMLDoc(dtd = dtd, isHTML = TRUE, addFinalizer = addFinalizer, node = node)
  doc
}

newXMLDoc <-
#
# Creates internal C-level libxml object for representing
# an XML document/tree of nodes.
#
function(dtd = "", namespaces = NULL, addFinalizer = TRUE, name = character(), node = NULL,
          isHTML = FALSE) 
{
  if(is(dtd, "XMLInternalNode")) {
    dtdNode = dtd
    dtd = character()
  } else
    dtdNode = NULL
  
  ans = .Call("R_newXMLDoc", dtd, namespaces, as.logical(isHTML), PACKAGE = "XML")
  class(ans) = oldClass(class(ans))
  
  addDocFinalizer(ans, addFinalizer)
  
  if(length(name))
     docName(ans) = as.character(name)

  if(length(dtdNode))
     addChildren(ans, dtdNode)

  if(length(node)) {
    if(is.character(node))
	## FIXME: there is no visible 'doc' here
       newXMLTextNode(node, addFinalizer = FALSE, parent = doc)
    else
       addChildren(ans, node)
  }

  ans
}

XMLOptions = new.env()
getOption =
function(name, default = NULL, converter = NULL)
{
  if(!exists(name, XMLOptions, inherits = FALSE)) 
    return(base::getOption(name, default))

   ans = get(name, XMLOptions)

   if(is.function(converter))
     converter(ans)
    else
      ans
}

setOption =
function(name, value)
{
   prev = getOption(name)
   assign(name, value, XMLOptions)
   prev
}  


newXMLNode <-
  ###XXX Note that there is another definition of this in dups.R
  # Which is now elided.
  
  # Create an internal C-level libxml node
  #
  #
  #  It is  possible to use a namespace prefix that is not defined.
  #  This is okay as it may be defined in another node which will become
  #  an ancestor of this newly created one.

  # XXX Have to add something to force the namespace prefix into the node
  # when there is no corresponding definition for that prefix.
  
function(name, ..., attrs = NULL,
         namespace = character(), namespaceDefinitions = character(),
         doc = NULL, .children = list(...), parent = NULL,
         at = NA,
         cdata = FALSE,
         suppressNamespaceWarning = getOption('suppressXMLNamespaceWarning', FALSE), #  i.e. warn.
         sibling = NULL, addFinalizer = NA,
         noNamespace = length(namespace) == 0 && !missing(namespace),
         fixNamespaces = c(dummy = TRUE, default = TRUE)
        )
{
   # determine whether we know now that there is definitely no namespace.
 
    # make certain we have a character vector for the attributes.
 if(length(attrs)) {
     ids = names(attrs)
     attrs = structure(as(attrs, "character"), names = ids)

       # Find any attributes that are actually namespace definitions.
     i = grep("^xmlns", names(attrs))
     if(length(i)) {
         warning("Don't specify namespace definitions via 'attrs'; use namespaceDefinitions")
         namespace = c(namespace, structure(attrs[i], names = gsub("^xmlns:", "", names(attrs)[i])))
         attrs = attrs[ -i]
     }
  } else
     attrs = character()

     # allow the caller to specify the node name as  ns_prefix:name
     # but we have to create it as name and the set the namespace.
  ns = character()  # the namespace prefix
  name = strsplit(name, ":")[[1]]
  if(length(name) == 2) {
     ns = name[1]
     name = name[2]
     noNamespace = FALSE
  }

 if(is.list(parent)) {
    if(length(parent) < 1 ||
        !(is(parent[[1]], "XMLInternalElementNode") || is(parent[[1]], "XMLInternalDocument")))
        stop("incorrect value for parent")
    
    parent = parent[[1]]
  }

   # if there is no doc, but we have a parent which is an XMLInternalDocument, use that.
 if(missing(doc) && !missing(parent) &&
     inherits(parent, "XMLInternalDocument")) {
   doc = parent
   parent = NULL
 }

    # Get the doc from the parent node/document.
 if(is.null(doc) && !is.null(parent))  {
   # doc = as(parent, "XMLInternalDocument")
   doc = if(inherits(parent, "XMLInternalDocument"))
            parent
         else
            .Call("R_getXMLNodeDocument", parent, PACKAGE = "XML")
 }


     # create the node. Let's leave the namespace definitions and prefix till later.

     # xmlSetProp() routine in R_newXMLNode() handles namespaces on the attribute names, even checking them.
  node <- .Call("R_newXMLNode", as.character(name), character(), character(), doc, namespaceDefinitions, 
                   addFinalizer, PACKAGE = "XML")

  if(!is.null(sibling))
     addSibling(sibling, node, after = as.logical(at))
  else if(!is.null(parent))
     addChildren(parent, node, at = at)


 if(TRUE) { # Create the name space definitions here rather than in C code.
      nsDefs = lapply(seq(along  = namespaceDefinitions),
                   function(i) 
                     newNamespace(node, namespaceDefinitions[[i]], names(namespaceDefinitions)[i], set = FALSE)
                   )
      if(length(namespaceDefinitions))
         names(nsDefs) = if(length(names(namespaceDefinitions))) names(namespaceDefinitions) else ""
  } else
      nsDefs = xmlNamespaceDefinitions(node)

       # Now that the namespaces are defined, we can define the attributes which _may_ use them.
  addAttributes(node, .attrs = attrs, suppressNamespaceWarning = suppressNamespaceWarning)

 if(is(namespace, "XMLNamespaceRef")) {
    setInternalNamespace(node, namespace)
 } else if(is.na(noNamespace) || !noNamespace)  {
    ns = getNodeNamespace(ns, nsDefs, node, namespace, noNamespace, namespaceDefinitions, parent, suppressNamespaceWarning)
    if(is.null(ns))
       !.Call("R_setNamespaceFromAncestors", node, PACKAGE = "XML")
#          .Call("R_getAncestorDefaultNSDef", node, TRUE, PACKAGE = "XML")
 }



     # Here is where we set the namespace for this node.
 if(length(ns) && (inherits(ns, c("XMLNamespaceRef", "XMLNamespaceDeclaration")) || (is.character(ns) && ns != ""))) 
     setXMLNamespace( node, ns) # should this be append = FALSE ?


    # Add any children to this node.
 if(length(.children))  {
   if(!is.list(.children))
      .children = list(.children)
   addChildren(node, kids = .children, cdata = cdata, addFinalizer = addFinalizer)
 }

 if(any(fixNamespaces)) { # !is.null(parent)) {
    xmlFixNamespaces(node, fixNamespaces)
   # fixDummyNS(node, suppressNamespaceWarning)
 }

 node
}

xmlFixNamespaces =
function(node, fix)
{

   if(length(fix) == 1)
      fix = structure(rep(fix, 2), names = c("dummy", "default"))

   if(length(names(fix)) == 0)
      names(fix)  = c("dummy", "default")

   
   if(fix["dummy"])
      xmlApply(node, function(x) .Call("R_fixDummyNS", x, TRUE, PACKAGE = "XML"))
   if(fix["default"])
      .Call("R_getAncestorDefaultNSDef", node, TRUE, PACKAGE = "XML")
}


FixDummyNS = 2L
FixDefaultNS = 4L


xmlNamespaceRef =
function(node)
  .Call("R_getXMLNsRef", node, PACKAGE = "XML")



if(FALSE) {
  # Quick check to see if the speed problem in newXMLNode above is in the extra processing
newXMLNode <-
function(name, ..., attrs = NULL,
         namespace = "", namespaceDefinitions = character(),
         doc = NULL, .children = list(...), parent = NULL,
         at = NA,
         cdata = FALSE,
         suppressNamespaceWarning = getOption('suppressXMLNamespaceWarning', FALSE) #  i.e. warn.
        )
{
  node = .Call("R_newXMLNode", name, as.character(attrs), character(), doc, character(), TRUE, PACKAGE = "XML")
  if(!is.null(parent))
     addChildren(parent, node, at = at)
  node
}
}


findNamespaceDefinition =
  #
  # Search up the node hierarchy looking for a namespace
  # matching that prefix.
  #
function(node, namespace, error = TRUE)
{
  ptr = node
  while(!is.null(ptr)) {
     tmp = namespaceDeclarations(ptr, TRUE)
     i = match(namespace, names(tmp))
     if(!is.na(i))
       return(tmp[[i]])
     ptr = xmlParent(ptr)
  }

  if(error)
    stop("no matching namespace definition for prefix ", namespace)

  NULL
}  

setXMLNamespace =
  #
  # Set the specified namespace as the namespace for this
  # node. 
  # namespace can be a prefix in which case we find it in the
  # definition in this node or its ancestors.
  # Otherwise, we expect a name = value character vector giving the
  # prefix and URI and we create a new namespace definition.
  # Alternatively, if you already have the namespace reference object
  # from earlier, you can pass that in.
  # Then we set the namespace on the node. 
function(node,  namespace, append = FALSE)
{
  if(is.character(namespace) && is.null(names(namespace))) 

     namespace = findNamespaceDefinition(node, namespace)

  else if(is.character(namespace))

     namespace = newNamespace(node, namespace)
  else if(!is.null(namespace) && !inherits(namespace, c("XMLNamespaceRef", "XMLNamespaceDeclaration")))
    stop("Must provide a namespace definition, a prefix of existing namespace or a reference to a namespace definition")

  .Call("R_xmlSetNs", node, namespace, FALSE, PACKAGE = "XML") 
}  


setAs("XMLNamespace", "character",
       function(from)
         unclass(from))

setAs("XMLNamespaceDefinition", "character",
       function(from)
         structure(from$uri, names = from$id))


setGeneric("xmlNamespace<-",
            function(x, ..., value)
              standardGeneric("xmlNamespace<-"))

setMethod("xmlNamespace<-", "XMLInternalNode",
            function(x, ..., value) {
               setXMLNamespace(x, value)
               x
            })


setGeneric("xmlNamespaces<-",
            function(x, append = TRUE, set = FALSE, value)
              standardGeneric("xmlNamespaces<-"))


setMethod("xmlNamespaces<-", "XMLNode",
            function(x, append = TRUE, set = FALSE, value) {

                if(inherits(value, "XMLNamespace"))
                  value = as(value, "character")
                else if(is.null(names(value)))
                   names(value) = ""

                    # check for duplicates?
                i = duplicated(names(value))
                if(any(i)) {
                    warning("discarding duplicated namespace prefixes ", paste(names(value)[i], collapse = ", "))
                    value = value[!i]
                }
                
                 if(append) {
                     cur = as(x$namespaceDefinitions, "character")
                     cur[names(value)] = value
                     value = cur
                 }

                x$namespaceDefinitions = as(value, "XMLNamespaceDefinitions")
               
                if(set) 
                  x$namespace = names(value)

               x
            })




setMethod("xmlNamespaces<-", "XMLInternalNode",
            function(x, append = TRUE, set = FALSE, value) {

                value = as(value, "character")
                
                if(is.null(names(value)))
                   names(value) = ""

                    # check for duplicates?
                i = duplicated(names(value))
                if(any(i)) {
                    warning("discarding duplicated namespace prefixes ", paste(names(value)[i], collapse = ", "))
                    value = value[!i]
                }

                if(append) {
                      # Work with existing ones
                   curDefs = namespaceDeclarations(x)
                   i = names(value) %in% names(curDefs)
                   if(any(i)) {
                       warning("discarding duplicated namespace prefixes ", paste(names(value)[i], collapse = ", "))
                       value = value[!i]
                   }
                }

                if(length(value) == 0)
                    # Should worry about the set.
                  return()

                if(length(set) == 1 && set == TRUE && length(value) > 1)
                   set = c(set, rep(FALSE, length(value) - 1))
                else
                   set = rep(set, length.out = length(value))

                
                for(i in seq(along = value))
                  newXMLNamespace(x, value[i], set = set[i])
                
                x
            })
 


newXMLNamespace = newNamespace =
  # Create a new namespace reference object.
function(node, namespace, prefix = names(namespace), set = FALSE)
{
   if(is.null(namespace))
      return(NULL) # XXX
   
   ns <- .Call("R_xmlNewNs", node, namespace, as.character(prefix), PACKAGE = "XML")
   if(set)
      setXMLNamespace(node, ns)
   ns     
}

checkNodeNamespace =
  #
  # can only be checked after we know the parent node, 
  # i.e. after it has been inserted.
  #
function(node, prefix = xmlNamespace(node))
{
  if(length(prefix) == 0 || prefix == "")
    return(TRUE)
  
       # XXX should check that namespace is defined
       # walk the parents.

  okay = FALSE
  p = xmlParent(node)
  while(!is.null(p)) {
    okay = prefix %in% names(xmlNamespaceDefinitions(p))
    if(okay)
      break
  }
  
  if(!okay)
    stop("using an XML namespace prefix '", prefix, "' for a node that is not defined for this node or its node's ancestors")

  TRUE
}  

# Still to do:
#   element, entity, entity_ref, notation
# And more in libxml/tree.h, e.g. the declaration nodes
# 

newXMLTextNode =
  #
  #  cdata allows the caller to specify that the text be 
  #  wrapped in a newXMLCDataNode
function(text,  parent = NULL, doc = NULL, cdata = FALSE, escapeEntities = is(text, "AsIs"),
          addFinalizer = NA)
{
  if(cdata) 
    return(newXMLCDataNode(text, parent, doc, addFinalizer = addFinalizer))

  a = .Call("R_newXMLTextNode", as.character(text), doc, addFinalizer, PACKAGE = "XML")
  if(escapeEntities)
    setNoEnc(a)
  
  if(!is.null(parent))
     addChildren(parent, a)

  a
}

newXMLPINode <-
function(name,  text,  parent = NULL, doc = NULL, at = NA, addFinalizer = NA)
{
  a = .Call("R_newXMLPINode", doc, as.character(name), as.character(text), addFinalizer, PACKAGE = "XML")
  if(!is.null(parent))
    addChildren(parent, a, at = at)
  a  
}

newXMLCDataNode <-
function(text, parent = NULL, doc = NULL, at = NA, sep = "\n", addFinalizer = NA)
{
  text = paste(as.character(text), collapse = "\n")
  a = .Call("R_newXMLCDataNode", doc,  text, addFinalizer, PACKAGE = "XML")
  if(!is.null(parent))
    addChildren(parent, a, at = at)
  a  
}

newXMLCommentNode <-
function(text, parent = NULL, doc = NULL, at = NA, addFinalizer = NA) 
{
  a = .Call("R_xmlNewComment",  as.character(text), doc, addFinalizer, PACKAGE = "XML")
  if(!is.null(parent))
    addChildren(parent, a, at = at)
  a  
}

replaceNodes =
function(oldNode, newNode, ...)
{
  UseMethod("replaceNodes")
}

replaceNodes.list =
function(oldNode, newNode, addFinalizer = NA, ...)
{
 mapply(replaceNodes, oldNode, newNode, MoreArgs = list(addFinalizer = addFinalizer, ...))
}

replaceNodes.XMLInternalNode =
function(oldNode, newNode, addFinalizer = NA, ...)
{
  oldNode = as(oldNode, "XMLInternalNode")
  #XXX deal with a list of nodes.
  newNode = as(newNode, "XMLInternalNode")
  
  .Call("RS_XML_replaceXMLNode", oldNode, newNode, addFinalizer, PACKAGE = "XML")
}  

#
if(FALSE) # This is vectorized for no reason
"[[<-.XMLInternalNode" =
function(x, i, j, ..., value)
{
   if(!is.list(value))
     value = list(value)

  if(is.character(i)) {
     if(length(names(x)) == 0)
         k = rep(NA, length(i))
     else
         k = match(i, names(x))
     
     if(any(is.na(k))) {
           # create a node with that name and text 
         value[is.na(k)] = mapply(function(name, val)
                                    if(is.character(val))
                                         newXMLNode(name, val)
                                    else
                                         val)
     }
     i = k
   }

   replace =  (i <= xmlSize(x))

   if(any(replace)) {
     replaceNodes(xmlChildren(x)[i[replace]], value[replace])
     value = value[!replace]
     i = i[!replace]
   }

   if(length(i))
      addChildren(x, kids = value, at = i)

   x
}  




"[[<-.XMLInternalNode" =
function(x, i, j, ..., value)
{
  if(is.character(i)) {
     if(length(names(x)) == 0)
         k = NA
     else
         k = match(i, names(x))
     
     if(is.na(k) && is.character(value) && !inherits(value, "AsIs")) {
           # create a node with that name and text
        value = newXMLNode(i, value)
     }
     i = k
   }

   replace = !is.na(i) & (i <= xmlSize(x))

   if(replace) 
     replaceNodes(xmlChildren(x)[[i]], value)
   else
     addChildren(x, kids = list(value), at = i)

   x
}  




setNoEnc =
function(node)
{
  if(!is(node, "XMLInternalTextNode"))
    stop("setNoEnc can only be applied to an native/internal text node, not ", paste(class(node), collapse = ", "))

  .Call("R_setXMLInternalTextNode_noenc", node, PACKAGE = "XML")
}



addChildren.XMLInternalNode =
addChildren.XMLInternalDocument =  
  #
  # XXX need to expand/recycle the at if it is given as a scalar
  # taking into account if the subsequent elements are lists, etc.
  #
  # Basically, if the caller specifies at as a scalar
  # we expand this to be the sequence starting at that value
  # and having length which is the total number of nodes
  # in kids.  This is not just the length of kids but
  # the number of nodes since some of the elements might be lists.
  #
function(node, ..., kids = list(...), at = NA, cdata = FALSE, addFinalizer = NA,
          fixNamespaces = c(dummy = TRUE, default = TRUE))
{
  kids = unlist(kids, recursive = FALSE)

  removeNodes(kids)[!sapply(kids, is, "character")]

  if(length(kids) == 1 && inherits(kids[[1]], "XMLInternalNode") && is.na(at)) {
     .Call("R_insertXMLNode", kids[[1]], node, -1L, FALSE, PACKAGE = "XML")
#     return(node)
  } else {

# if(all(is.na(at))) {
#    kids = lapply(kids, as, function(x) if(is.character(x)) newXMLTextNode(x) else as(x, "XMLInternalNode"))
#    .Call("R_insertXMLNodeDirectly", node, kids, PACKAGE = "XML")     
#    return(node)
# }
  
  
  if(!is.na(at)) {

       # if at is the name of a child node, find its index (first node with that name)
    if(is.character(at)) 
      at = match(at, names(node))

    
    if(length(at) == 1) 
       at = seq(as.integer(at), length = sum(sapply(kids, function(x) if(is.list(x)) length(x) else 1)))
    else  # pad with NAs
       length(at) = length(kids)
    
    return(lapply(seq(along = kids),
            function(j) {
               i = kids[[j]]

               if(is.character(i)) 
                 i = newXMLTextNode(i, cdata = cdata, addFinalizer = addFinalizer)

               if(!inherits(i, "XMLInternalNode")) #XX is(i, "XMLInternalNode")
                 i = as(i, "XMLInternalNode")

               if(.Call("R_isNodeChildOfAt", i, node, as.integer(at[j]), PACKAGE = "XML"))
                 return(i)

               if(is.na(at[j]))
                  .Call("R_insertXMLNode", i, node, -1L, FALSE, PACKAGE = "XML")
               else {
                  after = at[j] > 0
                  if(!after)
                     at[j] = 1

                  if(xmlSize(node) < at[j])
                    .Call("R_insertXMLNode", i, node, as.integer(NA), FALSE, PACKAGE = "XML")
                  else
                    .Call("RS_XML_xmlAddSiblingAt", node[[ at[j] ]], i, after, addFinalizer, PACKAGE = "XML") # if at = 0, then shove it in before the sibling.
               }
            }))
  }

  for(j in seq(along = kids)) {
      i = kids[[j]]
      
      if(is.list(i)) {  # can't happen now since we unlist()
         for(k in i)
            addChildren(node, k, addFinalizer = addFinalizer)
      } else {

        if(is.null(i))
           next
     
        if(is.character(i)) 
           i = newXMLTextNode(i, cdata = cdata, addFinalizer = FALSE)


        if(!inherits(i, "XMLInternalNode")) {
           i = as(i, "XMLInternalNode")
         }

        .Call("R_insertXMLNode", i, node, at[j], FALSE, PACKAGE = "XML")

        ns = attr(i, "xml:namespace")
        if(!is.null(ns)) {
           nsdef = findNamespaceDefinition(node, ns)
           if(!is.null(nsdef) && (inherits(nsdef, c("XMLNamespaceRef", "XMLNamespaceDeclaration")) || (is.character(nsdef) && nsdef != ""))) {
             setXMLNamespace( i, nsdef)
             attr(i, "xml:namespace") = NULL
           }
        }
     }
    }
  }

  
  if(!is(node, "XMLInternalDocument") && any(fixNamespaces))
    xmlFixNamespaces(node, fixNamespaces)

  node
}



addSibling =
function(node, ..., kids = list(...), after = NA)
{
  UseMethod("addSibling")
}

addSibling.XMLInternalNode =
function(node, ..., kids = list(...), after = TRUE, addFinalizer = NA)  
{
   #XXX Why add as children?
   if(FALSE && is.na(after))
     addChildren(node, kids = kids, at = NA)
   else {
     lapply(kids,
            function(x) {
              .Call("RS_XML_xmlAddSiblingAt", node, x, as.logical(after), addFinalizer, PACKAGE = "XML")
            })
   }
  
}  



removeNodes =
function(node, free = rep(FALSE, length(node)))
  UseMethod("removeNodes")

removeNodes.default =
function(node, free = rep(FALSE, length(node)))
 NULL

removeNodes.list = removeNodes.XMLNodeList = 
function(node, free = rep(FALSE, length(node)))
{
   if(!all(sapply(node, inherits, "XMLInternalNode"))) {
      warning("removeNode only works on internal nodes at present")
      return(NULL)
   }

    free = as.logical(free)
    free = rep(free, length = length(node))
    .Call("R_removeInternalNode", node, free, PACKAGE = "XML")
}

removeNodes.XMLNodeSet =
function(node, free = rep(FALSE, length(node)))
{
   removeNodes.list(node, free)
}
    

removeNodes.XMLInternalNode =
function(node, free = rep(FALSE, length(node)))
{  
    node = list(node)
    free = as.logical(free)
    .Call("R_removeInternalNode", node, free, PACKAGE = "XML")
}




removeChildren =
function(node, ..., kids = list(...), free = FALSE)
{
  UseMethod("removeChildren")
}

removeChildren.XMLNode =
  #
  #
function(node, ..., kids = list(...), free = FALSE)
{

  kidNames = names(node)
  w = sapply(kids,
              function(i)  {
                orig = i
                if(length(i) > 1)
                  warning("each node identifier should be a single value, i.e. a number or a name, not a vector. Ignoring ",
                           paste(i[-1], collapse = ", "))
                
                if(!inherits(i, "numeric"))
                    i = match(i, kidNames)

                if(is.na(i)) {
                  warning("can't find node identified by ", orig)
                  i = 0
                }
                i
              })

  node$children = unclass(node)$children[ - w ]
  node
}  

removeChildren.XMLInternalNode =
function(node, ..., kids = list(...), free = FALSE)
{  
   # idea is to get the actual XMLInternalNode objects
   # corresponding the identifiers in the kids list.
   # These are numbers, node names or node objects themselves
   # This could be fooled by duplicates, e.g. kids = list(2, 2)
   # or kids = list(2, "d") where "d" identifies the second node.
   # We can put in stricter checks in the C code if needed.
  nodes =  xmlChildren(node)
  nodeNames = xmlSApply(node, xmlName)
  v = lapply(kids,
             function(x)  {
                 if(inherits(x, "XMLInternalNode"))
                   x
                 else if(is.character(x)) {
                   i = match(x, nodeNames)
                   nodes[[i]]
                 } else
                   nodes[[as.integer(x)]]
               })
  
   free = rep(free, length = length(v))
   .Call("RS_XML_removeChildren", node, v, as.logical(free), PACKAGE = "XML")
   node
}



setGeneric("toHTML",
            function(x, context = NULL) standardGeneric("toHTML"))


setMethod('toHTML', 'vector',
            function(x, context = NULL) {
              tb = newXMLNode("table")
              if(length(names(x)) > 0) 
                addChildren(tb, newXMLNode("tr", .children = sapply(names(x), function(x) newXMLNode("th", x))))


              addChildren(tb, newXMLNode("tr", .children = sapply(x, function(x) newXMLNode("th", format(x)))))
              tb
            })

setMethod('toHTML', 'matrix',
            function(x, context = NULL) {
              tb = newXMLNode("table")
              if(length(colnames(x)) > 0) 
                addChildren(tb, newXMLNode("tr", .children = sapply(names(x), function(x) newXMLNode("th", x))))

              rows = sapply(seq(length = nrow(x)), 
                             function(i) {
                               row = newXMLNode("tr")
                               if(length(rownames(x)) > 0)
                                 addChildren(row, newXMLNode("th", rownames(x)[i]))
                               addChildren(row,  .children = sapply(x[i,], function(x) newXMLNode("th", format(x))))
                               row
                             })
              addChildren(tb, rows)

              tb
              })



SpecialCallOperators =
  c("+", "-", "*", "/", "%*%", "%in%", ":")

#XXX Not necessarily working yet! See RXMLDoc
setMethod('toHTML', 'call',
            function(x, context) {
                # handle special operators like +, -, :, ...
              if(as.character(v[[1]]) %in% SpecialCallOperators) {

              }

              v = newXMLNode(x[[1]], "(")
              for(i in v[-1]) 
                 addChildren(v, toHTML( i , context))

              v
            })

setAs("vector", "XMLInternalNode",
      function(from) {
          newXMLTextNode(as(from, "character"))
      })



print.XMLInternalDocument =
function(x, ...)
{
  cat(as(x, "character"), "\n")
}

print.XMLInternalNode =
function(x, ...)
{
  cat(as(x, "character"), "\n")
}  


setAs("XMLInternalNode", "character",
          function(from) saveXML.XMLInternalNode(from))

setAs("XMLInternalTextNode", "character",
          function(from) xmlValue(from))


checkAttrNamespaces =
function(nsDefs, .attrs, suppressNamespaceWarning)
{
      ns = sapply(strsplit(names(.attrs), ":"),
                   function(x)  if(length(x) > 1) x[1] else NA)
      i = which(!is.na(ns))
      m = match(ns[i], names(nsDefs))
      if(any(is.na(m))) {
         f = if(is.character(suppressNamespaceWarning))
                get(suppressNamespaceWarning, mode = "function")
             else
                warning

         f(paste("missing namespace definitions for prefix(es)", paste(ns[i][is.na(m)])))
      }
}

setGeneric("addAttributes",
           function(node, ..., .attrs = NULL, suppressNamespaceWarning = getOption('suppressXMLNamespaceWarning', FALSE), append = TRUE)
             standardGeneric("addAttributes"))

setMethod("addAttributes", "XMLInternalElementNode",
function(node, ..., .attrs = NULL, suppressNamespaceWarning = getOption('suppressXMLNamespaceWarning', FALSE), append = TRUE)
{
   if(missing(.attrs)) 
     .attrs = list(...)

   .attrs = structure(as.character(.attrs), names = names(.attrs))

   if(length(.attrs) == 0)
      return(node)
   
   if(is.null(names(.attrs)) || any(names(.attrs) == ""))
     stop("all node attributes must have a name")


   if(is.character(suppressNamespaceWarning) || !suppressNamespaceWarning) 
      checkAttrNamespaces(getEffectiveNamespaces(node), .attrs, suppressNamespaceWarning)

   if(!append)
     removeAttributes(node, .all = TRUE)
   
   .Call("RS_XML_addNodeAttributes", node, .attrs, PACKAGE = "XML")
   node
})

#if(!isGeneric("xmlAttrs<-"))
 setGeneric("xmlAttrs<-", function(node, append = TRUE, suppressNamespaceWarning = getOption('suppressXMLNamespaceWarning', FALSE), value)
                          standardGeneric("xmlAttrs<-"))

tmp =
function(node, append = TRUE, suppressNamespaceWarning = getOption('suppressXMLNamespaceWarning', FALSE), value)
{
   addAttributes(node, .attrs = value, suppressNamespaceWarning = suppressNamespaceWarning, append = append)
}

setMethod("xmlAttrs<-", "XMLInternalElementNode", tmp)
setMethod("xmlAttrs<-", "XMLNode", tmp)


setMethod("addAttributes", "XMLNode",
           function(node, ..., .attrs = NULL, suppressNamespaceWarning = getOption('suppressXMLNamespaceWarning', FALSE), append = TRUE) {
             if(missing(.attrs)) 
               .attrs = list(...)

             .attrs = structure(as.character(.attrs), names = names(.attrs))

             if(is.null(names(.attrs)) || any(names(.attrs) == ""))
               stop("all node attributes must have a name")

             if(is.character(suppressNamespaceWarning) || !suppressNamespaceWarning) 
                 checkAttrNamespaces(getEffectiveNamespaces(node), .attrs, suppressNamespaceWarning)             

             if(append) {
                i = match(names(.attrs), names(node$attributes))
                if(any(!is.na(i))) {
                  node$attributes[i[!is.na(i)]] =  .attrs[!is.na(i)]
                  .attrs = .attrs[is.na(i)]
                }
                node$attributes = c(node$attributes, .attrs)
             } else
                node$attributes = .attrs
             node
           })

setGeneric("removeAttributes", function(node, ..., .attrs = NULL, .namespace = FALSE,
                                        .all = (length(list(...)) + length(.attrs)) == 0)
                                 standardGeneric("removeAttributes"))


setGeneric("removeXMLNamespaces", 
             function(node, ..., all = FALSE, .els = unlist(list(...))) 
	         standardGeneric("removeXMLNamespaces"))


setMethod("removeXMLNamespaces", "XMLInternalElementNode",
          function(node, ..., all = FALSE, .els = unlist(list(...))) {

      	      if(all)
                .Call("RS_XML_removeAllNodeNamespaces", node, PACKAGE = "XML")
	      else {
                 if(is.character(.els))
                   .els = lapply(.els, function(x) x)
                 .Call("RS_XML_removeNodeNamespaces", node, .els, PACKAGE = "XML")
              }
          })

setMethod("removeAttributes", "XMLInternalElementNode",
#
# The idea here is to remove attributes by name
# We handle the case where these are a simple collection
# of character string identifiers given via the ... or as a character
# vector using, e.g., .attrs = c("a", "b")
#
# Each identifier can be of the form  "name" or "ns:name" giving
# the namespace prefix. We resolve the namespace and 
#

#  If we are dealing with regular attributes (no namespace attributes)
#  then we expect these as a character vector.
#
# The intent of the .namespace argument was originally to indicate that
# we wanted to remove the namespace definition. It appears that libxml2 does
# not support that. (And it would seem that this is a real pain as the xmlNsPtr
# objects can be shared across numerous places in a linked list, so it would 
# be very difficult to remove it from one node.)
#
#
#
function(node, ..., .attrs = NULL, .namespace = FALSE,
          .all = (length(list(...)) + length(.attrs)) == 0)
{
   if(missing(.attrs)) 
     .attrs = list(...)
   
   .attrs = as.character(.attrs)
  
   if(.all) {
     if(length(list(...)) || length(.attrs))
         stop(".all specified as TRUE and individual values specified via .../.attrs")

         # Use the integer indices to identify the elements.
     .Call("RS_XML_removeNodeAttributes", node, seq(along = xmlAttrs(node)), FALSE, PACKAGE = "XML")
     return(node)
   }
   

   if(is(.namespace, "XMLNamespaceDeclaration"))
     .namespace = list(.namespace)
#XXX   

   tmp = strsplit(.attrs, ":")
   prefix = sapply(tmp, function(x) if(length(x) > 1) x[1] else "")
   ids = sapply(tmp, function(x) if(length(x) == 1) x[1] else x[2])   

   if(any(prefix != "") && is.logical(.namespace))
     .namespace = TRUE
   
   if(is.logical(.namespace) && .namespace) {
     ns = namespaceDeclarations(node, TRUE)
     # need to create a list with the elements corresponding to the
     # (potentially repeated) ns elements

     i = match(prefix, names(ns))
     ns = ns[i]
     names(ns) = gsub("^.*:", "", .attrs) # or ids from above
     
     .attrs = ns
   }

   .Call("RS_XML_removeNodeAttributes", node, .attrs, .namespace, PACKAGE = "XML")
   node
})


setMethod("removeAttributes", "XMLNode",
function(node, ..., .attrs = NULL, .namespace = FALSE,
         .all = (length(list(...)) + length(.attrs)) == 0)
{
   a = node$attributes

   if(missing(.attrs)) 
     .attrs = list(...)
   
   .attrs = as.character(.attrs)

  if(.all) {
    if(length(.attrs))
      stop("Both individual attribute names and .all specified")
    node$attributes = character()
    return(node)
  }

  i = match(.attrs, names(a))
  if(any(is.na(i)) ) 
     warning("Can't locate attributes ", paste(.attrs[is.na(i)], collapse = ", "), "in XML node ", node$name)

  a = a[is.na(i)]
  
  node$attributes <- a
  node
})

#xmlNamespaceDefinitions =  # ??? added this but overrides other S3 generic.
namespaceDeclarations =
function(node, ref = FALSE, ...)
{
  .Call("RS_XML_getNsList", node,  as.logical(ref), PACKAGE = "XML")
}  


"xmlName<-" =
function(x, value)
{
   UseMethod("xmlName<-")
}  

"xmlName<-.XMLNode" <-
function(x, value)
{
   x$name <- value
   x
}

"xmlName<-.XMLInternalElementNode" <-
function(x, value)
{
   # we could handle a new namespace by accepting value as
   # a character vector with a name
   # e.g.   c(r:array = 'http://www.r-project.org')
   # Alternatively, just define the namespace on the node _before_
   # changing the name.
   id = names(value)
   if(!is.null(id) && length( (tmp <- strsplit(id, ":")[[1]])) > 1) {
       names(value) = tmp[1]
       newXMLNamespaces(x, .values = as(value, "character"))
       value = id
   }
  
   .Call("RS_XML_setNodeName", x, value, PACKAGE = "XML")

   x
}  



newXMLNamespaces =
  # allow for multiple namespaces
  # and also allow for "r:value"
  #
  #  newXMLNamespaces(node, r = "http://www.r-project.org", ...)
  #
function(node, ..., .values = list(...))
{
  ids = names(.values)
  ans = lapply(ids, function(id)
                      newNamespace(node, id, as.character(.values[[id]])))

  names(ans) = ids
  ans
}




xmlNodeMatch =
function(x, table, nomatch = NA_integer_)
{
  .Call("R_matchNodesInList", x, table, as.integer(nomatch), PACKAGE = "XML")
}  


setGeneric("xmlClone",
function(node, recursive = TRUE, addFinalizer = FALSE, ...)
           {
              oclass = class(node)
              ans = standardGeneric("xmlClone")
              if(!isS4(node))
                class(ans) = oclass
              ans
           })

setMethod("xmlClone", "XMLInternalDocument",           
function(node, recursive = TRUE, addFinalizer = NA, ...)
{
  ans = .Call("RS_XML_clone", node, as.logical(recursive), addFinalizer, PACKAGE = "XML")

  addDocFinalizer(ans, addFinalizer)
  ans
})

setMethod("xmlClone", "XMLInternalNode",           
function(node, recursive = TRUE, addFinalizer = FALSE, ...)
{  
  ans = .Call("RS_XML_clone", node, as.logical(recursive), addFinalizer, PACKAGE = "XML")
})








ensureNamespace =
  #
  # Idea is to make certain that the root node has definitions for the specified
  # namespaces.  The caller specifies the named vector of interest.
  # If the URL already exists, we return the corresponding prefix.
  # 
  #
  #  Returns the prefixes in the documents that correspond to the
  #  namespace definitions
  #
function(doc, what)
{
  if(is(doc, "XMLInternalDocument"))
     node = xmlRoot(doc)
  else
     node = doc
  
  defs = xmlNamespaceDefinitions(xmlRoot(doc), simplify = TRUE)
  i = match(what, defs)
  w = is.na(i)

  if(any(w)) {
     sapply(names(what)[w], function(id) newXMLNamespace(node, what[id], id))
     names(what)[w]
  } else
     names(defs)[i]
}


"xmlParent<-" =
 function(x, ..., value) {
  addChildren(value, ..., kids = list(x))                
}


setOldClass("XMLNamespaceRef")
setAs("XMLNamespaceRef", "character",
       function(from) {
 	.Call("R_convertXMLNsRef", from, PACKAGE = "XML")
       })




xmlSearchNs =
function(node, ns, asPrefix = TRUE, doc = as(node, "XMLInternalDocument"))
{

 .Call("R_xmlSearchNs", doc, node, as.character(ns), as.logical(asPrefix), PACKAGE = "XML")
}
