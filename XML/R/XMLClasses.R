#
# This file contains the definitions of methods
# for operating on the XMLNode objects to make
# the more user-friendly.  Specifically, these
# methods are 
#       print   displays the contents of a node and children
#               as XML text rather than R/S list
#
#       size    returns the number of children
#
#       name    retrieves the tag name
#
#       attrs   retrieves the attributes element of the XML node
#
#    [ and [[   access the children 
#                 (To get at the regular R/S fields in the object, use $
#                    e.g.  node$name, node$attributes, node$value)

#
# In S4/Splus5, we should use the new class mechanism.
#

setS3Method =
function(fun, class) {
  if(!useS4)
    return()
  
  cat("setting method for", fun, class, "\n")
  setMethod(fun, class, get(paste(fun, class, sep = ".")), where = topenv(parent.frame()))
}

if(FALSE)
setOldClass =
function(classes)
{
   ancestors = unique(sapply(classes[-1], oldClass))
   if(length(ancestors)) {
       classes = c(classes[1], ancestors)
       oldClassTable[[ classes[1] ]] <<- ancestors
   }
   methods::setOldClass(classes)
}

# For R 2.7.2 and older.  In 2.8.0, extends() for setOldClass() works 
# better.

oldClassTable =  list(
  "XMLNode" =  c("RXMLAbstractNode", "XMLAbstractNode"),
  "XMLTextNode" = c("XMLNode", "RXMLAbstractNode", "XMLAbstractNode"),
  "XMLPINode" = c( "XMLNode", "RXMLAbstractNode", "XMLAbstractNode") ,
  "XMLProcessingInstruction" = c( "XMLNode", "RXMLAbstractNode", "XMLAbstractNode") ,
  "XMLCommentNode" = c("XMLNode", "XMLTextNode", "RXMLAbstractNode", "XMLAbstractNode"),
  "XMLCDataNode" = c("XMLNode", "RXMLAbstractNode", "XMLAbstractNode"),
  "XMLHashTree" = c("XMLAbstractDocument"),
  "XMLHashTreeNode" = c("RXMLAbstractNode"),
  "XMLDocumentContent" = c(),
  "XMLDocument" = c("XMLAbstractDocument"),
  "XMLHashTree" = c("XMLAbstractDocument"),
  "XMLInternalDocument" = c("XMLAbstractDocument"),
  "HTMLInternalDocument" = c("XMLInternalDocument", "XMLAbstractDocument"),
  "XMLTreeNode" = c("RXMLAbstractNode")    
)



oldClass =
function(class)
{
  if(version$major == "2" && as.integer(version$minor) >= 8) 
     return(unique(c(class, extends(class))))

   c(class, oldClassTable[[ class ]])
}

###############################
# These were in xmlNodes, but need to be defined earlier.

setOldClass("XMLAbstractDocument")
setOldClass(c("XMLInternalDocument", "XMLAbstractDocument"))
setOldClass(c("XMLHashTree", "XMLAbstractDocument"))
setOldClass(c("XMLDocument", "XMLAbstractDocument"))

#XXXsetOldClass(c("HTMLInternalDocument", "XMLInternalDocument")) # , "XMLAbstractDocument"))
setOldClass(c("HTMLInternalDocument", "XMLInternalDocument", "XMLAbstractDocument"))

setOldClass("XMLAbstractNode")
setOldClass(c("RXMLAbstractNode", "XMLAbstractNode"))


# Why do we have to repeat this class inheritance information?
# We don't!
# setOldClass(c("XMLHashTreeNode", "RXMLAbstractNode", "XMLAbstractNode"))
# setOldClass(c("XMLNode", "RXMLAbstractNode", "XMLAbstractNode"))
# setOldClass(c("XMLTextNode", "XMLNode", "RXMLAbstractNode", "XMLAbstractNode"))
# setOldClass(c("XMLPINode", "XMLNode", "RXMLAbstractNode", "XMLAbstractNode"))
# setOldClass(c("XMLCommentNode", "XMLNode", "XMLTextNode", "RXMLAbstractNode", "XMLAbstractNode"))
# setOldClass(c("XMLProcessingInstruction", "XMLNode", "RXMLAbstractNode", "XMLAbstractNode"))
# setOldClass(c("XMLCDataNode", "XMLNode", "RXMLAbstractNode", "XMLAbstractNode"))


setOldClass(c("XMLHashTreeNode", "RXMLAbstractNode"))
setOldClass(c("XMLNode", "RXMLAbstractNode"))
###setOldClass(c("XMLTextNode", "XMLNode"))
methods::setOldClass(c("XMLTextNode", "XMLTextNode", "XMLNode", "RXMLAbstractNode", "XMLAbstractNode" ))
#setOldClass(c("XMLEntitiesEscapedTextNode", "XMLTextNode", "XMLNode", "RXMLAbstractNode", "XMLAbstractNode"))
#setOldClass(c("XMLEntitiesEscapedTextNode", "XMLTextNode"))
setOldClass(c("XMLPINode", "XMLNode"))
setOldClass(c("XMLCommentNode", "XMLNode"))
setOldClass(c("XMLProcessingInstruction", "XMLNode"))
setOldClass(c("XMLCDataNode", "XMLNode"))


setOldClass(c("XMLTreeNode", "XMLNode", "RXMLAbstractNode", "XMLAbstractNode" ))




setOldClass(c("XMLInternalNode", "XMLAbstractNode"))

setOldClass(c("XMLInternalCDataNode", "XMLInternalNode"))
setOldClass(c("XMLInternalPINode", "XMLInternalNode"))
setOldClass(c("XMLInternalCommentNode", "XMLInternalNode"))

setOldClass(c("XMLInternalElementNode", "XMLInternalNode"))
setOldClass(c("XMLInternalTextNode", "XMLInternalNode"))

setOldClass(c("XMLXIncludeStartNode", "XMLInternalNode"))
setOldClass(c("XMLXIncludeEndNode", "XMLInternalNode"))
setOldClass(c("XMLEntityDeclNode", "XMLInternalNode"))
setOldClass(c("XMLAttributeDeclNode", "XMLInternalNode"))
setOldClass(c("XMLDocumentNode", "XMLInternalNode"))
setOldClass(c("XMLDocumentTypeNode", "XMLInternalNode"))
setOldClass(c("XMLDocumentFragNode", "XMLInternalNode"))
setOldClass(c("XMLNamespaceDeclNode", "XMLInternalNode"))

setOldClass(c("XMLAttributeNode", "XMLInternalNode"))


setOldClass(c("XMLDTDNode", "XMLInternalNode"))

setOldClass("XMLNamespace")
setOldClass("XMLNamespaceDefinition")
setOldClass("XMLNamespaceDefinitions")

#setOldClass("XMLInternalDOM")
setOldClass(c("SimplifiedXMLNamespaceDefinitions", "XMLNamespaceDefinitions"))



#setClass("XPathNodeSet", representation(ref = "externalptr"))


setAs("XMLDocument", "XMLInternalDocument",
       function(from) {
         xmlParse(saveXML(from$doc$children$doc))
       })


############

#setMethod("[[", c("XMLInternalElementNode", "numeric") ,
"[[.XMLInternalElementNode" = 
function(x, i, j, ..., exact = NA, namespaces = xmlNamespaceDefinitions(x, simplify = TRUE), addFinalizer = NA)
{
   if(is(i, "numeric"))
     .Call("R_getNodeChildByIndex", x, as.integer(i), addFinalizer, PACKAGE = "XML")
   else
       NextMethod()
}


xmlChildren <-
function(x, addNames = TRUE, ...)
{
 UseMethod("xmlChildren")
}

setGeneric("xmlParent", 
            function(x, ...)
                standardGeneric("xmlParent"))
 



xmlChildren.XMLNode <-
#
# Retrieve the list of children (sub-nodes) within
# an XMLNode object.
#
function(x, addNames = TRUE, ...)
{
  structure(x$children, class = "XMLNodeList")
}


if(useS4) {
setGeneric("xmlChildren", function(x, addNames = TRUE, ...) standardGeneric("xmlChildren"))
setMethod("xmlChildren", "XMLNode", xmlChildren.XMLNode)
}

if(useS4) {
setGeneric("xmlName", function(node, full = FALSE) standardGeneric("xmlName"))
setMethod("xmlName", "XMLCommentNode", function(node, full = FALSE) "comment")

setMethod("xmlName", "XMLNode",
function(node, full = FALSE)
{
  ns  = unclass(node)$namespace
  if(!full || is.null(ns) || ns == "")
    return(node$name)
  
  # 
  if(!is.character(ns)) {
    tmp = ns$id
  } else if(inherits(ns, "XMLNamespace"))
    tmp = names(ns)
  else
    tmp = ns

  if(length(tmp))
     paste(tmp, unclass(node)$name, sep=":")
  else
     unclass(node)$name
})
} else {

xmlName <-
function(node, full = FALSE)
{
  UseMethod("xmlName", node)
}

xmlName.XMLComment <-
function(node, full = FALSE) {
 return("comment")
}

xmlName.XMLNode <-
#
# Get the XML tag name of an XMLNode object
#
function(node, full = FALSE)
{
  ns  = unclass(node)$namespace
  if(!full || is.null(ns) || ns == "")
    return(node$name)
  
  # 
  if(!is.character(ns)) {
    tmp = ns$id
  } else if(inherits(ns, "XMLNamespace"))
    tmp = names(ns)
  else
    tmp = ns

  if(length(tmp))
     paste(tmp, unclass(node)$name, sep=":")
  else
     unclass(node)$name
}

}

xmlAttrs <-
function(node, ...)
{
  UseMethod("xmlAttrs", node)
}


xmlAttrs.XMLNode <-
#
# Get the named list of attributes
# for an XMLNode object.
#
function(node, ...)
{
 node$attributes
}

if(useS4)
setMethod("xmlAttrs", "XMLNode", xmlAttrs.XMLNode)


"[.XMLNode" <-
#
# Extract the  children (sub-nodes) within
# the specified object identified by ...
# and return these as a list
#
function(x, ..., all = FALSE)
{
 obj <- xmlChildren(x) # x$children

 if(all) # "all" %in% names(list(...)) && list(...)[["all"]] == TRUE)
   structure(obj[ names(obj) %in% list(...)[[1]] ], class = "XMLNodeList")
 else
   structure(obj[...], class = "XMLNodeList") # NextMethod("[") 
}


"[[.XMLDocumentContent" <-
function(x, ...) 
{
  x$children[[...]]
}

"[[.XMLNode" <-
#
# Extract the  children (sub-nodes) within
# the specified object identified by ...
#
function(x, ...)
{
 xmlChildren(x)[[...]]
}

names.XMLNode <-
function(x)
{
# names(xmlChildren(x))
 xmlSApply(x, xmlName)
}

"names<-.XMLNode" <-
function(x, value)
{
 names(x$children) <- value
 x
}



length.XMLNode <-
function(x)
{
  xmlSize(x)
}

xmlSize <-
#
# The number of elements within (or length of) a collection
#
function(obj)
{
 UseMethod("xmlSize", obj)
}

xmlSize.XMLDocument <-
function(obj)
{
 return(length(obj$doc$children))
}

xmlSize.default <-
#
# The number of elements within (or length of) a collection
#
function(obj)
{
  length(obj)
}

xmlSize.XMLNode <-
#
# Determine the number of children (or sub-nodes) within an XML node.
#
function(obj)
{
  length(obj$children) 
}


print.XMLComment <- print.XMLCommentNode <-
function(x, ..., indent = "", tagSeparator = "\n")
{
  if(is.logical(indent) && !indent)
    indent <- ""
  
  cat(indent, "<!--", xmlValue(x), "-->", tagSeparator, sep="")
}

print.XMLTextNode <-
function(x, ..., indent = "", tagSeparator = "\n")
{
  if(is.logical(indent) && !indent)
    indent <- ""

  if(inherits(x, "EntitiesEscaped"))
    txt = xmlValue(x)
  else
    txt = insertEntities( xmlValue(x) )
  
  cat(indent, txt, tagSeparator, sep="")
}

setAs("XMLNamespaceDefinitions", "character",
       function(from) {

         if(length(from) == 0)
            return(character())

         ans = structure(sapply(from, function(x) x$uri), names = sapply(from, function(x) x$id))
         if(length(names(ans)) == 0)
           names(ans) = ""
         ans
       })

setAs("character", "XMLNamespaceDefinitions",
       function(from) {

         ids = names(from)
         if(is.null(ids))
           ids = rep("", length(from))

         structure(lapply(seq(along = from),
                           function(i)
                              structure(list(id = ids[[i]], uri = from[i], local = TRUE), class = "XMLNamespace")),
                   class = "XMLNamespaceDefinitions",
                   names = ids)
       })

print.XMLNode <-
#
# displays a node and attributes (and its children)
# in its XML format.
# 
function(x, ..., indent = "", tagSeparator = "\n")
{
 if(length(xmlAttrs(x))) {
   tmp <- paste(names(xmlAttrs(x)),paste("\"", insertEntities(xmlAttrs(x)), "\"", sep=""), sep="=", collapse=" ")
 } else 
   tmp <- ""

 if(length(x$namespaceDefinitions) > 0) {
    k = as(x$namespaceDefinitions, "character")
    ns = paste("xmlns", ifelse(nchar(names(k)), ":", ""), names(k), "=", ddQuote(k), sep = "", collapse = " ")
#   ns <- paste(sapply(x$namespaceDefinitions, 
#                       function(x) {
#                            paste("xmlns", if(nchar(x$id) > 0) ":" else "", x$id, "=", "\"", x$uri, "\"", sep="")
#                       }), collapse=" ")

 } else 
   ns <- ""


   # Add one space to the indentation level for the children.
   # This will accumulate across successive levels of recursion. 
  subIndent <- paste(indent, " ", sep="")
  if(is.logical(indent) && !indent) {
    indent <- ""
    subIndent <- FALSE
  }


    if (length(xmlChildren(x)) == 0) {
      ## Empty Node - so difference is <nodename />
      cat(indent, paste("<", xmlName(x, TRUE), ifelse(tmp != "", 
          " ", ""), tmp, ifelse(ns != "", " ", ""), ns, "/>", tagSeparator, 
          sep = ""), sep = "")
    } else if (length(xmlChildren(x))==1 &&
               inherits(xmlChildren(x)[[1]],"XMLTextNode")) {
      ## Sole child is text node, print without extra white space.
      cat(indent, paste("<", xmlName(x, TRUE), ifelse(tmp != "", 
          " ", ""), tmp, ifelse(ns != "", " ", ""), ns, ">",
          sep = ""), sep = "")
      kid = xmlChildren(x)[[1]]
      if(inherits(kid, "EntitiesEscaped"))
        txt = xmlValue(kid)
      else
        txt = insertEntities( xmlValue(kid) )
      
      cat(txt,sep="")
      cat(paste("</", xmlName(x, TRUE), ">", tagSeparator, 
          sep = ""), sep = "")
    } else {
      cat(indent, paste("<", xmlName(x, TRUE), ifelse(tmp != "", 
          " ", ""), tmp, ifelse(ns != "", " ", ""), ns, ">", tagSeparator, 
          sep = ""), sep = "")
      for (i in xmlChildren(x))
        print(i, indent = subIndent, tagSeparator = tagSeparator)
      cat(indent, paste("</", xmlName(x, TRUE), ">", tagSeparator, 
          sep = ""), sep = "")
    }
}

print.XMLEntityRef <-
function(x, ..., indent="", tagSeparator = "\n")
{
  if(is.logical(indent) && !indent)
    indent <- ""
  
  cat(indent, x$value)
}



print.XMLCDataNode <-
function(x, ..., indent="", tagSeparator = "\n")
{
  if(is.logical(indent) && !indent)
    indent <- ""
  
 cat(indent, "<![CDATA[", tagSeparator, sep = "")
   # Want new lines in value to be replaced by paste("\n", indent, sep="")
 cat(indent, x$value, sep = "")
 cat(indent, "]]>", tagSeparator, sep = "")
}


print.XMLProcessingInstruction <-
function(x, ..., indent="", tagSeparator = "\n")
{
  if(is.logical(indent) && !indent)
    indent <- ""
  
 cat(indent, paste("<?", x$name," ", x$value, "?>", tagSeparator, sep=""), sep = "")
}


xmlElementsByTagName <-
#
# Extract all the sub-nodes within an XML node
# with the tag name `name'.
#
function(el, name, recursive = FALSE)
{
    kids = xmlChildren(el)
    idx  =  (names(kids) == name)
    els  = kids[idx]    
#    idx <-  (names(el$children) == name)
#    els = el$children[idx]

    if(!recursive  || xmlSize(el) == 0)
      return(els)
    
    subs = xmlApply(el, xmlElementsByTagName, name, TRUE)
    subs = unlist(subs, recursive = FALSE)
    
    append(els, subs[!sapply(subs, is.null)])
  }


getDefaultNamespace =
function(doc, ns = xmlNamespaceDefinitions(doc, simplify = simplify), simplify = FALSE)
{
  if(length(ns) == 0)
      return(character())
  i = which(names(ns) == "")
  if(length(i))
     ns[i]
  else
     character()
  
#  val = unlist(sapply(ns, function(x) if(x$id == "") x$uri))
#  if(length(val))
#     val[1]
#  else   
#     character()
}  

matchNamespaces =
  # d = xmlTreeParse("data/namespaces.xml", useInternal = TRUE)
  #   "omg"
  #   c("ns", "omegahat", "r")
  #   c("ns", omegahat = "http://www.omegahat.net", "r")  
  #   c("ns" = "http://www.omegahat.net", "omg" = "http://www.omegahat.net/XML", "r")
  #
  #  Error because z and rs are not defined in the document.
  #  matchNamespaces(d, c("omg", "z", "rs"))
  #
  #
function(doc, namespaces,
          nsDefs = xmlNamespaceDefinitions(doc, recursive = TRUE, simplify = FALSE),
          defaultNs = getDefaultNamespace(doc, simplify = TRUE)
        )
{

    # 3 cases:
    #  i) we are given a single string (e.g. "r") which we use as a prefix for the default namespace
    # ii) we are given a vector of namespaces, but one has no name and we want to use that as the
    #  prefix for the default namespace
    #    e.g.  sum(names(namespaces) == "") == 1)
    # iii) given several elements with no name and we match these to those in the document
    #    if the first doesn't have a match, we use it as the default one.
    # iv) mixture of prefix = uri values and strings with no names.
    #

    # if it is a single "prefix" and we have a default namespace, then map the prefix to the default URI
    # and return.
  if(is.character(namespaces) && length(namespaces) == 1 &&
        is.null(names(namespaces)) && length(defaultNs) > 0) {
       tmp = defaultNs
       names(tmp)[names(tmp) == ""] = namespaces
         # make certain to convert to c(id = url) form from an XMLNamespaceDefinition
       tmp = as(tmp[[1]], "character")
         # if no name, so default namespace, then use the one in namespaces.
       if(length(names(tmp)) == 0 || names(tmp) == "")
         names(tmp) = namespaces
       return(tmp)
   }

    # fix the names so that we have empty ones if we have none at all.
  if(is.null(names(namespaces)))
    names(namespaces) = rep("", length(namespaces))

   # which need to be fixed up.
  i = (names(namespaces) == "")
  
  if(any(i)) {

     # from parameters now: nsDefs = xmlNamespaceDefinitions(xmlRoot(doc), recursive = TRUE)

      # deal with the first one as a special case. If this has no match,
      # we will map it to the default namespace's URI.
    if(i[1] && length(defaultNs) && is.na(match(namespaces[1], names(nsDefs)))) {
      names(namespaces)[1] = namespaces[1]
      namespaces[1] = defaultNs
      msg = paste("using", names(namespaces)[1], "as prefix for default namespace", defaultNs)
      e = simpleWarning(msg)
      class(e) = c("XPathDefaultNamespace", class(e))
      warning(e)
      i[1] = FALSE
    }
    
    if(sum(i) > 0) {
           # So there is at least one namespace without a name.

          # See if there are duplicates
        dups = names(nsDefs)[duplicated(names(nsDefs))]
        tmp = match(namespaces[i], dups)
        if(length(dups) > 0 && any(is.na(tmp)))
          stop("duplicate namespaces, so cannot match namespace prefix(es) ",
                paste(namespaces[i][is.na(tmp)], collapse = ", "),
                " in ", paste(unique(names(nsDefs)), collapse= ", "))

        idx = match(namespaces[i], names(nsDefs))
        if(any(is.na(idx)))
           stop("cannot find defined namespace(s) with prefix(es) ", paste(namespaces[i][is.na(idx)], collapse = ", "))
        names(namespaces)[i] = namespaces[i]
        namespaces[i] = sapply(nsDefs[idx], function(x) x$uri)
      
      #  warning("namespaces without a name/prefix are not handled as you might expect in XPath. Use a prefix")
    } else if(length(defaultNs) == 0)
        stop("There is no default namespace on the target XML document")
  }
  
  if(!is.character(namespaces) || ( length(namespaces) > 1 && length(names(namespaces)) == 0))
     stop("Namespaces must be a named character vector")

  if(length(namespaces) && (length(names(namespaces)) == 0 || any(names(namespaces) == "")))
     warning("namespaces without a name/prefix are not handled as you might expect in XPath. Use a prefix")

  namespaces
}  



getNodeSet =
function(doc, path, namespaces = xmlNamespaceDefinitions(doc, simplify = TRUE), fun = NULL, sessionEncoding = CE_NATIVE,
          addFinalizer = NA,  ...)
{
  xpathApply(doc, path, fun, ...,  namespaces = namespaces, sessionEncoding = sessionEncoding, addFinalizer = addFinalizer)
}


xpathSApply =
function(doc, path, fun = NULL, ... , namespaces = xmlNamespaceDefinitions(doc, simplify = TRUE),
          resolveNamespaces = TRUE, simplify = TRUE, addFinalizer = NA)
{
  answer = xpathApply(doc, path, fun, ..., namespaces = namespaces, resolveNamespaces = resolveNamespaces,
                         addFinalizer = addFinalizer)

    # Taken from sapply
    if (simplify && length(answer) && length(common.len <- unique(unlist(lapply(answer, 
        length)))) == 1) {
        if (common.len == 1) 
            unlist(answer, recursive = FALSE)
        else if (common.len > 1) 
            array(unlist(answer, recursive = FALSE), dim = c(common.len, 
                length(answer)), dimnames = if (!(is.null(n1 <- names(answer[[1]])) & 
                is.null(n2 <- names(answer)))) 
                list(n1, n2))
        else answer
    }
    else answer  
}

xpathApply =
  #
  # the caller can give the same prefixes of the namespaces defined in 
  # the target document as simple names.
  #
  #   xpathApply(d, "/o:a//c:c", fun = NULL, namespaces = c("o", "c"))
  #
  #
function(doc, path, fun = NULL, ... , namespaces = xmlNamespaceDefinitions(doc, simplify = TRUE),
          resolveNamespaces = TRUE, addFinalizer = NA)
{
  UseMethod("xpathApply")
}  

toXMLNode =
  #
  # For taking an internal node and converting it to an R-level node
  #
function(x, ...)
{
  txt = saveXML(x)
  xmlRoot(xmlTreeParse(txt, asText = TRUE))
}

xpathApply.XMLNode =
function(doc, path, fun = NULL, ... , namespaces = xmlNamespaceDefinitions(doc, simplify = TRUE),
          resolveNamespaces = TRUE, addFinalizer = NA, .node = NULL, noMatchOkay = FALSE)
{
  idoc = xmlParse(saveXML(doc), asText = TRUE)
  ans = xpathApply(idoc, path, fun, ..., namespaces = namespaces, resolveNamespaces = resolveNamespaces,
                        .node = .node, noMatchOkay = noMatchOkay)

  # Now convert the results
  if(length(ans)) 
     ans = lapply(ans, toXMLNode)

  ans
}


xpathApply.XMLInternalDocument =
function(doc, path, fun = NULL, ... , namespaces = xmlNamespaceDefinitions(doc, simplify = TRUE),
          resolveNamespaces = TRUE, addFinalizer = NA, .node = NULL, noMatchOkay = FALSE, 
           sessionEncoding = CE_NATIVE, noResultOk = FALSE) # native
{
  path = paste(path, collapse = " | ")

  if(is(namespaces, "list") && all(sapply(namespaces, is, "XMLNamespaceDefinition"))) {
     namespaces = structure(sapply(namespaces, `[[`, "uri"), names = names(namespaces))
  }
    
  if(resolveNamespaces && !inherits( namespaces, "XMLNamespaceDefinitions"))
    namespaces = matchNamespaces(doc, namespaces)
  
  if(!is.null(fun) && !is.call(fun))
    fun = match.fun(fun)

  # create an expression of the form fun(x, ...) and the C code will insert x for each node.
  args = list(...)
  if(length(args))  
    fun = as.call(c(fun, append(1, args)))


      #XXX Match the session encoding c("native" = 0, utf8 = 1, latin1 = 2)
  encoding = if(is.integer(sessionEncoding))
                sessionEncoding
             else
                getEncodingREnum(sessionEncoding)

  ans = .Call("RS_XML_xpathEval", doc, .node, as.character(path), namespaces, fun, encoding, addFinalizer, PACKAGE = "XML")

  if(!noMatchOkay && length(ans) == 0 && length(getDefaultNamespace(xmlRoot(doc))) > 0) {
    tmp = strsplit(path, "/")[[1]]
       # if they have a function call, ignore.
    tmp = tmp[ - grep("\\(", path) ]
    if(length(grep(":", tmp)) != length(tmp) && !noResultOk)
        warning("the XPath query has no namespace, but the target document has a default namespace. This is often an error and may explain why you obtained no results")
  }

  ans
}

xmlDoc =
function(node, addFinalizer = TRUE)
{
  if(!is(node, "XMLInternalElementNode"))
    stop("xmlDoc must be passed an internal XML node")
  
  doc = .Call("RS_XML_createDocFromNode", node, PACKAGE = "XML")
  addDocFinalizer(doc, addFinalizer)
  
  doc
}


# Used to use
#   getDefaultNamespace(doc)


if(FALSE)  {
xpathApply.XMLInternalNode =
  #
  #  There are several situations here.
  #  We/libxml2 needs a document to search.
  #  We have a node. If we use its document, all is fine, but the 
  #  search will be over the entire document.  We may get nodes from other sub-trees
  #  that do not pertain to our starting point (doc).
  #  Alternatively, we can create a new doc with this node as the top-level node
  #  and do the search. But then we end up with new nodes. So if you want to find
  #  nodes in the original document in order to change them rather than just read information
  #  from them, then you will be sorely mistaken when you think your changes have been applied
  #  to the original document.

  #
  # Regardless of what we do, we still have to figure out the adding of the doc attribute.
  #
function(doc, path, fun = NULL, ... , namespaces = xmlNamespaceDefinitions(doc, simplify = TRUE),
          resolveNamespaces = TRUE)
{
  addDocAttribute = FALSE

    # This is a wholesale copy.
  addDocAttribute = TRUE
  doc = xmlDoc(doc, TRUE)
  
     #
     # If the doc is already there, can't we just use that without copying it? Yes.
     # XXX maybe not. Looks like libxml2 starts at the top of the doc again.
     # But if there is no doc for this node, then we create a new doc and
     # put a finalizer on it. But we attach the document as an attribute to each of the
     # the resulting nodes. Then it will be protected from gc() and so will the nodes
     # until each of the nodes are released.
     #XXX???  What if we do a subsequent search on another of these nodes.
     # Then need to add it to the results.
if(FALSE) {
  tmp = as(doc, "XMLInternalDocument")
  addDocAttribute = is.null(tmp)    
  if(is.null(tmp)) {
     if(!is.null(attr(doc, "document")))
       doc = attr(doc, "document")
     else
       doc =  newXMLDoc(node = doc) # if we used  xmlDoc(doc), no finalizer.
  } else
     doc = tmp
}

  ans = xpathApply(doc, path, fun, ..., namespaces = namespaces, resolveNamespaces = resolveNamespaces)

  if(addDocAttribute && length(ans))
    ans = lapply(ans, function(x) { attr(x, "document") = doc; x})
  
  ans
}
} # end if if(FALSE)



getRootNode =
function(node)
{
  p = node
  while(!is.null(xmlParent(p))) 
    p = xmlParent(p)

  p
}

xpathApply.XMLInternalNode =
xpathSubNodeApply =  
  #
  # This allows us to use XPath to search within a sub-tree of the document, i.e. 
  # from a particular node.
  # This is slightly tricky because the libxml2 routines require a document to search.
  # We could copy the nodes to a new document, e.g.
  #    xmlDoc(node)
  # but then the results would be for new nodes, not the original ones.
  # So we would not be able to find nodes and then modify them as we would be modifying the
  # copies.
  #
  #  If this what is desired, use
  #      doc = xmlDoc(node)
  #      xpathApply(doc, xpath, ...)
  #  and then that is what you will get.

  #
  #  In the case that you want the original nodes in the result,
  #  then we have to do a little bit of work. We create a new document
  #  and set the source node as its root node.  We arrange to
  #   a) put the node back where it came from and b) free the document.
  #  So there is no memory leak.
  #
  #  The other thing we must do is to find the 
  #
  # 
  # This version avoids doing any copying of nodes when there is a document already
  # associated with the nodes.
  #
function(doc, path, fun = NULL, ...,
          namespaces = xmlNamespaceDefinitions(doc, simplify = TRUE),
           resolveNamespaces = TRUE, addFinalizer = NA)
{
  path = paste(path, collapse = " | ")
  
  node = doc

  addDocAttribute = FALSE
  createdNewDocument = FALSE
  tmp = as(doc, "XMLInternalDocument")

  putBack =
    function(node, info)  {
      if(!is.null(info$left))
        addSibling(info$left, node)
      else if(!is.null(info$right))
        addSibling(info$right, node, after = FALSE)
      else if(!is.null(info$parent))
        addChildren(info$parent, node)
      else if(!is.null(tmp))
        addChildren(tmp, node)        
    }
  info = list(parent = xmlParent(node),
              left = getSibling(node, after = FALSE),
              right = getSibling(node, after = TRUE))     
  
 
       # The approaches here are to create a new empty document and then set the node
       # to be its root. We don't set the document for each of the sub-nodes but just this
       # top-level node. Then we arrange that when we end this function, we discard the
       # newly created document and put the node back into the original tree in its
       # original position.
       # This involves knowing the parent and the position at which to put the node back into the tree.
       # If it is the top most node, i.e. no parent, then it is simple - just set the parent back to NULL.
       # If it has a parent, but no siblings, just set the parent.
       # And if it has a sibling, put it back next to that sibling.
       #     If it is the first child, put to the left of the sibling.
       #     If it is not, put to the right.
       # Need to make certain the resulting nodes have the original document
       # Use xmlSetTreeDoc rather than node->doc = node as this is recursive.
       # And so this is now all done in the on.exit() via the call to RS_XML_setDocEl()  
#    doc = newXMLDoc(node = node, addFinalizer = getNativeSymbolInfo("R_xmlFreeDocLeaveChildren")$address)

    doc = newXMLDoc(addFinalizer = FALSE)
    parent = xmlParent(node)
    .Call("RS_XML_setRootNode", doc, node, PACKAGE = "XML")
    on.exit({ .Call("RS_XML_unsetDoc", node, unlink = TRUE, parent, TRUE, PACKAGE = "XML")
              .Call("RS_XML_freeDoc", doc, PACKAGE = "XML")
              if(!is.null(tmp)) {
                                        # Need to create a new document with the current node as the root.
                                        # When we are finished, we have to ensure that we put the node back into the original document
                                        # We can use the same mechanism as when we have to create the document from scratch.
               .Call("RS_XML_setDocEl", node, tmp, PACKAGE = "XML") 
              }
              putBack(node, info)
            })

    docName(doc) = paste("created for xpathApply for", path, "in node", xmlName(node))
   

  ans = xpathApply(doc, path, NULL, namespaces = namespaces, resolveNamespaces = resolveNamespaces, addFinalizer = addFinalizer)

  if(length(ans) == 0)
    return(ans)
  
    # now check if the result was actually a descendant of our top-level node for this
    # query. It is possible that it arose from a different sub-tree.

  w = sapply(ans, function(el) .Call("RS_XML_isDescendantOf", el, node, strict = FALSE, PACKAGE = "XML"))

  ans = ans[w]

#  if(FALSE && addDocAttribute && length(ans))
#     ans = lapply(ans, function(x) { attr(x, "document") = doc; x})

#  if(createdNewDocument)
#       # Need to remove the links  from these nodes to the parent.
#    lapply(ans, function(x) .Call("RS_XML_unsetDoc", x, unlink = FALSE, TRUE))

  
  if(!is.null(fun))
     lapply(ans, fun, ...)
  else
     ans
}



if(TRUE)
xpathApply.XMLInternalNode =
function(doc, path, fun = NULL, ...,
          namespaces = xmlNamespaceDefinitions(doc, simplify = TRUE),
           resolveNamespaces = TRUE, addFinalizer = NA)
{
  ndoc = as(doc, "XMLInternalDocument")
  if(is.null(ndoc))
    xpathSubNodeApply(doc, path, fun, ..., namespaces = namespaces, resolveNamespaces = resolveNamespaces, addFinalizer = addFinalizer)
  else
    xpathApply.XMLInternalDocument(ndoc, path, fun, ...,
                                 namespaces = namespaces, resolveNamespaces = resolveNamespaces,
                                 .node = doc, addFinalizer = addFinalizer)
}


xpathApply.XMLDocument =
#xpathApply.XMLNode =  
function(doc, path, fun = NULL, ... , namespaces = xmlNamespaceDefinitions(doc, simplify = TRUE),
          resolveNamespaces = TRUE, .node = NULL, addFinalizer = NA)
{
  txt = saveXML(doc)
  doc = xmlParse(txt, asText = TRUE)
  ans = xpathApply(doc, path, fun, ..., namespaces = namespaces, resolveNamespaces = resolveNamespaces, .node = .node,
                     addFinalizer = addFinalizer)

  if(length(ans)) 
     ans = lapply(ans, toXMLNode)

  ans  
 # stop("XPath expressions cannot be applied to R-level nodes. Use xmlParse() to process the document and then use xpathApply()")
}


# d = xmlTreeParse("data/book.xml", useInternal = TRUE)
# ch = getNodeSet(d, "//chapter")
# xpathApply(ch[[1]], "//section/title", xmlValue)


# d = xmlTreeParse("data/mtcars.xml", useIntern = TRUE); z = getNodeSet(d, "/dataset/variables")
#  xpathApply(z[[1]], "variable[@unit]", NULL, namespaces = character())

getXMLPath =
function(node, defaultPrefix = "ns")
{
  paste(unlist(c("", xmlAncestors(node, xmlName, defaultPrefix))), collapse = "/")
}

xmlAncestors =
function(x, fun = NULL, ..., addFinalizer = NA, count = -1L)
{
  ans = list()
  tmp = x
  while(!is.null(tmp)) {
    if(!is.null(fun))
      ans = c(fun(tmp, ...), ans)
    else
      ans = c(tmp, ans)

    if(count > 0 && length(ans) == count)
      break
    
    tmp = xmlParent(tmp, addFinalizer = addFinalizer)    
  }
  ans
}  
