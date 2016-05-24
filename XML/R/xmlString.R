setClass("XMLString", contains = "character")
xml =
function(x)  
{
  new("XMLString", x)
}

isXMLString =
function(str)
{
   is(str, "XMLString") ||  length(grep("<([a-zA-Z]+:)?[a-zA-Z]+(/?>| [a-zA-Z]+=[\"'])", str)) > 0
}


RXMLNamespaces = 
  c(r = "http://www.r-project.org", 
    rh = "http://www.r-project.org/help")

xmlParseString =
  #
  # have to do some trickery with garbage collection to avoid parsing
  # the tree and handing back a node within that  that will be freed
  # when the top-level document is GC'ed
  #
  # If the caller gives us the target document, we parent the nodes
  # into that and all is well.
  # Otherwise, we have trouble and potential leak at present as we
  # explicitly kill off the finalizer.
  #
  # This should be cured now in the XML package.
  #
  
function(content, doc = NULL, namespaces = RXMLNamespaces, clean = TRUE, addFinalizer = NA)
{
  f =
    function(cdata = FALSE)
       newXMLNode("para", newXMLTextNode(content, cdata = cdata, doc = doc), doc = doc)
  
    # If the user has told us explicitly that this is not XML but raw text
    # then put it into a para enclosed within CDATA, just in case.
  if(inherits(content, "AsIs"))
       return(f(TRUE))
  if(!isXMLString(content))
      return(f(TRUE))
  content = as(content, "XMLString")

  ns = paste(paste("xmlns", names(RXMLNamespaces), sep = ":"),
              sprintf('"%s"', RXMLNamespaces), sep = "=", collapse = " ")
  txt = paste('<para ', ns, '>', content, "</para>", sep = "")

  local.doc = tryCatch(xmlParse(txt, addFinalizer = addFinalizer),   # addFinalizer = !inherits(doc, "XMLInternalDocument")
                       error = function(e) e)

  if(inherits(local.doc, "condition"))
    return(f(TRUE))

  tmp = xmlRoot(local.doc)
  if(xmlSize(tmp) == 1)  # inherits(tmp[[1]], "XMLInternalElementNode") && xmlName(tmp[[1]]) == "para")
    tmp = tmp[[1]]

  if(clean) 
    removeXMLNamespaces(tmp, .els = names(namespaces))

  
   # XXX
  if(inherits(doc, "XMLInternalDocument")) {
    manageMemory = manageMemory_p(addFinalizer)
    .Call("RS_XML_copyNodesToDoc", tmp, doc, addFinalizer, PACKAGE = "XML")
  } else
    tmp
}  
