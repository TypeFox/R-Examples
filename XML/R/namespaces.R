setGeneric("simplifyNamespaces",
             function(doc, ...)
                standardGeneric("simplifyNamespaces"))

setMethod("simplifyNamespaces", "character",
            function(doc, ...) {
               pdoc = xmlParseDoc(doc, NSCLEAN)
               simplifyNamespaces(pdoc, ...)
            })   

xmlCleanNamespaces =
  #
  # @eg xmlParse("~/GitWorkingArea/XML/inst/exampleData/redundantNS.xml")
  #
  # ?Should we write the result to a file if we are given a file?
  #
  #
function(doc, options = integer(), out = docName(doc), ...)
{
   if(is(doc, "XMLInternalDocument"))
      doc = saveXML(doc)
   options = unique(c(options, NSCLEAN))
   newDoc = xmlParse(doc, ..., options = options)

   if(is.logical(out))
     out =  if(out)  docName(doc) else character()
   
   if(is.character(out) && length(out))
     saveXML(newDoc, out)
   else
     newDoc
}

setMethod("simplifyNamespaces", "XMLInternalDocument",
            function(doc, alreadyCleaned = FALSE, ...) {           

                 # find all the nodes, but discard the root node.
               allNodes = getNodeSet(doc, "//node()") # [-1]
               root = xmlRoot(doc)

                 # For each node, get its namespace definitions,
                 # and then zoom in on the nodes that have namespace definitions.
               nsDefs = lapply(allNodes, xmlNamespaceDefinitions, simplify = TRUE)
               w = sapply(nsDefs, length) > 0
               tmp = structure(unlist(nsDefs[w]), names = sapply(nsDefs[w], names))


               d = data.frame(uri = tmp, prefix = names(tmp), stringsAsFactors = FALSE)
   
               multi = unlist(by(d, d$prefix, function(x) if(length(unique(x$uri)) == 1) character() else x$prefix[1]))
               if(length(multi))
                 d = d[ ! (d$prefix %in% multi), ]

               #  Now we can move these namespace definitions to the top.
               #
               #
               #
               by(d, nsDefs,
                               function(x) {
                                    u = unique(x$prefix)
                               })
             
                  # remove the 
               sapply(allNodes[w], removeXMLNamespaces)
               
               nsDefs
})



getNodeNamespace =
 # Figure out what namespace to use for this node and return a reference to that
 # namespace definition object in C (a xmlNsPtr)
function(ns, nsDefs, node, namespace, noNamespace, namespaceDefinitions = NULL, parent = NULL,
          suppressNamespaceWarning = FALSE)
{
  if(noNamespace)
    return(NULL)

 if(is.character(namespace) && length(namespace) && !is.na(namespace) && namespace == "") {
   if(length(namespaceDefinitions) == 0)
       return(findNamespaceDefinition(node, ""))
 }
  
 if((is.list(namespace) || is.character(namespace)) && length(namespace) > 0) {
           # a single element with no name so this is the prefix.
    if(length(namespace) == 1 && length(names(namespace)) == 0) {
      if(namespace %in% namespaceDefinitions) {
         i = match(namespace, namespaceDefinitions)
         ns = nsPrefix = names(namespaceDefinitions)[i]
      } else  if(namespace != "") {
        ns = nsPrefix = namespace
      }
    } else {
           # we have names and/or more than one element. So these are namespace definitions
      if(length(names(namespace)) == 0) 
        names(namespace) <- rep("", length(namespace))

      if(length(namespace) > 1 && !is.na(match(namespace[1], names(namespace)[-1]))) {
        if(length(ns)) 
          warning("ignoring first element of namespace and using prefix from node name, ", ns)
        else {
          ns = namespace[1]
          namespace = namespace[-1]
        }
      }

      if(length(namespace) > 1 && sum(names(namespace) == "") > 1)
        warning("more than one namespace to use as the default")

      nsDefs = lapply(seq(along = namespace),
                      function(i) {
                        prefix = names(namespace)[i]

                        newNamespace(node, namespace[[i]], prefix)
                                # Don't set the namespace. This is just a definition/declaration for
                                # this node, but not necessarily the namespace to use for this node.
                                # We set this below
                      })
      names(nsDefs) = names(namespace)
    }
  }

          # Now handle the prefix for this node.
    if(length(ns)) {
      i = match(ns, names(nsDefs))
      if(is.na(i)) {
         if(!is.null(parent)) 
           ns = findNamespaceDefinition(node, ns)
         else {
#	     raiseNsWarning(ns, suppressNamespaceWarning)
#            attr(node, "xml:namespace") = ns
#            ns = NULL
   	    ns = newNamespace(node, character(), ns)
          }
          if(!inherits(ns, "XMLNamespaceRef"))
             ns <- newNamespace(node, ns, "")
      } else
         ns <- nsDefs[[i]]
    } else  {
      i = match("", names(nsDefs))
      ns = if(is.na(i)) NULL else nsDefs[[i]]

         # if now namespace and we have a parent, use its namespace
         # if it has a namespace
      if(!noNamespace && length(ns) == 0 && length(parent) > 0) {
          ns = xmlNamespaceRef(parent)
          if(!is.null(ns) && names(as(ns, "character")) != "")
            ns = NULL
      }
    }

   ns
}

raiseNsWarning = 
function(ns, suppressNamespaceWarning)
{
   if(is.character(suppressNamespaceWarning))
       f = get(suppressNamespaceWarning, mode = "function")
   else if(is.logical(suppressNamespaceWarning)) {
     if(!suppressNamespaceWarning) 
        f = warning 
     else
        return(NULL)
 } else
    f = function(...) {}
                
    f("cannot find namespace definition for '", ns, "' because the node is not in a document and there are no matching local namespace definitions for this node")
}


fixDummyNS = 
function(node, suppressNamespaceWarning = getOption('suppressXMLNamespaceWarning', FALSE))
{

return(NULL)
   nodes = getNodeSet(node, "//*[./namespace::*[. = '<dummy>']]", addFinalizer = FALSE)
   lapply(nodes, completeDummyNS, suppressNamespaceWarning = suppressNamespaceWarning)
}

completeDummyNS = 
function(node, suppressNamespaceWarning = getOption('suppressXMLNamespaceWarning', FALSE))
{
  if(is.null(xmlParent(node)))
     return(FALSE)

  prefix = names(xmlNamespace(node))
  ns = findNamespaceDefinition(xmlParent(node), prefix, error = FALSE)
  if(is.null(ns))
    raiseNsWarning(prefix, suppressNamespaceWarning)
#    (if(suppressNamespaceWarning) warning else stop)("can't find namespace definition for prefix ", prefix)
  else {
      # remove the current namespace definition and kill it.
    .Call("R_replaceDummyNS", node, ns, prefix, PACKAGE = "XML")
    # setXMLNamespace(node, ns)
  }
}
