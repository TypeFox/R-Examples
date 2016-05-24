namespaceNodeHandlers =
  #
  # This is to manage a collection of node handlers 
  # which have handlers for nodes with the same name but in an
  # different namespace.
  # For example, suppose we have a node array  in an R representation
  # and another in Matlab but these are differentiated by the namespace.
  # r:array and m:array.
  # This function arranges to invoke the correct handler when a node is encountered.
  #
  #  namespaceNodeHandlers("r:array" = function(...) ...,
  #                        "m:array" = function(...) ...,
  #                        other = function(...) ...)
  #
  #
  #  namespaceNodeHandlers("r:array" = function(...) ...,
  #                        "m:array" = function(...) ...,
  #                        other = function(...) ...,
  #                        nsDefs= c(r = "http://www.r-project.org",
  #                                  m = "http://www.mathworks.com")  )
  #
  #
  # If there is no handler for a given node, then call the default one
  # i.e. .startElement or startElement
  #
function(..., .handlers = list(...), nsDefs = NULL, useDotNames = TRUE)
{
   # Get the node name and namespace prefix/alias for each of the handler names
  tmp = strsplit(names(.handlers), ":")

  prefix = sapply(tmp, function(x) if(length(x) > 1) x[1] else "")
  nodeNames = sapply(tmp, function(x) if(length(x) > 1) x[2] else x[1])

   # Now, find out which ones are duplicated.
  w = duplicated(nodeNames)
  if(!any(w))
    return(.handlers)

  dups = nodeNames[w]
  
    # Now, take out the handler functions that have the same node name
    # as any other.
  w = nodeNames %in% dups 
  nsHandlers = .handlers[ w ]
    # and remove them from .handlers
  .handlers = .handlers[ !w ]


    # This function will act as the proxy for doing the dispatch for a particular node.
  generalHandler =
  function(node, ...) {

       # Get the node name and the namespace prefix.
      id = xmlName(node)

      ns = xmlNamespace(node)

      if(is.null(ns))
        ns = ''
      
      if(length(nsDefs)) {
          # get the namespace definition from the node
          # and its URI and then match this to the nsDefs
          # That gives us the prefix to use

        i = (ns == nsDefs)
        if(!any(i) && ns != '')
          ns = character() #stop("can't match namespace '", as.character(ns), "' in ", paste(nsDefs, collapse = ", "))
        else
          ns = names(nsDefs)[i]
      }

      if(length(ns) && ns != "")
        tmp = paste(ns, id, sep = ":")
      else
        tmp = id

      f = nsHandlers[[ tmp ]]

        # if we didn't find a handler, use the startElement one
      if(is.null(f)) 
        f = .handlers[[ if(useDotNames) '.startElement' else 'startElement' ]]

        # if we have a handler, call it.
        # Otherwise, just return the node... after all that!
      if(!is.null(f))       
        f(node, ...)
      else
        node
    }


   .handlers[ dups ] = rep(list(generalHandler), length(dups))

   class(.handlers) <- "XMLNamespaceNodeHandlers"

   if(length(nsDefs))
     class(.handlers) <- c(class(.handlers), "RequiresNamespaceInfo")

   .handlers
}  
