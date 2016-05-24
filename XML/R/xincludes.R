
setGeneric("findXIncludeStartNodes", 
function(doc, ...)  
{
  standardGeneric("findXIncludeStartNodes")
})

setMethod("findXIncludeStartNodes", "character", 
function(doc, ...)  
{
  findXIncludeStartNodes(xmlParse(doc), ...)
})

setMethod("findXIncludeStartNodes", "XMLInternalDocument", 
function(doc, ...)  
{
  findXIncludeStartNodes(xmlRoot(doc), ...)
})

setMethod("findXIncludeStartNodes", "XMLInternalElementNode", 
function(doc, ...)  
{
  nodes = .Call("R_findXIncludeStartNodes", xmlRoot(doc), PACKAGE = "XML")
  names(nodes) = sapply(nodes, xmlGetAttr, "href", NA)
  nodes
})



findXInclude =
function(x, asNode = FALSE, recursive = FALSE)
{
  while(!is.null(x)) {
    tmp = getSiblingXIncludeStart(x, TRUE)
    if(!is.null(tmp))
      return(fixFindXInclude(tmp, asNode, recursive))

     sib = x
     if(is(sib, "XMLXIncludeStartNode"))
        return(fixFindXInclude(sib, asNode, recursive)) # if(asNode) sib else xmlAttrs(sib))

     x = xmlParent(x)
  }

  fixFindXInclude(x, asNode, recursive)
}

bad.findXInclude = 
 # This version just looks in the left sibling, not all siblings to the left.
function(x, asNode = FALSE, recursive = FALSE)
{
  ans = NULL
  while(!is.null(x)) {
     prev = getSiblingXIncludeStart(x, FALSE)
     if(inherits(prev, "XMLXIncludeStartNode")) {
        ans = prev
        break
     }

     x = xmlParent(x)
  }

  fixFindXInclude(ans, asNode, recursive)
}

fixFindXInclude = 
function(ans, asNode = FALSE, recursive = FALSE)
{
  if(is.null(ans))
    return(NULL)

  if(recursive) {
    tmp = getXIncludePath(ans)
    if(FALSE && grepl(sprintf("^(%s|http:|ftp:)", .Platform$file.sep), tmp))
      tmp
    else
      sprintf("%s%s%s",
               paste(dirname(unique(tmp)), collapse = .Platform$file.sep),
               .Platform$file.sep,
               xmlAttrs(ans))
  } else
    if(asNode) ans else xmlAttrs(ans)["href"]
}

getXIncludePath =
function(node)
{
  x = xmlParent(node)
  ans = character()
  while(!is.null(x)) {
    ans = c(ans, findXInclude(x))
    prev = x
    x = xmlParent(x)
  }
  c(docName(prev), ans)
}

getSiblingXIncludeStart =
function(x, asNode = FALSE)
{
     sib = x
     while(!is.null(sib)) {
       if(inherits(sib, "XMLXIncludeEndNode"))
         return(NULL)
       
       if(inherits(sib, "XMLXIncludeStartNode"))
         return(if(asNode) sib else xmlAttrs(sib))
       
       sib <- getSibling(sib, FALSE)
     }

     NULL
}


getNodePosition =
function(x) {
   if(is.list(x))
     return(sapply(x, getNodePosition))
   
    tmp = getNodeLocation(x)
    sprintf("%s:%d", tmp$file[1], tmp$line)
}


getNodeLocation =
function(node, recursive = TRUE, fileOnly = FALSE)
{
   if(is.list(node))
     return(lapply(node, getNodeLocation, recursive, fileOnly))
            
   fil = findXInclude(node, recursive = recursive)
   if(is.null(fil))
     fil = docName(node)

   if(fileOnly)
      fil[1]
   else
      list(file = fil, line = getLineNumber(node))
}


getLineNumber =
function(node, ...)
{
  if(!is(node, "XMLInternalNode"))
      stop("This must be an C-level/native/internal XML node, i.e. of class 'XMLInternalNode'. Got ", paste(class(node), collapse = ", "))

  .Call("R_getLineNumber", node, PACKAGE = "XML")
}
  
