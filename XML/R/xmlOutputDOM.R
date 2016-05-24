setClass("XMLOutputStream", "namedList")

setClass("XMLOutputDOM", contains = "XMLOutputStream")
setClass("XMLOutputBuffer", contains = "XMLOutputStream")
setClass("XMLInternalDOM", contains = "XMLOutputStream")

xmlOutputDOM <-
function(tag = "doc", attrs = NULL, dtd = NULL, nameSpace = NULL, nsURI = character(0),
          xmlDeclaration = NULL)
{
 buf <- NULL
 current <- NULL

 startingNode = 1

 if(is.logical(xmlDeclaration) && xmlDeclaration)
    xmlDeclaration = xmlPINode("xml", 'version = "1.0"')
 else if(is.character(xmlDeclaration)) {
     if(length(grep('version *=', xmlDeclaration)) == 0)
       xmlDeclaration = paste(xmlDeclaration, "version='1.0'")
     xmlDeclaration = xmlPINode("xml", xmlDeclaration)
 }

 if(length(dtd)) 
   dtd = paste("<!DOCTYPE", tag, "SYSTEM", ddQuote(dtd[1]), if(length(dtd) > 1) paste("PUBLIC", ddQuote(dtd[2])), ">")

 
 reset <-
 function() {
  buf <<- xmlNode(tag, attrs = attrs, namespace = nameSpace)   

  if(length(nsURI) > 0) {
   names(nsURI) <- paste("xmlns", names(nsURI), sep=":")
   buf$attributes <<- nsURI
  }
  current <<- integer(0)
  invisible(buf)
 }

 reset()


 addTag <- 
 function(tag, ..., attrs=NULL, close=TRUE, namespace=NULL, .children = list(...)) {
   if(missing(namespace))
     namespace <- nameSpace

   addNode(n <- xmlNode(tag, attrs= attrs, namespace = namespace, .children = .children))
   if(close == FALSE) {
     current <<- c(current, xmlSize(getCurrent()))
   }

  invisible(n)
 }

 getCurrentExpr <- 
 function() {
   if(length(current) > 0) {
     p <- seq(2, length=length(current)-1)
     kall <- call("[[", as.name("buf"), current[1])
     for(i in p) {
       kall <- call("[[", kall, current[i])
     }
   } else
     kall <- as.name("buf")
     
   kall
 }

 getCurrent <- function() {
    eval(getCurrentExpr())
 }

  # We want to append this to the currently open (or active)
  # node as defined by `current'
  #   d[[1]] <- append.xmlNode(d[[1]], xmlNode("phone"))
 addNode <- 
 function(node) {
   kall <- getCurrentExpr()
  if(length(current) > 0){
     lhs <- kall
     kall <- call("append.xmlNode", kall, node)
     kall <- call("<<-", lhs, kall)
  } else {
     kall <- call("append.xmlNode", kall, node)
  }

  val <- eval(kall) 
  if(length(current) == 0)
     buf <<- val

  invisible(node)
 }

 addComment <- function(...) {
    addNode(xmlCommentNode(paste(sapply(list(...), as.character), sep="")))
 }

 addCData <- function(text) {
    addNode(xmlCDataNode(text))
 }

 addPI <- function(name, text) {
    addNode(xmlPINode(name, text))  
 }

 addText <- function(text, namespace = "") {
    addNode(xmlTextNode(text, namespace))  
 }    

 closeTag <-
 function(name="", namespace=NULL)  {
    # namespace is ignored since we already have the tag name!
   current <<- current[-length(current)]
 }

 getValue =
  function() {
      # Add DOCTYPE if we have a dtd.
    if(!is.null(xmlDeclaration))
      structure(list(xmlDeclaration = xmlDeclaration, root = buf, doctype = dtd), class = "XMLRDocument")  
    else
       buf
  }
 
 con <- list( value= getValue,
              addTag = addTag,
              addEndTag = function(name){ closeTag(name)},
              closeTag = closeTag,
              reset = reset,
              addNode = addNode,
              add = function(...) {},
              addComment = addComment,
              addPI = addPI,
              addCData = addCData,             
              current = function(){current}
            ) 
  #class(con) <- c("XMLOutputDOM", "XMLOutputStream")
  # con
  ans = new("XMLOutputDOM", con)
  names(ans) = names(con)
  ans
}

xmlRoot.XMLRDocument =
function(x, skip = TRUE, ...)
  x$root

print.XMLRDocument =
function(x, ...)
{
  if(!is.null(x$xmlDeclaration))
    print(x$xmlDeclaration)

  if(!is.null(x$doctype))
    cat(x$doctype, "\n")

  print(x$root, ...)
}
