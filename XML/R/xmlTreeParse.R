xmlSchemaParse =
function(file, asText = FALSE, xinclude = TRUE, error = xmlErrorCumulator())
{
  xmlParse(file, asText = asText, isSchema = TRUE, xinclude = xinclude, error = error)
}

BOMRegExp = "(\\xEF\\xBB\\xBF|\\xFE\\xFF|\\xFF\\xFE)"

xmlTreeParse <- 
   #
   # XML parser that reads the entire `document' tree into memory
   # and then converts it to an R/S object. 
   # Uses the libxml from Daniel Veillard at W3.org. 
   #
   # asText  treat the value of file as XML text, not the name of a file containing
   #       the XML text, and parse that.
   #
   #
function(file, ignoreBlanks = TRUE, handlers = NULL,
           replaceEntities = FALSE, asText = FALSE, trim = TRUE, validate = FALSE, getDTD = TRUE,
           isURL = FALSE, asTree = FALSE, addAttributeNamespaces = FALSE,
           useInternalNodes = FALSE, isSchema = FALSE,
           fullNamespaceInfo = FALSE, encoding = character(),
           useDotNames = length(grep("^\\.", names(handlers))) > 0, 
           xinclude = TRUE, addFinalizer = TRUE, error = xmlErrorCumulator(), isHTML = FALSE, options = integer(),
           parentFirst = FALSE)
{
  isMissingAsText = missing(asText)
  
  if(length(file) > 1) {
    file = paste(file, collapse = "\n")
    if(!missing(asText) && !asText) 
      stop(structure(list(message = "multiple URLs passed to xmlTreeParse. If this is the content of the file, specify asText = TRUE"),
                     class = c("MultipleURLError", "XMLParserError", "simpleError", "error", "condition")))
    asText = TRUE
  }
  
  if(missing(isURL) && !asText) 
    isURL <- length(grep("^(http|ftp|file)://", file, useBytes = TRUE, perl = TRUE))


  if(isHTML) {
    validate = FALSE
    getDTD = FALSE
    isSchema = FALSE
    docClass = "HTMLInternalDocument"
  } else
    docClass = character()

  checkHandlerNames(handlers, "DOM")

  if(missing(fullNamespaceInfo) && inherits(handlers, "RequiresNamespaceInfo"))
    fullNamespaceInfo = TRUE
  

  oldValidate = xmlValidity()
  xmlValidity(validate)
  on.exit(xmlValidity(oldValidate))
  
    # check whether we are treating the file name as
    # a) the XML text itself, or b) as a URL.
    # Otherwise, check if the file exists and report an error.
 if(!asText && isURL == FALSE) {
  if(file.exists(file) == FALSE)
    if(!missing(asText) && asText == FALSE) {
     e = simpleError(paste("File", file, "does not exist"))
     class(e) = c("FileNotFound", class(e))
     stop(e)
    }
    else
     asText <- TRUE
 }

 if(asText && length(file) > 1)
   file = paste(file, collapse = "\n")

 old = setEntitySubstitution(replaceEntities)
 on.exit(setEntitySubstitution(old), add = TRUE)

     # Look for a < in the string.
  if(asText && length(grep(sprintf("^%s?\\s*<", BOMRegExp), file, perl = TRUE, useBytes = TRUE)) == 0) {  # !isXMLString(file) ?
    if(!isHTML || (isMissingAsText && !inherits(file, "AsIs"))) {
      e = simpleError(paste("XML content does not seem to be XML:", sQuote(file)))
      class(e) = c("XMLInputError", class(e))
      (if(isHTML) warning else stop)(e)
    }
  }
  

 if(!is.logical(xinclude)) {
   # if(is(xinclude, "numeric"))
   #  xinclude = bitlist(xinclude) # see bitList.R
   # else
     xinclude = as.logical(xinclude)
 }

 if(!asText && !isURL)
   file = path.expand(as.character(file))

  if(useInternalNodes && trim) {
    prevBlanks = .Call("RS_XML_setKeepBlanksDefault", 0L, PACKAGE = "XML")
    on.exit(.Call("RS_XML_setKeepBlanksDefault", prevBlanks, PACKAGE = "XML"), add = TRUE)
  }

  .oldErrorHandler = setXMLErrorHandler(error)
  on.exit(.Call("RS_XML_setStructuredErrorHandler", .oldErrorHandler, PACKAGE = "XML"), add = TRUE)

  if(length(options))
     options = sum(options)  #XXX coerce to parser options
  
 ans <- .Call("RS_XML_ParseTree", as.character(file), handlers, 
              as.logical(ignoreBlanks), as.logical(replaceEntities),
              as.logical(asText), as.logical(trim), as.logical(validate), as.logical(getDTD),
              as.logical(isURL), as.logical(addAttributeNamespaces),
              as.logical(useInternalNodes), as.logical(isHTML), as.logical(isSchema),
              as.logical(fullNamespaceInfo), as.character(encoding), as.logical(useDotNames),
              xinclude, error, addFinalizer, as.integer(options), as.logical(parentFirst), PACKAGE = "XML")


  if(!missing(handlers) && length(handlers) && !as.logical(asTree))
    return(handlers)

  if(!isSchema && length(class(ans)))
    class(ans) = c(docClass, oldClass(class(ans)))

  if(inherits(ans, "XMLInternalDocument"))
    addDocFinalizer(ans, addFinalizer)
  else if(!getDTD && !isSchema) {
       #??? is this a good idea.
     class(ans) =  oldClass("XMLDocumentContent")
  } 

  ans
}



xmlNativeTreeParse = xmlInternalTreeParse = xmlTreeParse
formals(xmlNativeTreeParse)[["useInternalNodes"]] = TRUE
formals(xmlInternalTreeParse)[["useInternalNodes"]] = TRUE
xmlParse = xmlNativeTreeParse

if(FALSE) {
   # Another approach is to just change the call, as below, but this is tricky
   # to get evaluation of arguments, etc. right.
tmp.xmlInternalTreeParse =
function(file, ignoreBlanks = TRUE, handlers=NULL,
           replaceEntities=FALSE, asText=FALSE, trim=TRUE, validate=FALSE, getDTD=TRUE,
           isURL=FALSE, asTree = FALSE, addAttributeNamespaces = FALSE,
           isSchema = FALSE,
           fullNamespaceInfo = FALSE, encoding = character(),
           useDotNames = length(grep("^\\.", names(handlers))) > 0,  # will be switched to TRUE in the future.
           xinclude = TRUE, addFinalizer = TRUE)
{
  e = sys.call()
  e[[1]] = as.name("xmlTreeParse")
  e[[length(e) + 1]] = FALSE
  names(e)[length(e)] = "useInternalNodes"
  eval(e, parent.env())
}

 # Could try adding this to the top of xmlTreeParse
    # But it won't work with, e.g. lapply(fileNames, xmlInternalTreeParse)
#  if(missing(useInternalNodes) && as.character(sys.call()[[1]]) == "xmlInternalTreeParse")
#     useInternalNodes = FALSE
}


setGeneric("getEncoding",
function(obj, ...)
{
  standardGeneric("getEncoding")
})

setMethod("getEncoding", "ANY", function(obj, ...) NA)

setMethod("getEncoding", "XMLInternalDocument",
            function(obj, ...) {
              .Call("R_getDocEncoding", obj, PACKAGE = "XML")
            })

setMethod("getEncoding", "XMLInternalNode",
            function(obj, ...) {
              .Call("R_getDocEncoding", obj, PACKAGE = "XML")
            })


if(FALSE) {
setMethod("getEncoding", "XMLInternalDOM",
            function(obj, ...) {
               getEncoding(obj)
            })
}

xmlValidity =
function(val = integer(0))
{
  .Call("RS_XML_getDefaultValiditySetting", as.integer(val), PACKAGE = "XML")
}




processXInclude =
function(node, flags = 0L)
  UseMethod("processXInclude")

processXInclude.list =
function(node, flags = 0L)
{
  lapply(node, processXInclude, flags)
}

processXInclude.XMLInternalDocument =
function(node, flags = 0L)

{
  .Call("RS_XML_xmlXIncludeProcessFlags", node, as.integer(flags), PACKAGE = "XML")
}  

processXInclude.XMLInternalElementNode =
function(node, flags = 0L)
{
#  if(xmlName(node) != "include")  # Should check name space also
#    stop("can only process XInclude on include nodes")

  .Call("RS_XML_xmlXIncludeProcessTreeFlags", node, as.integer(flags), PACKAGE = "XML")
} 



