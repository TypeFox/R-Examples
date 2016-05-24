

GeneralHandlerNames =
  list(SAX =  c("text", "startElement", "endElement", "comment",
                 "startDocument", "endDocument",
                 "processingInstruction", "entityDeclaration",  "externalEntity"),
       DOM =  c("text", "startElement", "comment", "entity", "cdata",
                 "processingInstruction"))

checkHandlerNames =
function(handlers, id = "SAX")
{
  if(is.null(handlers))
    return(TRUE)

  ids = names(handlers)
  i = match(ids, GeneralHandlerNames)
  prob = any(!is.na(i))  
  if(prob) {
    warning("future versions of the XML package will require names of general handler functions to be prefixed by a . to distinguish them from handlers for nodes with those names.  This _may_ affect the ", paste(names(handlers)[!is.na(i)], collapse = ", "))
  }

  if(any(w <- !sapply(handlers, is.function)))
     warning("some handlers are not functions: ", paste(names(handlers[w]), collapse = ", "))  

  !prob
}

xmlEventParse <- 
#
# Parses an XML file using an event parser which calls user-level functions in the
# `handlers' collection when different XML nodes are encountered in the parse stream.
#
# See also xmlParseTree()
#
function(file, handlers = xmlEventHandler(), ignoreBlanks = FALSE, addContext = TRUE,
          useTagName = TRUE, asText = FALSE, trim=TRUE, useExpat = FALSE,
          isURL=FALSE, state = NULL,
          replaceEntities = TRUE, validate = FALSE, saxVersion = 1,
          branches = NULL,  useDotNames =  length(grep("^\\.", names(handlers))) > 0,
          error = xmlErrorCumulator(), addFinalizer = NA) 
{  
  if(libxmlVersion()$major < 2 && !is.character(file))
    stop("Without libxml2, the source of the XML can only be specified as a URI.")


  i = grep("^/", names(handlers))
  if(length(i)) {
    endElementHandlers = handlers[i]
    names(endElementHandlers) = gsub("^/", "", names(endElementHandlers))
    handlers = handlers[ - i]
  } else
    endElementHandlers = list()
    
  
  checkHandlerNames(handlers, "SAX")

  if(validate)
    warning("Currently, libxml2 does support validation using SAX/event-driven parsing. It requires a DOM.")
  else {
      oldValidate = xmlValidity()
      xmlValidity(validate)
      on.exit(xmlValidity(oldValidate))
  }

  if(!any(saxVersion == c(1, 2))) {
     stop("saxVersion must be 1 or 2")
  }

  
  if(inherits(file, "connection")) {
    con = file    
    if(!isOpen(file)) {
      open(file, "r")
      on.exit(close(con))
    }
    file = function(len) {
              txt = readLines(con, 1)
              if(length(txt) == 0) return(txt)
              paste(txt, "\n", sep = "")
           }
  } else if(is.function(file)) {
      # call with -1 to allow us to close the connection
      # if necessary.
    on.exit(file(-1))
  } else {
   if(!asText && missing(isURL)) { 
        # check if this is a URL or regular file.
     isURL <- length(grep("http://",file)) | length(grep("ftp://",file)) | length(grep("file://",file))
   }

   if(isURL == FALSE && asText == FALSE) {
    file = path.expand(file)
    if(file.exists(file) == FALSE)
     stop(paste("File", file, "does not exist "))
   }
   file = as.character(file)
 }

 branches = as.list(branches)
 if(length(branches) > 0 && (length(names(branches)) == 0 || any(names(branches) == "")))
    stop("All branch elements must have a name!")

  old = setEntitySubstitution(replaceEntities)
  on.exit(setEntitySubstitution(old))

  if(!is.function(error))
    stop("error must be a function")
  
  .oldErrorHandler = setXMLErrorHandler(error)
  on.exit(.Call("RS_XML_setStructuredErrorHandler", .oldErrorHandler, PACKAGE = "XML"), add = TRUE)
  
 state <- .Call("RS_XML_Parse", file, handlers,  endElementHandlers, 
                    as.logical(addContext), as.logical(ignoreBlanks),  
                     as.logical(useTagName), as.logical(asText), as.logical(trim), 
                      as.logical(useExpat), state, as.logical(replaceEntities),
                       as.logical(validate), as.integer(saxVersion), branches, as.logical(useDotNames), error,
                        addFinalizer,
                 PACKAGE = "XML")

  if(!is.null(state))
     return(state)
  else
     return(invisible(handlers))
}



xmlStopParser =
function(parser)
{
  if(!inherits(parser, "XMLParserContext"))
    stop("Need an XMLParserContext object for xmlStopParser")

  .Call("RS_XML_xmlStopParser", parser, PACKAGE = "XML")
}  


xmlParserContextFunction =
function(f, class = "XMLParserContextFunction")
{
  class(f) = c(class, class(f))

  f
}



