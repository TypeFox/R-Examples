xmlErrorCumulator =
function(class = "XMLParserErrorList", immediate = TRUE)
{
  messages = character()
  function(msg, ...) {
         # curently discards all the extra information.
    if(length(grep("\\\n$", msg)) == 0)
      paste(msg, "\n", sep = "")
    
     if(immediate)
       cat(msg)
    
     if(length(msg) == 0) {
          # collapse into string. Probably want to leave as separate elements of a character vector.
          # Make into real objects with the ... information.
       e = simpleError(paste(1:length(messages), messages, sep = ": ",collapse = ""))
       class(e) = c(class, class(e))
       stop(e)
     }

     messages <<- c(messages, msg)
  }
}


xmlStop =
  #
  # Never used anymore.
  # Related to the non-structed error handling.
function(msg, class = "XMLParserError")
{
  err = simpleError(msg)
  class(err) = c(class , class(err))
  stop(err)
}

makeXMLError = 
function(msg, code, domain, line, col, level, filename, class = "XMLError")
{
  err = simpleError(msg)
  err$code = getEnumValue(code, xmlParserErrors)
  err$domain = getEnumValue(domain, xmlErrorDomain)
  err$line = line
  err$col = col
  err$level = getEnumValue(level, xmlErrorLevel)
  err$filename = filename
  
  class(err) = c(class, class(err))
  err
}

htmlErrorHandler = 
function(msg, code, domain, line, col, level, filename, class = "XMLError")
{
  e = makeXMLError(msg, code, domain, line, col, level, filename, class)
  dom = names(e$domain)
  class(e) = c(names(e$code),
               sprintf("%s_Error", gsub("_FROM_", "_", dom)), 
               class(e))

  if(e$code == xmlParserErrors["XML_IO_LOAD_ERROR"])
    stop(e)
}

xmlStructuredStop =
function(msg, code, domain, line, col, level, filename, class = "XMLError")
{
  err = makeXMLError(msg, code, domain, line, col, level, filename, class)

  stop(err)
}  


xmlErrorFun =
function()
{
  errors = list()
  h = function(msg, code, domain, line, col, level, filename) {
    if(length(msg) == 0)
      return(TRUE)


   err = list(msg = msg, code = code,
              domain  = domain, line = line,
              col = col, level = level, filename = filename)

    err = fixXMLError(err)
    errors[[length(errors) + 1]] <<- err

  }

  structure(list(handler = h, errors = function() structure(errors, class = "XMLStructuredErrorList"), reset = function() errors <<- list),
              class = "XMLStructuredErrorCumulator")
}

setOldClass("XMLStructuredErrorList")

print.XMLStructuredErrorList =
function(x, ...) {
   if(length(x) == 0)
     print(NULL)
   else
     print(t(sapply(x, function(x) unlist(x[c("line", "msg")]))))
}

getXMLErrors=
  #
  #  This attempts to read the specified file using the function given in parse
  # and then returns a list of the errors in the document.
  # This a somewhat convenient mechanism for fixing up, e.g., malformed HTML 
  # pages or other XML documents.
  
function(filename, parse = xmlParse, ...)
{
  f = xmlErrorFun()
  opts = options()
  options(error = NULL)
  on.exit(options(opts))
  tryCatch(parse(filename, ..., error = f$handler), error = function(e){})
  f$errors()                
}      



# Low level error handler
setXMLErrorHandler =
function(fun)
{
  prev = .Call("RS_XML_getStructuredErrorHandler", PACKAGE = "XML")

  sym = getNativeSymbolInfo("R_xmlStructuredErrorHandler", "XML")$address

  .Call("RS_XML_setStructuredErrorHandler", list(fun, sym), PACKAGE = "XML")
  
  prev
}


fixXMLError =
function(err)
{
  err$domain = getEnumValue(err$domain, xmlErrorDomain)
  err$code = getEnumValue(err$code, xmlParserErrors)
  err$level = getEnumValue(err$level, xmlErrorLevel)

  class(err) = "XMLError"
  
  err
}


getEnumValue =
function(value, defs)
{
    # might use for the class.
  name = substitute(defs)

  i = which(value == defs)

  defs[i]
}
