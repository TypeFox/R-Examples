  startElement =
  function(ctxt, name, attrs, ...)  {
      print(ctxt)
      cat("<startElement>", name, attrs, "\n")
      if(name == "rewriteURI") {
           cat("Terminating parser\n")
	   xmlStopParser(ctxt)
      }
  }
  class(startElement) = "XMLParserContextFunction"  
  endElement =
  function(name, ...) 
    cat("<endElement>", name, "\n")

  fileName = system.file("exampleData", "catalog.xml", package = "XML")
  xmlEventParse(fileName, handlers = list(startElement = startElement, endElement = endElement))
