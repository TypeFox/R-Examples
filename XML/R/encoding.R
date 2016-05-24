CE_NATIVE = 0L
CE_UTF8 = 1L
CE_LATIN1 = 2L

  # Map an encoding or document's encoding to the corresponding R internal enum value
setGeneric("getEncodingREnum", 
              function(doc, ...)
	        standardGeneric("getEncodingREnum"))

setMethod("getEncodingREnum", "XMLInternalDocument",
           function(doc, ...) 
              getEncodingREnum( getEncoding(doc) ))

setMethod("getEncodingREnum", "XMLInternalElementNode", # was XMLInternalElement, but no such class?
           function(doc, ...) 
              getEncodingREnum( as(doc, "XMLInternalDocument") ))

setMethod("getEncodingREnum", "character",
           function(doc, ...) {
             if(length(doc) == 0 || is.na(doc))
               return(CE_NATIVE)
             
             str = tolower(doc)
             if(any(str == c("utf8", "utf-8")))
                 CE_UTF8
             else if(any(str == c("latin1", "iso-8859-1")))
                 CE_LATIN1
             else
                 CE_NATIVE # or NA?
           })
