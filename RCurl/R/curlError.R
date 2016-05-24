curlError =
function(type, msg, asError = TRUE)
{
   if(!is.character(type)) {
     i = match(type, CURLcodeValues)
     typeName = if(is.na(i))
                   character()
                else
                   names(CURLcodeValues)[i]
   } 

   typeName = gsub("^CURLE_", "", typeName)

   fun = (if(asError) stop else warning)

   fun( structure(list(message = msg, call = sys.call()), class = c(typeName, "GenericCurlError", "error", "condition")) )
}


getCurlErrorClassNames =
function()
{
  gsub("^CURLE_", "", names(CURLcodeValues))
}
