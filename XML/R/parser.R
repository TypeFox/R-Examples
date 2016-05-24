# cat(paste(names(parserOptions), paste(2, 0:(length(parserOptions) - 1), sep = "^"), sep = " = ", collapse = "\n"))

#setClass("XMLParserOption", "EnumValue")

parserOptions =
  structure(c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 
4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576
), .Names = c("RECOVER", "NOENT", "DTDLOAD", "DTDATTR", "DTDVALID", 
"NOERROR", "NOWARNING", "PEDANTIC", "NOBLANKS", "SAX1", "XINCLUDE", 
"NONET", "NODICT", "NSCLEAN", "NOCDATA", "NOXINCNODE", "COMPACT", 
"OLD10", "NOBASEFIX", "HUGE", "OLDSAX"))

RECOVER = 2^0
NOENT = 2^1
DTDLOAD = 2^2
DTDATTR = 2^3
DTDVALID = 2^4
NOERROR = 2^5
NOWARNING = 2^6
PEDANTIC = 2^7
NOBLANKS = 2^8
SAX1 = 2^9
XINCLUDE = 2^10
NONET = 2^11
NODICT = 2^12
NSCLEAN = 2^13
NOCDATA = 2^14
NOXINCNODE = 2^15
COMPACT = 2^16
OLD10 = 2^17
NOBASEFIX = 2^18
HUGE = 2^19
OLDSAX = 2^20

xmlParseDoc =
function(file, options = 1L, encoding = character(), asText = !file.exists(file),
          baseURL = file)
{
   if(is.character(options)) {
     i = pmatch(options, names(parserOptions))
     if(any(is.na(i)))
         stop("unrecognized XML parser options: ", paste(options[is.na(i)], collapse = ", "))

     options = parserOptions[i]
   } else {
     if(!all(options %in% parserOptions))
        stop("unrecognized XML parser options: ", paste(options[!(options %in% parserOptions)], collapse = ", "))
   }
   options = as.integer(sum(options))

   if(asText)
     .Call("R_xmlReadMemory", file, nchar(file), as.character(encoding), options, as.character(baseURL), PACKAGE = "XML")
   else
     .Call("R_xmlReadFile", path.expand(file), as.character(encoding), options, PACKAGE = "XML")
}
