# This writes an XML file using the specified top level node
writeXML <- function(node, file)
{
   if (is.character(file))
   {
      file <- file(file, open='wb')
      on.exit(close(file))
   }

   writeRaw <- function(s, file) writeBin(charToRaw(s), file)

   writeXML.recurse <- function(node)
   {
      if (inherits(node, 'XMLTextNode'))
      {
         if (!is.null(node$raw))
         {
            writeRaw(xmlValue(node), file=file)
         } else {
            writeRaw(escape(xmlValue(node)), file=file)
         }
      } else {
         writeRaw(sprintf('<%s', xmlName(node, full=TRUE)), file=file)
         atts <- xmlAttrs(node)

         # This won't be necessary for "minixml" package
         ns <- xmlNamespaceDefinitions(node)
         if (! is.null(ns) && length(ns) > 0)
         {
            xatts <- sapply(ns, function(n) n$uri)
            nms <- names(xatts)
            if (! is.null(nms) && length(nms) > 0) {
               names(xatts) <- paste('xmlns:', nms, sep='')
               atts <- c(xatts, atts)
            }
         }

         writeRaw(genXMLAttributes(atts), file=file)

         if (xmlSize(node) == 0)
         {
            writeRaw('/>', file=file)
         } else {
            writeRaw('>', file=file)
            for (child in xmlChildren(node))
            {
               writeXML.recurse(child)
            }
            writeRaw(sprintf('</%s>', xmlName(node, full=TRUE)), file=file)
         }
      }
   }

   writeRaw('<?xml version="1.0" encoding="UTF-8"?>\n', file=file)
   writeXML.recurse(node)
   writeRaw('\n', file=file)
}
