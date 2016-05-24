getXIncludes =
function(filename, recursive = TRUE, skip = character(),
         omitPattern = "\\.(js|html?|txt|R|c)$",
          namespace = c(xi = "http://www.w3.org/2003/XInclude"),
          duplicated = TRUE)
{
   doc = xmlParse(filename, xinclude = FALSE)

   if(missing(namespace)) {
     ns = xmlNamespaceDefinitions(doc, simplify = TRUE)
     if("http://www.w3.org/2001/XInclude" %in% ns)
       namespace = c(xi = "http://www.w3.org/2001/XInclude")
   }

   nodes = getNodeSet(doc, "//xi:include", namespaces = namespace)
   files = sapply(nodes, xmlGetAttr, "href")
   nonRecursive = as.logical(sapply(nodes, xmlGetAttr, "text", FALSE))

   if(length(omitPattern))
      nonRecursive = grepl(omitPattern, files) | nonRecursive

   if(recursive) {
     processed = c(filename, skip)
     for(f in unique(files[!nonRecursive])) {
           # path name relative to the base document of the XInclude
        f = getRelativeURL(f, filename) # dirname(filename))
        if(file.exists(f))
           files = c(files, getXIncludes(f, TRUE, skip = processed))
        else
           warning(f, " doesn't exist")
        processed = c(processed, f)
     }
   }
   files = unlist(files)
   if(!duplicated)
     unique(files)
   else
     files
}
