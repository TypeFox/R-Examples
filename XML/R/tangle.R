#
# tangle code from an XML file to a collection of files
#

getXPathExpr =
function(language, nodeNames)
{
 paste(paste("//", unlist(outer(language, nodeNames, FUN = "paste", sep = ":")), sep = ""), collapse = "|")
}

getTargetFiles =
function(doc, language = names(xmlNamespaceDefinitions(doc)),
          nodeNames = c("code", "function", "plot", "class", "method"),
           xpath = getXPathExpr(language, nodeNames))
{
  if(is.character(doc))
    doc = xmlParse(doc)
  
  nodes = getNodeSet(doc, xpath)
  ans = structure(sapply(nodes, xmlGetAttr, "file"),
                    names = sapply(nodes, function(x) names(xmlNamespace(x)) ))

  ans = tapply(ans, names(ans), function(x) unique(unlist(x)))
  ans[ sapply(ans, length) != 0 ]
}

xmlTangle =
function(doc, files = getTargetFiles(doc, xpath = xpath), dir = ".",
         language = names(xmlNamespaceDefinitions(doc)),
         nodeNames = c("code", "function", "plot", "class", "method"),
           xpath = getXPathExpr(language, nodeNames))
{
  if(is.character(doc))
    doc = xmlParse(doc)

  if(length(files) == 0 && "r" %in% language) {
      return(tangleR(doc, out = NA))
  }
  
  files =
   structure(lapply(names(files),
                   function(ns) {
                     xp = paste("//", ns, ":", nodeNames, sep = "")
                     structure(sapply(files[[ns]],
                                      function(file) {
                                        expr = paste(xp, "[@file=", sQuote(file), "]", collapse = "|")
                                        paste(xpathSApply(doc, expr, xmlValue), collapse = "\n")
                                      }), names = files[[ns]])
                   }), names = names(files), class = "FileContentsList")

  if(!is.na(dir))
    save.FileContentsList (files, dir)
  else
    files
}

save.FileContentsList =
function(x, dir = ".")
{
  x = structure(unlist(x, recursive = FALSE), names = unlist(lapply(x, names)))
  files = paste(dir, names(x), sep = .Platform$file.sep)
  sapply(seq(along = files),
          function(i) cat(x[[i]], file = files[i]))
  files
}
  
