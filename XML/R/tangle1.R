tangleR = xxx_getRCode =  # conflicts with getRCode in xmlInternalSource.R
function(doc, tags =  c("code", "plot", "function"), out = gsub("\\.[a-zA-Z]+$", ".R", docName(doc)))
{
   if(is.character(doc))
      doc = xmlParse(doc)

   xp = sprintf("//r:%s[not(@eval='false') and not(ancestor::section[@eval='false'])]",
                    tags)
   code = xpathSApply(doc, paste(xp, collapse = " | "), xmlValue, namespaces = c("r" = "http://www.r-project.org"))

   if(length(out) && !is.na(out)) {
      cat(code, sep = "\n", file = out)
      out
   } else
     code
}

