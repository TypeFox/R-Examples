if(FALSE)  {
library(Rffi)
free.cif = CIF(voidType, list(pointerType))
freeDoc =
function(doc)
{
   callCIF(free.cif, "xmlFreeDoc", doc)
}
}

library(XML)
ff = list.files("~/Data/Kiva/loans", pattern = 'xml$', full.names = TRUE)
docs = lapply(ff, function(x) {
               print(x);
               doc = xmlParse(gsub("&#[0-9]+;", "", readLines(x, warn = FALSE)), asText = TRUE);
               amt = as.numeric(xpathSApply(doc, "//loan/funded_amount", xmlValue))
#               freeDoc(doc)
                gc()
               amt
             })
