library(XML)

getRefCount =
function(node)
   .Call("R_getXMLRefCount", node)

doc = xmlParse(system.file("exampleData", "mtcars.xml", package = "XML"))

els = getNodeSet(doc, "//record")
sapply(els, getRefCount)
rm(els)
gc()

els = getNodeSet(doc, "//record", addFinalizer = FALSE)
sapply(els, getRefCount)


r = xmlRoot(doc, addFinalizer = FALSE)
getRefCount(r)
r = xmlRoot(doc, addFinalizer = NA)
getRefCount(r)
rm(r); gc()

r = xmlRoot(doc, addFinalizer = FALSE)
kids = xmlChildren(r, addFinalizer = FALSE)
sapply(kids, getRefCount)
rm(kids)

r = xmlRoot(doc, addFinalizer = FALSE)
kids = xmlChildren(r, addFinalizer = FALSE)
sapply(kids, getRefCount)

a = kids[[1]]
b = getSibling(a, addFinalizer = FALSE)
getRefCount(b)

xmlSize(b)




