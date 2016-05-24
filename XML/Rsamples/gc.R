library(XML)

gctorture(TRUE)
node = newXMLNode("top")
doc = newXMLDoc(node = node)
print(doc)
print(node)
attr(node, "document") = doc
rm(doc)
gc()



# Node sets.

top = newXMLNode("top")
nodes = lapply(c("A", "B", "C"), newXMLNode, parent = top)
lapply(c("x", "y", "z"), newXMLNode, parent = nodes[[1]])
lapply(c("x", "y", "z"), newXMLNode, parent = nodes[[2]])

cat(saveXML(top))

gctorture(TRUE)
z = xpathApply(top, "//y|//z")
rm(z)

a = xpathApply(top, "//A")
print(a)
print(nodes[[1]])

x = xpathApply(a[[1]], "//x")
print(x)
rm(a)
gc()
rm(x)
gc()


# From Martin.

xml <- xmlTreeParse("<A><B></B></A>", useInternal=TRUE, asText=TRUE)
nds <- xpathApply(xml, "/*/*")
xpathApply(nds[[1]], "/*/*", xmlName)
#[[1]]
#[1] "B"

#vs. in 1.93-2.2:

xml <- xmlTreeParse("<A><B></B></A>", useInternal=TRUE, asText=TRUE)
nds <- xpathApply(xml, "/*/*")
xpathApply(nds[[1]], "/*/*", xmlName)
#list()



