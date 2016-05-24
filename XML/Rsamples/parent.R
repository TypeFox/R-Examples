library(XML)
a = newXMLNode("foo")
newXMLNode("bar", parent = a)
p = newXMLNode('top')
xmlParent(a) = p
rm(a)
gc()
rm(p)
gc()
