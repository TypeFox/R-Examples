library(XML)
gctorture()
replicate(4, {
b = xmlParse('data/book.xml')
node = getNodeSet(b, "//chapter[2]/section[1]")[[1]]

p = newXMLNode("para", "This is ", newXMLNode("em", "emphasized text"), " within a paragraph", parent = node)
newXMLNode("para", "This is another paragraph", parent = node)
newXMLTextNode("This is second sentence", parent = p)
rm(b)
rm(p)
})

# above works with no problem.




