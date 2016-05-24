library(XML)
gctorture(TRUE)
a = newXMLNode("a") ; b = newXMLNode("b", "boo", parent = a); c = newXMLNode("c", "hoo", parent = a)
print(a)

#xmlValue(a[["b"]]) <- "mo"
#print(a)

xmlValue(a[["c"]]) <- "jo"
print(a)
