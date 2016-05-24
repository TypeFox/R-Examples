library(XML)
gctorture(TRUE)
b = newXMLNode("bob", namespace = c(r = "http://www.r-project.org", omg = "http://www.omegahat.net"))

cat(saveXML(b), "\n")

addAttributes(b, a = 1, b = "xyz", "r:version" = "2.4.1", "omg:len" = 3)
cat(saveXML(b), "\n")

removeAttributes(b, "a", "r:version")
cat(saveXML(b), "\n")

removeAttributes(b, .attrs = names(xmlAttrs(b)))

addChildren(b, newXMLNode("el", "Red", "Blue", "Green",
                           attrs = c(lang ="en")))

k = lapply(letters, newXMLNode)
addChildren(b, kids = k)

cat(saveXML(b), "\n")

removeChildren(b, "a", "b", "c", "z")

  # can mix numbers and names
removeChildren(b, 2, "e")  # d and e

cat(saveXML(b), "\n")

i = xmlChildren(b)[[5]]
xmlName(i)

 # have the identifiers
removeChildren(b, kids = c("m", "n", "q"))


x <- xmlNode("a", 
               xmlNode("b", "1"),
               xmlNode("c", "1"),
               "some basic text")

v = removeChildren(x, "b")

  # remove c and b
v = removeChildren(x, "c", "b")

  # remove the text and "c" leaving just b
v = removeChildren(x, 3, "c")

## Not run: 
##D     # this won't work as the 10 gets coerced to a 
##D     # character vector element to be combined with 'w'
##D     # and there is no node name 10.
##D  removeChildren(b, kids = c(10, "w"))
## End(Not run)

 # for R-level nodes (not internal)

z = xmlNode("arg", attrs = c(default="TRUE"),
              xmlNode("name", "foo"), xmlNode("defaultValue","1:10"))

o = addChildren(z,
                "some text",
                xmlNode("a", "a link", attrs = c(href = "http://www.omegahat.net/RSXML")))
o
