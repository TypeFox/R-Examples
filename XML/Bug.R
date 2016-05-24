library(XML)
z = newXMLNode("foo")
save(z, "/tmp/z.rda")
load("/tmp/z.rda")
z
