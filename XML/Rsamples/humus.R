library(XML)
u = "http://humus101.com/?p=2737"
doc = htmlParse(u, encoding = "UTF-8")
v = getNodeSet(doc, "//strong")[[1]]
getEncoding(v)
# [1] "UTF-8"
xmlValue(v, encoding = "UTF-8")
Encoding(xmlValue(v, encoding = "UTF-8"))

