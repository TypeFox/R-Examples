# DTL - Dec 14, 2011
libary(XML)
doc = xmlParse("tracing.xml")

nsURI = c(n = "http://www.mbfbioscience.com/2007/neurolucida")

contours = getNodeSet(doc, "//n:contour", nsURI)
length(contours) # 2

getXY =
function(node)
{
  as(xmlAttrs(node)[c("x", "y")], "numeric")
}

c.xy = sapply(contours, function(x) t(sapply(x[names(x) == "point"], getXY)))

markers = getNodeSet(doc, "//n:marker/n:point", nsURI)
markers.xy = as.data.frame(t(sapply(markers, getXY)))
names(markers.xy) = c("x", "y")
markers.xy$Site = factor(xpathSApply(doc, "//n:marker/n:property", xmlValue, namespaces = nsURI))

barplot(table(markers.xy$Site))

plot(c.xy[[1]], type = "l", xlab = "x", ylab = "y")
points(markers.xy, col = "red")
