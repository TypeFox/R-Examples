library(XML)

SVG.ns = c("http://www.w3.org/2000/svg") 
svg = newXMLNode("svg", namespaceDefinitions = SVG.ns)

nn = newXMLNode("omg:x", namespaceDefinitions = c(omg = "http://www.omegahat.net"))
nn = newXMLNode("omg:x", namespaceDefinitions = c(omg = "http://www.omegahat.net",
                                                   r = 'http://www.r-project.org'))

svg = newXMLNode("svg", namespaceDefinitions =  c(omg = "http://www.omegahat.net"))

svg = newXMLNode("omg:svg", namespaceDefinitions =  c(omg = "http://www.omegahat.net"))
svg = newXMLNode("no:svg", namespaceDefinitions =  c(omg = "http://www.omegahat.net"))


SVG.ns = c("http://www.w3.org/2000/svg") 
svg = newXMLNode("top", namespaceDefinitions =  SVG.ns)
g = newXMLNode("g", parent = svg)
r = newXMLNode("r:code", namespaceDefinitions = c(r = "http://www.r-project.org"), parent = g)
k = newXMLNode("output",  parent = r)
xmlNamespace(k)


XML:::xmlNamespaceRef(svg)




doc = newXMLDoc()
SVG.ns = c("http://www.w3.org/2000/svg") 
svg = newXMLNode("top", namespaceDefinitions =  SVG.ns, doc = doc)
g = newXMLNode("g", parent = svg)
r = newXMLNode("r:code", namespaceDefinitions = c(r = "http://www.r-project.org"), parent = g)
k = newXMLNode("output",  parent = r)
xmlNamespace(k)


getNodeSet(doc, "//namespace::*")
