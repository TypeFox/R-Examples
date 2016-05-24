library(RXMCDA)

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.0.0"), parent=tree)

root<-getNodeSet(tree, "/xmcda:XMCDA")

parameters<-newXMLNode("methodParameters", parent=root[[1]], namespace=c())

# value: integer=3
parameter <- newXMLNode("parameter",attrs = c(name="numIt"), parent=parameters, namespace=c())

value <- newXMLNode("value", parent=parameter, namespace=c())

newXMLNode("integer", value=3, parent=value, namespace=c())

# value: boolean="true"
parameter <- newXMLNode("parameter",attrs = c(name="true"), parent=parameters, namespace=c())

value <- newXMLNode("value", parent=parameter, namespace=c())

newXMLNode("boolean", value="true", parent=value, namespace=c())

# value: boolean="0"
parameter <- newXMLNode("parameter",attrs = c(name="falseAsInt"), parent=parameters, namespace=c())

value <- newXMLNode("value", parent=parameter, namespace=c())

newXMLNode("boolean", value=0, parent=value, namespace=c())

# value: invalid boolean
parameter <- newXMLNode("parameter",attrs = c(name="invalidBoolean"), parent=parameters, namespace=c())

value <- newXMLNode("value", parent=parameter, namespace=c())

newXMLNode("boolean", value="3", parent=value, namespace=c())


y<-getNodeSet(tree,"//methodParameters")

stopifnot(getParameters(y[[1]])$numIt == 3)

stopifnot(getParameters(y[[1]])$true == TRUE)
stopifnot(getParameters(y[[1]])$falseAsInt == FALSE)
stopifnot(identical(getParameters(y[[1]])$invalidBoolean,NA))


