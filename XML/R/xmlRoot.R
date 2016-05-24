setGeneric("xmlRoot<-",
             function(x, ..., value)
                standardGeneric("xmlRoot<-"))

setMethod("xmlRoot<-", c("XMLInternalDocument", value = "character"),
function(x, ..., value)
{
   newXMLNode(value, doc = x)
   x
})

setMethod("xmlRoot<-", c("XMLInternalDocument", value = "XMLInternalNode"),
function(x, ..., value)
{
  #XXX check that this does the reference counting correctly
  # specifically, d = newXMLDoc(); xmlRoot(d) = "bar"; xmlRoot(d) = newXMLNode("foo")

  .Call("RS_XML_setRootNode", x, value, PACKAGE = "XML")
  x
})

setMethod("xmlRoot<-", "XMLHashTree",
function(x, ..., value)
{
  x$.addNode(value)
  x
})

