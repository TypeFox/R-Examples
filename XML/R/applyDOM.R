xmlDOMApply <- 
function(dom, func)
{
 .Call("RS_XML_RecursiveApply", dom, func, NULL, PACKAGE = "XML")
}
