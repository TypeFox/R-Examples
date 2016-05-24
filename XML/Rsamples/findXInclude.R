library(XML)
doc = xmlParse("~/Books/XMLTechnologies/book.xml")
rc = getNodeSet(doc, "//r:code | //r:function", "r")

# ParseXML.xml line 923
findXInclude(rc[[161]])
getNodeLocation(rc[[161]])




# But this should be in ParseXML.xml at line 923.
# So the line number is correct, but the file name is incorrect.



# gdb break RS_XML_xmlNodeName
xmlName(rc[[161]])
# print *node->parent->parent->parent->parent->prev

# rc[[134]] is in convenienceFcns.xml

findXInclude(rc[[135]], recursive = TRUE)



doc = xmlParse("~/Books/XMLTechnologies/Rpackages/XDocTools/inst/sampleDocs/simpleSections.xml")
ss = getNodeSet(doc, "/*/title | //section")
lapply(ss, getNodeLocation)


#

myXInc =
  function(node)
  {
      x = node
        while(!is.null(x)) {
               prev = getSibling(x, FALSE)
                    if(inherits(prev, "XMLXIncludeStartNode"))
                              return(xmlAttrs(prev))
                    x = xmlParent(x)
             }

        NULL
    }


