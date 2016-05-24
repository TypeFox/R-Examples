library(XML)

doc = xmlParse('<doc xmlns:r="http://www.r-project.org">
                 <el r:width="10" width="72"/>
                 <el width="46"/>
               </doc>')

a = xmlAttrs(xmlRoot(doc)[[1]])
a[1]

