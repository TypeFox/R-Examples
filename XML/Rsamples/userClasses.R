library(XML)

fileName<- system.file("exampleData", "test.xml", package="XML")

fdoc<- xmlTreeParse(fileName, useInternalNodes = TRUE)
new<-structure(fdoc, class=c("new", class(fdoc)))

print.new <- function(x,..){
  cat("Say we add data specific info\n")
  NextMethod(XML:::print.XMLInternalDocument, x)
}
setOldClass(c("new",  "XMLInternalDocument", "XMLAbstractDocument"))

setAs("new", "character", function(from) saveXML(from))
new
