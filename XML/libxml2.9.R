library(XML)
f <- system.file("exampleData", "gnumeric.xml", package = "XML")
startElement = function(x, ...) cat(x, "\n")
xmlEventParse(file(f), handlers = list(startElement = startElement))

if(FALSE) {
 con = file(f, "r")
 g = function(...) { browser(); cat("here\n"); readLines(con, 1) }
 xmlEventParse(g, handlers = list(startElement = startElement))
 close(con)
}


