library(qtbase)


# Tests dynamic slots. Crashed R due to type-handlers.cpp incorrectly releasing memory.
a <- Qt$QLineEdit()
b <- Qt$QLineEdit()
qconnect(a, "textChanged", function(x) { b$text <- x })

b$text <- "2"
a$text <- "1"
if(!identical(a$text, b$text))
{
	stop("Encountered error in testing dynamic slots, ", a$text, " != ", b$text)
}


# Tests multiple inheritance. Crashed R due to environment holding wrong pointer.
gscene <- Qt$QGraphicsScene()
rtxt <- gscene$addText("some text")
rtxt
(rtxt$toPlainText()) ## fine

item.list <- gscene$items()
item.list[[1]]
(item.list[[1]]$toPlainText()) ## crash

