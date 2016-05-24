test.glabel <- function() {
  w <- gwindow()
  g <- ggroup(cont = w, horiz = FALSE)

  text <- "label text"; newText <- "new"
  l <- glabel(text, cont = g)

  # svalue
  checkEquals(svalue(l), text)

  # svalue<-
  svalue(l) <- newText
  checkEquals(svalue(l), newText)

  # font<-
  font(l) <- c("weight"="bold", color="red")

  # editable
  l1 <- glabel(text, editable=TRUE, cont = g)

}
