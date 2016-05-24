test.gradio <- function() {
  w <- gwindow()
  g <- ggroup(cont = w, horiz = FALSE)

  items <- letters[1:4]
  
  l <- gradio(items, selected=1, cont = g)
  
  ## svalue
  checkEquals(svalue(l), items[1])
  checkEquals(svalue(l, index=TRUE), 1)

  ## svalue<-
  svalue(l) <- "b"
  checkEquals(svalue(l), "b")

  svalue(l, index=TRUE) <- 3
  checkEquals(svalue(l, index=TRUE), 3)

  ## [
  checkEquals(l[], items)

  ## [<-
  l[] <- items[1:4]
  checkEquals(l[], items[1:4])
}
