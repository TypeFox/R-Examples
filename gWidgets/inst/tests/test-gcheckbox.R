test.gcheckbox <- function() {
  w <- gwindow()
  g <- ggroup(cont = w, horiz = FALSE)

  text <- "text"; newText <- "newtext"

  
  l <- gcheckbox(text, checked=TRUE, cont = g)
  
  ## svalue
  checkEquals(svalue(l), TRUE)

  ## svalue<-
  svalue(l) <- FALSE
  checkEquals(svalue(l), FALSE)

  ## [
  checkEquals(l[], text)

  ## [<-
  l[] <- newText
  checkEquals(l[], newText)
}
