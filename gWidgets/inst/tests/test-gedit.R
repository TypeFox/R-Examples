test.gedit <- function() {
  w <- gwindow()
  g <- ggroup(cont = w, horiz = FALSE)

  text <- "label text"; newText <- "new"
  l <- gedit(text, cont = g)

  # svalue
  checkEquals(svalue(l), text)

  # svalue<-
  svalue(l) <- newText
  checkEquals(svalue(l), newText)

  # [<-
  l[] <- state.name
  
}
