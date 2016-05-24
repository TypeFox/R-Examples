test.gbutton <- function() {
  w <- gwindow()
  g <- ggroup(cont = w, horiz = FALSE)
  
  text <- "label text"; newText <- "quit"
  
  ## plain vanilla
  l <- gbutton(text, cont = g)
  
  ## svalue
  checkEquals(svalue(l), text)
  
  ## svalue<-
  svalue(l) <- newText
  checkEquals(svalue(l), newText)
  
  ## font<-
  font(l) <- c(weight="bold", color="red")

  # gaction
  a <- gaction(text, handler = function(h,...) print("hi"))
  l <- gbutton(action = a, cont =g)

  checkEquals(svalue(a), text)


  ## enabled
  b <- gbutton("asdf", cont=g)
  enabled(b) <- FALSE
  checkEquals(enabled(b), FALSE)
  

  
}
