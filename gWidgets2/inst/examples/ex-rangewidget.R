if(interactive()) {
  ## a range widget uses either a slider or a linked spinbutton to select a value
  w <- gwindow("Range widget", visible=FALSE)
  g <- ggroup(cont=w, horizontal=TRUE)
  sl <- gslider(from=0, to=100, by=1, value=0, cont=g, expand=TRUE, fill="both")
  sp <- gspinbutton(from=0, to=100, by=1, value=0, cont=g)

  ## Two ways to do this:
  ##  addHandlerChanged(sl, function(...) svalue(sp) <- svalue(sl))
  ##  addHandlerChanged(sp, function(...) svalue(sl) <- svalue(sp))

  f <- function(h, ...) svalue(h$action) <- svalue(h$obj)
  addHandlerChanged(sl, f, action=sp)
  addHandlerChanged(sp, f, action=sl)
  
  visible(w) <- TRUE
}
