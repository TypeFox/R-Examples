## example of gaction
if(interactive()) {
  ## means to make a simple handler
  f <- function(x)  function(...) print(x)
  
  ## list of action objects
  lst <- list(new=gaction("new", icon="new", handler=f("new")),
              open=gaction("open", icon="open", handler=f("open")),
              save=gaction("save", icon="save", handler=f("save"))
              )
  
  ## can enable or disable an action -- all instances reflect the state. Svalue works too.
  enabled(lst$save) <- FALSE
  
  w <- gwindow("Example of actions in buttons and toolbar", visible=FALSE)
  tb <- gtoolbar(lst, cont=w)
  b <- gbutton(action=lst$open, cont=w)
  
  visible(w) <- TRUE
}
