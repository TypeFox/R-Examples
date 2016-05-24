if(interactive()) {

  w <- gwindow("Buttons", visible=FALSE)
  g <- ggroup(cont=w, horizontal=FALSE)

  ## various buttons

  ## with icons
  b1 <- gbutton("open", cont=g)

  ## without icon
  b2 <- gbutton("ouvrir", cont=g)

  ## by an action
  act <- gaction("open", tooltip="open", icon="open", handler=function(...) {})
  b3 <- gbutton(action=act, cont=g)

  ## with a handler
  b4 <- gbutton("click me", cont=g, handler=function(h,...) {
    if(svalue(b2) == "open")
      svalue(b2) <- "ouvrir"
    else
      svalue(b2) <- "open"
  })


  ## handlers can be blocked/unblocked
  b5 <- gbutton("Click me for a message", cont=g)
  id <- addHandlerClicked(b5, function(h,...) print("Ouch"))
  b6 <- gcheckbox("toggle handler message", cont=g, use.togglebutton=TRUE, handler=function(h,...) {
      if (svalue(b6)) {
          blockHandler(b5, id)
      } else {
          unblockHandler(b5, id)
      }
  })
  
  visible(w) <- TRUE
}
