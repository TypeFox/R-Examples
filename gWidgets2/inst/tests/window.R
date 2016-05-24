w <- gwindow("windows and dialogs")
sb <- gstatusbar("hello", cont=w)
g <- ggroup(cont=w, horizontal=FALSE)


## handlers

## Called when window is closed
addHandlerDestroy(w, function(h,...) print("Closing ..."))

## called when close. Return TRUE if you want the window to close.
addHandlerUnrealize(w, function(h,...) gconfirm("Really close", parent=h$obj))

## test

## svalue
expect_equal(svalue(w), "windows and dialogs")
## svalue<-
svalue(w) <- "new title"
expect_equal(svalue(w), "new title")


## position window by geometry (x+y pixels from upper left)
gwindow("100, 200 down", parent=c(100, 200))


w2 <- gwindow("parent")
w3 <- gwindow("child", parent=w2)
dispose(w2)
Sys.sleep(1)
expect_equal(isExtant(w3), FALSE)




## subwindow
gbutton("subwindow", cont=g, handler=function(h, ...) {
  w1 <- gwindow("Subwindow", parent=w)
  g <- ggroup(cont=w1, horizontal=FALSE)
  txt <- gtext("", cont=g, expand=TRUE)
  gbutton("dismiss", cont=g, handler=function(h,...) dispose(w1))
})


## dialogs

## alert
gbutton("alert", cont=g, handler=function(h,...) galert('alert', parent=w))


## messasge
gbutton("message", cont=g, handler=function(h, ...) {
  gmessage("message", parent=w)
})

## confirm
gbutton("confirm", cont=g, handler=function(h,...) {
  out <- gconfirm("Okay", parent=w)
  print(out)
})

## input
gbutton("input", cont=g, handler=function(h,...) {
  out <- ginput(msg="message", text="initial", title="title", parent=w)
  print(out)
})





