if(interactive()) {
w <- gwindow("dialogs example", visible=FALSE)
g <- ggroup(cont = w, horizontal = FALSE)

## filebrowse
gfilebrowse("select a file", cont = g)

## calendar
gcalendar("date", cont = g)

## gfile
gbutton("gfile", cont = g, handler = function(...) {
  out <- gfile()
  print(out)
})

## gmessage
gbutton("gmessage", cont = g, handler = function(h,...) {
  gmessage("message", title="title", icon = "warning", parent = h$obj)
})

## gconfirm
gbutton("gconfirm", cont = g, handler = function(h,...) {
  out <- gconfirm("confirm", title="title", icon = "warning", parent = h$obj)
  print(out)
})

## ginput
gbutton("ginput", cont = g, handler = function(h,...) {
  out <- ginput("input", title="title", icon = "warning", parent = h$obj)
  print(out)
})


visible(w) <- TRUE

}
