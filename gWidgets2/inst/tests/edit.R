w <- gwindow("gedit", visible=FALSE)
g <- gformlayout(cont=w)

## vanilla
e1 <- gedit("vanilla", cont=g, label="vanilla")

## initial message
e2 <- gedit("", initial.msg="initial", cont=g, label="initial")

## handler
e5 <- gedit("", initial.msg="type a value", label="handler", cont=g, handler=function(h,...) {
  print("call handler")
})


## add a handler
e6 <- gedit("",cont=g, label="change handler" )
addHandlerChanged(e6,  handler=function(h,...) {
  print("call handler")
})

## blur handler
e7 <- gedit("",cont=g, label="blur handler")
addHandlerBlur(e7, handler=function(h,...) print("call handler"))

## keystroke
e8 <- gedit("", cont=g, label="keystroke handler")
addHandlerKeystroke(e8, handler=function(h,...) print(h$key))


visible(w) <- TRUE

## tests

## svalue
expect_equal(svalue(e1), "vanilla")

## svalue<-
svalue(e2) <- "added"
expect_equal(svalue(e2), "added")

## [<- typeahead. This is in combobox

## handler
expect_output(e6$invoke_change_handler(), "call handler")
