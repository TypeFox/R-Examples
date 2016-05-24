w <- gwindow("buttons")
g <- ggroup(cont=w, horizontal=FALSE)

## vanilla
b1 <- gbutton("label", cont=g)

## with icon
b2 <- gbutton("ok", cont=g)

## with handler
b3 <- gbutton("click me", cont=g, handler=function(h,...) message("ouch"))

## from action
a <- gaction("action", icon="help", tooltip="tooltip", handler=function(h,...) message("action"), parent=w)
b4 <- gbutton(action=a, cont=g)

## add handler
addHandlerChanged(b2, handler=function(h,...) print("add handler"))

## tests

## svalue
expect_equal(svalue(b1), "label")

## svalue <-
svalue(b1) <- "new label"
expect_equal(svalue(b1), "new label")

## invoke handler
expect_output(b2$invoke_change_handler(), "add handler")
