w <- gwindow("gbutton", visible=FALSE)
g <- ggroup(horizontal = FALSE, cont = w)

b <- gbutton("button", cont = g)

# svalue
print(svalue(b))

#svalue<-
svalue(b) <- "new button text"

## addHandlerChanged (click)
addHandlerChanged(b, function(h,...) print("hi"))

# button with icon
b1 <- gbutton("quit", cont = g)

# action button
a <- gaction("action", icon = "quit", handler = function(h,...) print("action"))
gbutton(action=a, cont = g)

             
visible(w) <- TRUE
