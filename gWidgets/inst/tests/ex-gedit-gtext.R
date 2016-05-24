w <- gwindow("edit, text", visible=FALSE)
g <- ggroup(cont = w, horizontal=FALSE)

## gedit

e <- gedit("edit", cont = g)

# svalue
print(svalue(e))

# svalue<-
svalue(e) <- "new text"

# handler
addHandlerChanged(e, function(h,...) print("changed"))
addHandlerKeystroke(e, function(h,...) print(h$key))
##?? others

## gtext
t <- gtext("edit", cont = g)

#svalue
print(svalue(t))

# svalue<-
svalue(t) <- "new text"

# insert
insert(t, "more new text") # ?\n?
insert(t, "even more with font.attr", font.attr = c("color"="red"))

# handler
addHandlerKeystroke(t, handler=function(h,...) print(h$key))



visible(w) <- TRUE
