w <- gwindow("gtext example", visible=FALSE)

## constrouctor -- font.attr sets for buffer
widget <- gtext("test text", cont = w, font.attr=list(size=24L, color="blue"))
insert(widget, "new text", font.attr=c(family="monospace"))

                
# svalue
print(svalue(widget))

# svalue<-
svalue(widget) <- "new label"

# font<-
# sets for buffer if no selection
#font(widget) <- list(family="monospace", "weight"="bold", "color"="red", size="xx-large")

visible(w) <- TRUE
