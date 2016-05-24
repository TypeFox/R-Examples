w <- gwindow("glabel example", visible=FALSE)

widget <- glabel("test label", cont = w)

# svalue
print(svalue(widget))

# svalue<-
svalue(widget) <- "new label"

# font
font(widget) <- c("weight"="bold", "color"="red")

visible(w) <- TRUE
