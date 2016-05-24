w <- gwindow("slider, spinbutton", visible=FALSE)
g <- ggroup(cont = w, horiz = FALSE)


## slider
sl <- gslider(from=0, to = 100, by =1, value=50, cont = g)

## svalue
print(svalue(sl))

## svalue <-
svalue(sl) <- 75

## handler
addHandlerChanged(sl, function(h,...) print(svalue(h$obj)))

## spinbutton
sb <- gspinbutton(from=0, to = 100, by =1, value=50, cont = g)

## svalue
print(svalue(sb))

## svalue <-
svalue(sb) <- 75

## handler
addHandlerChanged(sb, function(h,...) print(svalue(h$obj)))


visible(w) <- TRUE
