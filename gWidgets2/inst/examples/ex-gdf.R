library(gWidgets2)

## An example of using some controls to modify which rows are displayed using visible<-

X <- mtcars

## Layout
w <- gwindow()
g <- gvbox(cont=w)

d <- gdf(X, cont=g)
b <- ggroup(cont=g)

nms <- names(X)
combo <- gcombobox(nms, cont=b)
glabel("==", cont=b)                    # could easily generalize
val <- gedit("", cont=b)

## Handlers.
## if entry widget is non-empty, find matches to display within selected variable
addHandlerChanged(val, handler=function(h,...) {
    if (svalue(combo, index=TRUE) < 1)
        return()
    value <- svalue(val)    
    if(nchar(value) == 0)
        return()

    var <- X[[svalue(combo)]]
    
    ind <- var == value
    visible(d) <- ind
})

## if a new variable is changed, update the dropdown list and
## reset value
addHandlerChanged(combo, handler=function(h,...) {
    vars <- sort(unique(X[[svalue(combo)]]))
    val[] <- vars
    svalue(val) == ""
})
