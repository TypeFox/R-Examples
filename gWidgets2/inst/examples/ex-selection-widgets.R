about <- "
The selection widgets in `gWidgets`
-----------------------------------

The selection widgets are: `gcheckbox`, `gradio`, `gcheckboxgroup`,
`gcombobox`, and `gtable`.

These have similar features:

* The values to select from are passed to the `items` argument and are accessible
  through `[` and `[<-`.
* The selected value is available from `svalue`.
* The selected index is found by `svalue(obj, index=TRUE)`.
* The `drop` argument is sometimes useful (`gtable`).
* The `names` method refers to column headers or for `gcheckbox` the label
* The `length` method gives the length of the items to select from
"

library(gWidgets2)
w <- gwindow("Selection widgets", visible=FALSE)
gtoolbar(list(
              gaction("About", icon="help", handler=function(...) {
                w1 <- gwindow("About this example", visible=FALSE)
                g <- ggroup(cont=w1) 
                l <- glabel(about, cont=g); font(l) <- list(style="monospace")
                visible(w1) <- TRUE
              },parent=w),
              gseparator(parent=w),
              gaction("Quit", icon="cancel", handler=function(...) dispose(w), parent=w)              
              ),
         cont=w)

nb <- gnotebook(cont=w)

## gcheckbox
g <- gvbox(label="checkbox", cont=nb)
cb <- gcheckbox("checkbox", cont=g)
cb1 <- gcheckbox("toggle", use.toggle=TRUE, cont=g)
addSpring(g)

## gcheckboxgroup
g <- gvbox(label="checkboxgroup", cont=nb)
cbg <- gcheckboxgroup(c("check", "box", "group"), cont=g)
cbg <- gcheckboxgroup(c("check", "box", "group", "using", "table"), use.table=TRUE, cont=g)
addSpring(g)

## gradio
g <- gvbox(label="gradio", cont=nb)
rb <- gradio(c("horizontal", "radio"), horizontal=TRUE, cont=g)
rb <- gradio(c("vertical", "radio"), vertical=TRUE, cont=g)
addSpring(g)

## gcombobox
g <- gvbox(label="gcombobox", cont=nb)
x <- data.frame(a = state.name[1:5],
                b = c("ok", "help", "quit", "plot", "cancel"),
                c = state.abb[1:5])

## a vector
combo1 <- gcombobox(x$a, cont=g)
## a data frame
combo2 <- gcombobox(x['a'], cont=g)
## with icon (possibly) and tooltip (possibly)
combo2 <- gcombobox(x, cont=g)
addSpring(g)

## gtable
g <- gnotebook(label="gtable", cont=nb)
## wotj a vector
gtable(x$a, cont=g, label="from vector")
## with an icon
gtable(x[1:2], cont=g, icon.col=2, label="with icon")
## with icon and tooltip
gtable(x, cont=g, icon.col=2, tooltip.col=3, label="icon/tooltip")

svalue(nb) <- 1
visible(w) <- TRUE

