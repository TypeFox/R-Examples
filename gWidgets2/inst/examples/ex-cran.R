about <- "
This example shows a simple interface to select a CRAN mirror.
It extends the `chooseCRANmirror` functionality by adding
filtering and giving a reasonable default size.
"

## our data

mirrors <- getCRANmirrors()

## Basic layout
w <- gwindow("Select a CRAN mirror", visible=FALSE)
g <- gvbox(container=w)
g$set_borderwidth(10)                   # svalue is between widget spacing
g1 <- ggroup(container=g)
glabel("Filter:", container=g1)
filter_by <- gedit("", initial.msg="regular expression", container=g1)

## main widget
tbl <- gtable(mirrors['Host'], container=g, expand=TRUE)

## button group
g2 <- ggroup(container=g)
addSpring(g2)
gbutton("About", container=g2, handler=function(h,...) {
  w1 <- gwindow("About", parent=w)
  box <- gvbox(container=w1, spacing=10)
  box$set_borderwidth(10)
  glabel(about, container=box)
})

select <- gbutton("Select", container=g2)
enabled(select) <- FALSE

## Set size and show ..
size(w) <- c(300, 400)
visible(w) <- TRUE


## Put in interactivity

addHandlerKeystroke(filter_by, function(h,...) {
  re <- svalue(h$obj)
  if(nchar(re) == 0) {
    visible(tbl) <- TRUE
  } else {
    visible(tbl) <- grepl(re, mirrors$Host)
    if(sum(visible(tbl)) == 1) { ## unique item
      svalue(tbl, index=TRUE) <- which(visible(tbl))
    }
  }
})


addHandlerSelectionChanged(tbl, handler=function(h,...) {
  enabled(select) <- length(svalue(h$obj)) > 0
})


mirror_selected <- function(...) {
  ind <- svalue(tbl, index=TRUE)
  if(length(ind) == 0) {
    galert("No mirror selected", parent=w)
  } else {
    message("Set CRAN to ", mirrors[ind, 1])
    dispose(w)
  }
}

sapply(list(select, tbl), "addHandlerChanged", handler=mirror_selected)

