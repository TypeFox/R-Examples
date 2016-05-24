##' A means to generate a GUI automatically from a function's arguments
##'
##' In gWidgets the ggenericwidget function provided some of this. This is an alternative for the
##' super quick but super limited GUI.
library(gWidgets2)


## Set up an editor for each type of argument. More work here could be done
arg_editor <- function(x, cont, ...) UseMethod("arg_editor")
arg_editor.default <- function(x, cont, ...) gedit(x, cont=cont, ...)
arg_editor.call <- function(x, cont, ...) arg_editor(eval(x), cont, ...)
arg_editor.name <- function(x, cont, ...) {
  gedit("", cont=cont, ...)
}
arg_editor.numeric <- function(x, cont, ...) gcombobox(sort(unique(x)), editable=TRUE, coerce.with=as.numeric, cont=cont)
arg_editor.character <- function(x, cont, ...) gcombobox(sort(unique(x)), editable=TRUE, coerce.with=as.character, cont=cont)
arg_editor.factor <-  function(x, cont, ...) gcombobox(sort(unique(x)), cont=cont)
arg_editor.logical <- function(x, cont, ...) gcombobox(c(TRUE, FALSE), selected=2 - as.numeric(x[1]), cont=cont)

##' For each value from formals, this places an editor in to a table
make_gui <- function(FUN, parent) {
  args <- formals(FUN)
  if(missing(parent)) {
    parent <- gwindow("Some GUI", visible=FALSE)
  }

  g <- ggroup(cont=parent, horizontal=FALSE)

  lyt <- glayout(cont=g)
  n <- length(args)

  f <- function(i, nm, x) {
    lyt[i,1] <- nm
    lyt[i,2] <- arg_editor(args[[i]], cont=lyt)
  }
  mapply(f, i=seq_along(args), nm=names(args), x=args)

  bg <- ggroup(cont=g); addSpring(bg)
  gbutton("ok", cont=bg, handler=function(h, ...) {
    values <- lapply(lyt[,2], svalue)
    nms <- sapply(lyt[,1], svalue)
    names(values) <- nms
    do.call(FUN, values)
  })
  addSpace(bg, 12)
  gbutton("dismiss", cont=bg, handler=function(h,...) dispose(parent))

  visible(parent) <- TRUE

  invisible(lyt)
}

## Test it out.
s <- function(x=1:4, y = letters, z=TRUE) {
  print(list(x, y, z))
}

make_gui(s)



