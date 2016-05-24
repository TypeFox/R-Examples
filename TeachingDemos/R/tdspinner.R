# these are utility functions, possibly will be replaced by true internals


tdspinner <- function(parent, ...) {
    # this is a quick hack to provide spinboxes without loading tcltk2
    tcltk::tkwidget(parent, "spinbox", ...)
}


have.ttk <- function() {
    # based on e-mail from Prof. Brian Ripley
    # will work until version 8.10 or 10.0, then may need to update
    as.character(tcltk::tcl("info","tclversion")) >= "8.5"
}
