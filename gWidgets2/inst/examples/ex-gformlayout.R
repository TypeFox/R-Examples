## Simple example of gformlayout
w <- gwindow("Form layout example")

g <- ggroup(cont=w, horizontal=FALSE)
fl <- gformlayout(cont=g)

## arguments for a t.test
gedit("", initial="A variable name", label="x", cont=fl)
e <- gedit("", initial="(optional) variable name", label="y", cont=fl)
gcombobox(c("two.sided", "less", "greater"), label="alternative", cont=fl)
gedit("", initial="mean", label="mu", coerce.with=as.numeric, cont=fl)
cb <- gcheckbox("", checked=FALSE, label="paired", cont=fl)

## link up components
enabled(cb) <- FALSE
addHandlerChanged(e, handler=function(h, ...) enabled(cb) <- nchar(svalue(h$obj)) > 0)

gbutton("run t.test", cont=g, handler=function(h,...) {
  values <- svalue(fl)
  values <- try(check_values(values), silent=TRUE)
  if(!inherits(values, "try-error")) {
    out <- capture.output(do.call("t.test", values))
    cat(paste(out, collapse="\n"))
  } else {
    galert("Error with input values", parent=w)
  }
})


## check that values are correct. Throw error if not
check_values <- function(values) {
  values$x <- get(values$x, .GlobalEnv)  #
  if(nchar(values$y) == 0)
     values$y <- NULL
  if(!is.null(values$y))
     values$y <- get(values$y, .GlobalEnv)
  if(is.na(values$mu))
    stop("Mean is not specified")

  values
}
