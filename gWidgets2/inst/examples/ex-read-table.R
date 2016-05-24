## read.csv

about <- "
This example shows a mock up for an interface to `read.table`. The basic
layout follows a similar dialog from RStudio.

As an example, we see the text widget and table widget being used to
present the pre and post reading of the data. The text widget is
configured with no wrapping and a monospace font. Otherwise, the
example is fairly standard.
"

f <- tempfile()
write.csv(mtcars, file=f)

## layout
w <- gwindow("Read csv", visible=FALSE)
pg <- ggroup(cont=w)
pg$set_borderwidth(10L)

lg <- gvbox(cont=pg)
rg <- gvbox(cont=pg, expand=TRUE, fill=TRUE)

## put main text widget and table widget into paned group so user can
## adjust allocated space
pg <- gpanedgroup(cont=rg, expand=TRUE, horizontal=FALSE)
fr <- gvbox(cont=pg, expand=TRUE, fill=TRUE, spacing=5)
l <- glabel(gettext("Input file:"), cont=fr, anchor=c(-1,0))
font(l) <- list(weight="bold")
raw_input <- gtext("", wrap=FALSE,
                   font.attr=list(family="monospace"),
                   container=fr, expand=TRUE, fill=TRUE)

fr <- gvbox(cont=pg, expand=TRUE, fill=TRUE, spacing=5)
l <- glabel(gettext("Data frame:"), cont=fr, anchor=c(-1,0))
font(l) <- list(weight="bold")
output <- gtable(data.frame(X1=""), expand=TRUE, fill=TRUE, cont=fr)    

## buttons in a button group
bg <- ggroup(cont=rg)
addSpring(bg)
import_btn <- gbutton("import", cont=bg)
cancel_btn <- gbutton("cancel", cont=bg, handler=function(...) dispose(w))
about_btn <- gbutton("about", cont=bg, handler=function(...) {
  w1 <- gwindow("About",  parent=w)
  g <- gvbox(cont=w1); g$set_borderwidth(10L)
  glabel(about, cont=g)
  gseparator(cont=g)
  bg <- ggroup(cont=g); addSpring(bg)
  gbutton("dismiss", cont=bg, handler=function(...) dispose(w1))
})
  
## use a form layout for ease in laying out the controls to adjust the
## arguments for `read.table`
flyt <- gformlayout(cont=lg)
addSpring(lg)

name <- gedit("", initial.msg="Variable name", cont=flyt, label="Name")

heading <- gradio(c("Yes", "No"), horizontal=TRUE, cont=flyt, label="Heading")

seps <- c("Whitespace"=" ", "Comma" = ",", Semicolon=";", "Tab"="\t")
gcombobox(names(seps), selected=2, cont=flyt, label="Separator")

decs <- c("Period"=".", "Comma"=",")
gcombobox(names(decs), cont=flyt, label="Decimal")

quotes <- c("Double quote (\")" = '"',
            "Single quote (')" = "'",
            "No quote" = "")
gcombobox(names(quotes), cont=flyt, label="Quote")


## interactions
update_output <- function(...) {
  vals <- svalue(flyt)
  l <- list(file=f,
            header=vals$Heading == "Yes",
            sep=seps[vals$Separator],
            quote=quotes[vals$Quote],
            dec=decs[vals$Decimal])
  out <- do.call(read.table, l)
  output[] <- out
  out
}
sapply(flyt$children, addHandlerChanged, update_output)


addHandlerChanged(import_btn, handler=function(h,...) {
  DF <- update_output()
  nm <- svalue(name)
  if(exists(nm))
    if(!gconfirm(c("Variable exists", "really overwrite?"), parent=w))
      return()
  assign(nm, DF, .GlobalEnv)
  galert(sprintf("Saved data to %s", nm), parent=w)
})


## initialize
svalue(raw_input) <-  paste(readLines(f), collapse="\n")
update_output()
visible(w) <- TRUE
## set size after, not before. Otherwise gWidgetstclk computes the
## size allocation for the button group incorrectly
size(w) <- c(800, 400)

