## Test whether a main window works with
## * menubar
## * toolbar
## * content area
## * statusbar


## toolbar style
style="both"


f <- function(...) print("hi")

a <- gaction(label="action", icon = "quit", handler = f)

mbl <- list(File=list(
              save=gaction("save", icon="save",handler=f),
              test=a)
            )
tbl <- list(file=gaction("save",icon="save",handler=f),
            test = a,
            stop = gaction("stop", icon="stop",handler=function(...) dispose(w)),
            quit1 = gaction("quit", icon="quit",handler=f)
            )

w <- gwindow("test window", visible=FALSE)
mb <- gmenu(mbl, cont=w)
tb <- gtoolbar(tbl, cont=w, style=style)
## main content
txt <- gtext(cont=w)
sb <- gstatusbar("status", cont=w)


## test statusbar
svalue(sb) <- "This was added to status bar"

## test window
svalue(w) <- "title added via svalue"


visible(w) <- TRUE


## tests for RUnit
test.gwindow <- function() {
  title <- "test"
  w <- gwindow(title)
  checkEquals(svalue(w), title)
  svalue(w) <- toupper(title)
  checkEquals(svalue(w), toupper(title))
  dispose(w)
  checkException(svalue(w) <- title)
}
