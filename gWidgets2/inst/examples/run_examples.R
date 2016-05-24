
exs <- list.files(system.file("examples", package="gWidgets2"), full=FALSE)
exs <- Filter(function(i) !grepl("run_examples", i), exs)

exs_full <- list.files(system.file("examples", package="gWidgets2"), full=TRUE)
exs_full <- Filter(function(i) !grepl("run_examples", i), exs_full)

w <- gwindow("Run examples", visible=FALSE)
g <- gvbox(cont=w, use.scrollwindow=TRUE)
tbl <- glayout(cont=g)

f <- function(nm, fname) {
  force(nm); force(fname)
  n <- dim(tbl)[1] + 1
  tbl[n, 1] <- gbutton("See source...", cont=tbl, handler=function(h,...) {
    w1 <- gwindow(sprintf("Source of %s", nm))
    size(w1) <- c(600, 500)
    g <- gvbox(cont=w1); g$set_borderwidth(10)
    x <- readLines(fname)
    txt <- gtext(paste(x, collapse="\n"), wrap=FALSE,
                 font.attr=list(style="monospace"),
                 cont=g, expand=TRUE, fill=TRUE)
    gseparator(cont=g)
    bg <- ggroup(cont=g); addSpring(bg)
    gbutton("dismiss", cont=bg, handler=function(h, ...) {
      dispose(w1)
    })
  })

  tbl[n, 2] <- gbutton(sprintf("Run %s example", nm), cont=tbl, handler=function(h,...) {
    source(fname, local=FALSE)
  })
}

mapply(f, exs, exs_full)
addSpring(g)                            # push to top
visible(w) <- TRUE

