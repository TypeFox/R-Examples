
## @knitr setup1, echo=FALSE
options(replace.assign=TRUE)
options(width=80)

opts_chunk$set(tidy=FALSE,
               highlight=TRUE,  
               #fig.width=5, fig.height=5,
               fig.show='hold', cache=FALSE, par=TRUE) 
knit_hooks$set(par=function(before, options, envir){
if (before && options$fig.show!='none') 
  par(mar=c(5,4,2,2) + 0.1, #c(4,4,.1,.1),
      cex.lab=.95,cex.axis=.9
      ,mgp=c(2,.7,0),tcl=-.3
                )
}, crop=hook_pdfcrop)



## @knitr setup, echo=FALSE, results="hide"
options(width=72)
##options(prompt='  ',continue='  ')  # remove prompt characters at start of lines


## @knitr hello_world, results="hide"
library(gWidgets2)
options(guiToolkit="RGtk2")
## containers
win <- gwindow("Basic example", visible=FALSE)
gp <- gvbox(container=win)
## control
btn <- gbutton("click me for a message", container=gp)
## interactivity
addHandlerClicked(btn, handler=function(h,...) {
  galert("Hello world!", parent = win)  # a dialog call
})
## a method call
visible(win) <- TRUE


## @knitr nested_container, results="hide"
## Some filler
lorem <- "Lorem ipsum dolor sit amet, consectetur adipiscing elit."
##
win <- gwindow("Nested groups")
g <- gvbox(container=win)
g$set_borderwidth(10L)
##
txt <- gtext(lorem, container=g,        # text widget
             expand=TRUE, fill=TRUE) 
##
bg <- ggroup(cont=g)
addSpring(bg)
gbutton("dismiss", container=bg, 
        handler=function(h,...) dispose(win))
gbutton("about", container=bg, handler=function(h,...) {
  gmessage("Shows lorem ipsum text", parent=win)
})


## @knitr gformlayout_example, results="hide"
win <- gwindow("t-test", visible=FALSE)
g <- gvbox(container=win)
g$set_borderwidth(10L)
##
flyt <- gformlayout(container=g, expand=TRUE)
## 
gedit("", initial.msg="variable", 
      label="x", container=flyt)
gcombobox(c("two.sided", "less", "greater"), 
          label="alternative", container=flyt)
gedit("0", coerce.with=as.numeric,  
      label="mu", container=flyt)
gcheckbox("", checked=FALSE, 
          label="paired", container=flyt)
gslider(from=0.5, to = 1.0, by=.01, value=0.95, 
        label="conf.level",  container=flyt)
##
bg <- ggroup(container=g)
addSpring(bg)
gbutton("values...", container=bg, handler=function(h,...) {
  print(svalue(flyt))                   # replace me...
})
addSpring(g)                            # better for Qt
##
size(win) <- c(400, 250)
visible(win) <- TRUE


## @knitr pick_your_race, results="hide"
win <- gwindow("handler example", visible=FALSE)
g <- gvbox(container=win)
f <- gframe("Ethnicity", container=g)
cb <- gcheckboxgroup(c("White", 
                  "American Indian and Alaska Native", 
                  "Asian", 
                  "Black or African American", 
                  "Native Hawaiian and Other Pacific Islander"),
                container=f)
bg <- ggroup(cont=g); addSpring(bg)
b <- gbutton("Go", container=bg)
enabled(b) <- FALSE
##
addHandlerChanged(cb, handler=function(h,...) {
  enabled(b) <- length(svalue(h$obj)) > 0
})
##
visible(win) <- TRUE


## @knitr remove, eval=FALSE
## if(gconfirm(c("Remove x", "this can't be undone")))
##   rm("x")


## @knitr about
about <- "
A simple GUI to simplify the loading and unloading of packages.
This GUI uses `gcheckboxgroup`, with its `use.table` argument, to
present the user with familiar checkboxes to indicate selection.
Some indexing jujitsu is needed to pull out which value is checked to
trigger the event.
"


## @knitr installed
installed <- installed.packages() ## matrix
installed_packages <- installed[, "Package"]


## @knitr 
package_status <- function() {
  ## Return if package is loaded
  installed_packages %in% loadedNamespaces()
}


## @knitr 
w <- gwindow("package manager", visible=FALSE)
g <- gvbox(cont=w)
g$set_borderwidth(10L)


## @knitr 
a <- package_status()
tbl <- gcheckboxgroup(installed_packages, checked=package_status(),
                      use.table=TRUE,
                      expand=TRUE, container=g)


## @knitr 
bg <- ggroup(cont=g)
addSpring(bg)
gbutton("About", container=bg, handler=function(...) {
  w1 <- gwindow("About", parent=w, visible=FALSE)
  g <- gvbox(container=w1); g$set_borderwidth(10)
  glabel(about, container=g, expand=TRUE)
  gseparator(container=g)
  bg <- ggroup(cont=g)
  addSpring(bg)
  gbutton("dismiss", cont=bg, handler=function(h,...) {
    dispose(w1)
  })
  visible(w1) <- TRUE
})


## @knitr 
visible(w) <- TRUE


## @knitr 
update_tbl <- function(...) {
  blockHandlers(tbl)
  on.exit(unblockHandlers(tbl))
  
  svalue(tbl, index=TRUE) <- package_status()
}


## @knitr 
addHandlerChanged(tbl, handler=function(h, ...) {
  ind <- svalue(h$obj, index=TRUE)
  old_ind <- which(package_status())

  if(length(x <- setdiff(old_ind, ind))) {
    message("detach ", installed_packages[x])
    pkg <- sprintf("package:%s", installed_packages[x])
    detach(pkg, unload=TRUE, character.only=TRUE)
  } else if (length(x <- setdiff(ind, old_ind))) {
    require(installed_packages[x], character.only=TRUE)
  }
  update_tbl()
})


## @knitr 
about <- "GUI to upgrade installed packages"


## @knitr setCRAN



## @knitr 
repos <- getOption("repos")
repos["CRAN"] <- "http://streaming.stat.iastate.edu/CRAN/"
options(repos = repos)
#
pkg <- old.packages()[,c("Package", "Installed", "ReposVer")]


## @knitr 
w <- gwindow("Upgrade installed packages", visible=FALSE)
g <- gvbox(container=w)
g$set_borderwidth(10)


## @knitr 

fg <- ggroup(container=g)
glabel("Filter by:", container=fg)
fltr <- gedit("", initial.msg="Filter by regexp", 
              container=fg)
tbl <- gtable(pkg, chosen.col=1, multiple=TRUE, 
              container=g, expand=TRUE)


## @knitr 
bg <- ggroup(container=g); addSpring(bg)
gbutton("About", container=bg, handler=function(h,...) {
  w1 <- gwindow("About", parent=w, visible=FALSE)
  g <- gvbox(container=w1); g$set_borderwidth(10)
  glabel(about, container=g, expand=TRUE)
  bg <- ggroup(container=g); addSpring(bg)
  gbutton("dismiss", container=bg, 
          handler=function(h,...) dispose(w1))
  visible(w1) <- TRUE
})


## @knitr 
update_btn <- gbutton("Update selected", container=bg, 
                      handler=function(h,...) {
  pkgs <- svalue(tbl)
  if(length(pkgs) == 0) return()
  
  sapply(pkgs, install.packages)
  ## update pkg, update table, clear filter
  tbl[] <- pkg <<-  
    old.packages()[,c("Package", "Installed", "ReposVer")]
  svalue(fltr) <- ""
})
enabled(update_btn) <- FALSE
#
visible(w) <- TRUE


## @knitr 
addHandlerKeystroke(fltr, handler=function(h,...) {
  regexp <- svalue(h$obj)
  if(nchar(regexp) > 0 && regexp != "") {
    ind <- grepl(regexp, pkg[, 'Package'])
    visible(tbl) <- ind
  } else {
    visible(tbl) <- rep(TRUE, nrow(pkg))
  }
})


## @knitr 
addHandlerSelectionChanged(tbl, handler=function(h,...) {
  enabled(update_btn) <- length(svalue(h$obj))
})


