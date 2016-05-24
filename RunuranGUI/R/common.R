#############################################################################
##
##  Auxiliary functions
##
#############################################################################


## --------------------------------------------------------------------------
## Message

RunuranGUI.TESTING <- FALSE

mygmessage <- function(msg, title, icon) {

  ## For running tests non-interactively we need non-modal dialog widgets.
  ## Thus we use the global switch 'RunuranGUI.TESTING'.

  if (exists("RunuranGUI.TESTING") && isTRUE(RunuranGUI.TESTING))
    ## non-modal box for testing
    galert(msg, title=title, delay=1)
  else
    ## modal dialog box for installed package
    gmessage(msg, title=title, icon=icon)

}


## --------------------------------------------------------------------------
## Error message

error.message <- function(msg, title="UNU.RAN - Error") {
  ## galert(message,title=title,delay=100)
  mygmessage(msg, title=title, icon="error")
}


## --------------------------------------------------------------------------
## Internal error

internal.error <- function(msg="unknown") {
  text <- paste("Internal error!\n\n", msg, "\n\nPlease report.",sep="")
  mygmessage(text, title="Internal error",icon="error")
  stop (text)
}


## --------------------------------------------------------------------------
## Toggle widget:  added <--> deleted

toggle.group <- function(obj,group,add) {
  if ( tag(group,"added") != add ) {
    if (isTRUE(add)) {
      add(obj,group)
      tag(group,"added") <- TRUE
    }
    else {
      delete(obj,group)
      ## svalue(group,index=TRUE) <- 1
      tag(group,"added") <- FALSE
    }
  }
}


## --------------------------------------------------------------------------
## List of arguments

function.args <- function(func) {

  ## parse arguments of function
  fargs <- formals(func)

  ## store in list
  args <- list()
  
  for (p in names(fargs)) {
    args[[p]] <-
      paste( deparse(fargs[[p]]), sep="", collaps="")
  }

  return (args)
}

## --------------------------------------------------------------------------
## Create help page

show.help.pages <- function(pages=list()) {

  ## check arguments
  if (length(pages)<=0)
    return (NULL)

  ## create a new window of given size
  win <- gwindow("Runuran - Help")
  size(win) <- c(600, 400)
  group <- ggroup(horizontal=FALSE, spacing=10, container=win)

  ## create notebook
  nb <- gnotebook(container=group, expand=TRUE, closebuttons=FALSE)

  ## add page for R code
  for (p in pages) {
    ## get help text
    help.txt <- get.help.text(p)
    ## create page
    if (! is.null(help.txt)) 
      pcode <- gtext(text=help.txt, label=p, container=nb)
  }
  svalue(nb) <- 1
  
  ## add 'close' button
  closeGroup = ggroup(container=group)
  addSpring(closeGroup)
  gbutton("close", handler=function(h,...){dispose(win)},
          container=closeGroup)
}

## Obtain help text in relevant packages.
get.help.text <- function(topic) {
  ## search for help message in relevant packages
  pkgname <- NULL
  for (pkg in c("Runuran", "RunuranGUI")) {
    out <- do.call("help", list(topic=topic, package=pkg))
    if (length(out) != 0) {
      pkgname <- pkg
      break
    }
  }
  if (is.null(pkgname)) return(NULL)

  ## read help text 
  help.txt <- ""  ## keep R CMD check happy
  help.con <- textConnection("help.txt", "w", local = TRUE)
  tools::Rd2txt(utils:::.getHelpFile(out), out=help.con, package=pkgname,
                width=80L)
  close(help.con)

  ## convert output
  help.txt <- gsub("_\b","", paste(help.txt, collapse="\n"))
  
  return (help.txt)
}


## --------------------------------------------------------------------------
## Handler for help page

show.help <- function(h,...) {
  show.help.pages(h$action)
}
       

## --------------------------------------------------------------------------
## Wraps given string with "function(x) { ... }";
## to be used as argument for 'coerce.with' in gedit().

wrap.function.body.NA <- function(x) {
  ## returns NA in case of an empty field
  if (gsub("\\s+", "", x, perl=TRUE) == "")
    NA_character_
  else
    paste("function(x){",x,"}", sep="")
}

wrap.function.body.NULL <- function(x) {
  ## returns NULL in case of an empty field
  if (gsub("\\s+", "", x, perl=TRUE) == "")
    NULL
  else
    paste("function(x){",x,"}", sep="")
}


## --------------------------------------------------------------------------
## Create check box with edit field which is disabled when the checkbox
## is not checked

gcheckedit <- function(label, checked=TRUE, text, ttip=NULL, coerce.with=NULL, width=25, container) {
  ## label     ... label for checkbox
  ## checked   ... whether checkbox is checked by default
  ## text      ... initial text in edit widget
  ## ttip      ... tooltip for edit widget
  ## coerce.with . see gedit()
  ## width     ... width of gedit field
  ## container ... container to attach widget to

  ## edit widget 
  edt <- gedit(text=text, coerce.with=coerce.with, width=width, container=container)
  tag(edt,"label") <- label
  enabled(edt) <- tag(edt,"enabled") <-
    if(isTRUE(checked)) TRUE else FALSE
  if(! is.null(ttip)) tooltip(edt) <- ttip

  ## check box
  cbx <- gcheckbox(label, checked=checked, container=container,
                   handler=function(h,...) {
                     enabled(edt) <- tag(edt,"enabled") <-
                       if(isTRUE(svalue(h$obj))) TRUE else FALSE
                   })

  return (list(cbx=cbx,edt=edt))
}

## --------------------------------------------------------------------------
