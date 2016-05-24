##' A class to record the toolkit a gui object uses
##'
##' @exportClass guiWidgetsToolkit
##' @rdname S4-classes
##' @name guiWidgetsToolkit-class
setClass("guiWidgetsToolkit",
         representation(toolkit="character"),
         prototype(toolkit="")
         )


## Subclasses are done in the toolkit packages. For example, this is done in gWidgets2RGtk2:
## setClass("guiWidgetsToolkitRGtk2",
##          contains="guiWidgetsToolkit")

##################################################


##' set or get the current toolkit for gWidgets
##'
##' @param name name of toolkit (e.g. "tcltk", "RGtk2", "Qt" (not
##' qtbase)). If NULL, then we search for it in a) an iherited toolkit
##' object b) the "guiToolkit" option (which can be set via
##' \code{options("guiToolkit"="RGtk2")}, say. If that fails, we
##' prompt for a selection for any installed toolkit.  In the typical
##' usage, this all happens in the background, except perhaps once.
##'
##' In design this is to allow different toolkits to be used with
##' different GUIs, but due to differences in event loops, this often
##' leads to lockups, so is not recommended.
##' @return an instance of guiWidgetsToolkit sub class.
##' @export
guiToolkit <- function(name=NULL) {
  ## plan, if name is NULL, and options("guiToolkit") NULL then we popup a menu
  ## with choices coming from all installed packages named gWidgetsXXXX
  ## when a name is selected, we require the package gWidgets+name

  if(missing(name) || is.null(name)) {
    ## try to get from inheritance, then get from option

    x = try(get("toolkit", inherits=TRUE), silent=TRUE)
    if(is(x, "guiWidgetsToolkit")) {
      name <- x
    } else {
      name <- getOption("guiToolkit")
    }
  }

  if(!is.null(name) && is.na(name)) {
    message("Choice overridden")
    return(NULL)          # use NA to override choice
  }

  
  ## no if it is null, we have to find the possible choices
  if(is.null(name)) {

    ## A list of possible packages
    poss_packages <- c("gWidgets2RGtk2",
                       "gWidgets2tcltk",
                       "gWidgets2Qt",
                       "gWidgets2rJava",
                       "gWidgets2wxWidgets",
                       "gWidgets2WWW")

    
    f <- function(x) !inherits(try(find.package(x), silent=TRUE), "try-error")
    choices <- Filter(f, poss_packages)

    
    if(interactive()) {
      if(length(choices) == 0) {
        message("No toolkit packages are installed.")
        f <- system.file("install/installing_toolkits.txt", package="gWidgets2")
        cat(paste(f, "\n"))
        return(NULL)
      } else if(length(choices) == 1) {
        theChoice <- choices
      } else {
        theChoice <- menu(choices, title="Select a GUI toolkit")
        if(theChoice == 0) {
          warning("No toolkit loaded")
          return(NULL)
        } else {
          theChoice <- choices[theChoice]
        }
      }
      ## go with theChoice
      name <- gsub("^gWidgets2","",theChoice)
      options("guiToolkit"=name)

    } else {
      ## not interactive 
      return(NULL)
    }
  }

  ## override
  if(name == "qtbase")
    name <- "Qt"

  if(grepl("^2", name)) {
    message("Why is there a 2 in name")
    name <- gsub("^2", "", name)
  }

  
  ## require the package
  require(sprintf("gWidgets2%s", name, sep=""), character.only=TRUE)


  ## check for headless Gtk
  if (name == "RGtk2") {
      gtk_initialized <- eval(parse(text=sprintf("RGtk2:::%s", ".gtkInitCheck()")))
      if (!gtk_initialized)
          stop("Can't load RGtk2")
  }

  
  ## we return an instance of the toolkit class
  obj <- new(sprintf("guiWidgetsToolkit%s", name), toolkit = name)
  return(obj)
}

##' Which toolkit are we using?
##'
##' @return string of toolkit (RGtk2, tcltk, Qt, ???)
##' @export
gtoolkit <- function() {
  obj <- guiToolkit()
  obj@toolkit
}
