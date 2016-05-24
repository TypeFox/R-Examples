##' A class to record the toolkit a gui object uses
setClass("guiWidgetsToolkit",
         representation(toolkit="character"),
         prototype(toolkit="")
         )

##' toolkit class for RGtk2
setClass("guiWidgetsToolkitRGtk2",
         contains="guiWidgetsToolkit",
         prototype=prototype(new("guiWidgetsToolkit"))
         )

##' toolkit class for rJava
setClass("guiWidgetsToolkitrJava",
         contains="guiWidgetsToolkit",
         prototype=prototype(new("guiWidgetsToolkit"))
         )

##' toolkit class for SJava
setClass("guiWidgetsToolkitSJava",
         contains="guiWidgetsToolkit",
         prototype=prototype(new("guiWidgetsToolkit"))
         )
##' toolkit class for  tcltk
setClass("guiWidgetsToolkittcltk",
         contains="guiWidgetsToolkit",
         prototype=prototype(new("guiWidgetsToolkit"))
         )

##' toolkit class for  RwxWidgets
setClass("guiWidgetsToolkitRwxWidgets",
         contains="guiWidgetsToolkit",
         prototype=prototype(new("guiWidgetsToolkit"))
         )


##' toolkit class for qtbase
setClass("guiWidgetsToolkitQt",
         contains="guiWidgetsToolkit",
         prototype=prototype(new("guiWidgetsToolkit"))
         )


##################################################


##' all packages that are registered with gWidgets. Used if guiToolkit not specified
registered_packages <- c("gWidgetsRGtk2", "gWidgetstcltk", "gWidgetsQt", "gWidgetsrJava",
                         "gWidgetsRwxWidgets")

##' set or get the current toolkit for gWidgets
guiToolkit <- function(name=NULL) {
  ## plan, if name is NULL, and options("guiToolkit") NULL then we popup a menu
  ## with choices coming from all installed packages named gWidgetsXXXX
  ## when a name is selected, we require the package gWidgets+name

  if(is.null(name)) {
    ## try to get from inheritance, then get from option

    x = try(get("toolkit", inherits=TRUE), silent=TRUE)
    if(!inherits(x,"try-error")) {
      ## check that toolkit is of guiWidgets type
      x = try("x@toolkit", silent=TRUE)
      if(!inherits(x,"try-error"))
        name = x
      else
        name = getOption("guiToolkit")
    } else {
      name = getOption("guiToolkit")
    }
  }
  if(!is.null(name) && is.na(name)) return(NULL)          # use NA to override choice
  ## no if it is null, we have to find the possible choices
  if(is.null(name)) {

    f <- function(x) !inherits(try(find.package(x), silent=TRUE), "try-error")
    choices <- Filter(f, registered_packages)
    
    if(interactive()) {
      if(length(choices) == 0) {
        warning("No toolkits installed")
        return(NULL)
      } else if(length(choices) == 1) {
        theChoice = choices
      } else {
        theChoice = menu(choices, title="Select a GUI toolkit")
        if(theChoice == 0) {
          warning("No toolkit loaded")
          return(NULL)
        } else {
          theChoice = choices[theChoice]
        }
      }
      ## go with theChoice
      name = gsub("^gWidgets","",theChoice)
      options("guiToolkit"=name)
    } else {
      ## not interactive 
      return(NULL)
    }
  }

  ## require the package
  require(paste("gWidgets",name,sep=""), character.only=TRUE)
  ## we return an instance of the toolkit class
  obj = new(paste("guiWidgetsToolkit",name,sep=""), toolkit = name)
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
