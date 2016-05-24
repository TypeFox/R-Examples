## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

### DATA

#    quickTool(playState,
#              label = "Data...",
#              icon = "gtk-cdrom",
#              tooltip = "View and edit attached data objects",
#              f = data_handler,
#              show = length(playState$env) > 0)

data_handler <- function(widget, playState)
{
    lstObjects <- function(envir= .GlobalEnv) {
        objlist <- ls(envir=envir)
        objclass <- sapply(objlist, function(objName) {
            obj <- get(objName, envir=envir)
            class(obj)[1]
        })
        data.frame(Name = I(objlist), Class = I(objclass))
    }
    browseEnv2 = function(envir = .GlobalEnv) {
        listOfObjects <- lstObjects(envir=envir)
        gtable(listOfObjects, container = gwindow("Data objects (double-click to edit)"),
               handler = function(h,...) {
                   myName <- svalue(h$obj)
                   oldObj <- get(myName, envir=envir)
                   newObj <- edit(oldObj)
                   if (identical(newObj, oldObj)) return()
                   assign(myName, newObj, envir=envir)
                   gmessage(paste("Edited object", myName,
                                  "-- you might want to reload the plot."))
               })
    }
    browseEnv2(playState$env)
}
