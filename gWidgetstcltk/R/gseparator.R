##################################################
## add a separator to a container. Needs the container

##  inspired by
##  http://search.cpan.org/src/WGDAVIS/Tk-Separator-0.50/Separator.pm

setClass("gSeparatortcltk",
         contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )

## should this return object?
setMethod(".gseparator",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   horizontal = TRUE, container = NULL, ...) {

            force(toolkit)

            ## if null, we return a stub. Useful for gmenu, gtoolbar
            if(is.null(container)) {
              gp <- ttkframe(.TkRoot)   # empty
              obj <- new("gSeparatortcltk", block=gp, widget=gp,
                toolkit=toolkit, ID=getNewID(), e = new.env())
              return(obj)
            }

            
            if(is(container,"logical") && container)
              container = gwindow()
            if(!is(container,"guiWidget")) {
              warning("Container is not correct. No NULL containers possible\n" )
              return()
            }

            
            theArgs = list(...)
            if(!is.null(theArgs$col))
              col = theArgs$col
            else
              col = "black"

            tt <- getWidget(container)
##            gp <- ttkframe(tt)

            if(horizontal)
              orient <- "horizontal"
            else
              orient <- "vertical"
            sep <- ttkseparator(tt, orient=orient)

            ## if(horizontal)
            ##   tkpack(sep)#, expand=TRUE, fill="x")
            ## else
            ##   tkpack(sep)#, expand=TRUE, fill="y")
            
            obj = new("gSeparatortcltk", block=sep, widget=sep,
              toolkit=toolkit, ID=getNewID(), e = new.env())


            ## add gp to container. Fixe expand argument to be TRUE
#            theArgs$expand = TRUE
            theArgs$obj <- container
            theArgs$value <- obj
            do.call("add", theArgs)
#            add(container, obj, ...)


            invisible(obj)
            
          })

## no size method
setReplaceMethod(".size", 
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gSeparatortcltk"),
                 function(obj, toolkit, ..., value) {
                   gwCat("No size<- method for separators")
                   return(obj)
                 })


## setMethod(".add",
##           signature(toolkit="guiWidgetsToolkittcltk", obj="gLayouttcltk",
##                     value="gSeparatortcltk"),
##           function(obj, toolkit, value, ...) {
##           })

.isgSeparator <- function(obj) {
  (is(obj,"guiComponent") && is(obj@widget,"gSeparatortcltk") ) ||
    is(obj,"gSeparatortcltk")
}



