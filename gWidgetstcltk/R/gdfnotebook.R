setClass("gDfNotebooktcltk",
         representation = representation(
           gnotebook="guiWidget"
           ),
         contains="gNotebooktcltk"
         )

setMethod(".gdfnotebook",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   items = NULL,
                   container = NULL,
                   ... # passed to Group, gnotebook = nb,
                                        # notebook = nb$notebook)
    ) {

            force(toolkit)

            return(glabel("gdfnotebook not available", container=container)@widget)
  
          })

##################################################
##
## gWidgetMethods (inherits others from gnotebook
## object is name of R object *or* file

## REWRITE me to dispatch on value. This first part is ugly and broken
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gDfNotebooktcltk"),
          function(obj, toolkit, value, ...) {
          })
          

