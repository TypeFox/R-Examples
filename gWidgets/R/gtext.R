##' @include guiComponents.R

##' multi-line text edit class
setClass("gText",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' Multiline text edit constructor
##'
##' @exports
gtext <- function(
                  text = NULL, width = NULL, height = 300, font.attr = NULL,
                  wrap = TRUE,
                  handler = NULL, action = NULL, container = NULL,      ... ,
                  toolkit=guiToolkit()){
  widget <- .gtext(toolkit,
                   text=text, width=width, height=height, font.attr=font.attr, wrap=wrap,
                   handler=handler, action=action, container=container ,...
                   )
  obj <- new( 'gText',widget=widget,toolkit=toolkit) 
  return(obj)
}
 

##' generic for toolkit dispatch
##' @alias gtext
setGeneric( '.gtext' ,
           function(toolkit,
                    text = NULL, width = NULL, height = 300, font.attr = NULL,
                    wrap = TRUE,
                    handler = NULL, action = NULL, container = NULL,... )
           standardGeneric( '.gtext' ))


################### methods ###############################

################### insert ###############################
##' generic to insert text into gtext widget
setGeneric("insert",function(obj,value,
                             where = c("end","beginning","at.cursor"),
                             font.attr=NULL,
                             do.newline=TRUE,
                             ...) standardGeneric("insert"))

##' insert text method
setMethod("insert",signature(obj="gText"),
          function(obj, value, where = c("end","beginning","at.cursor"), font.attr = NULL,
                   do.newline = TRUE, ...) {
            toolkit = obj@toolkit
            where = match.arg(where)
            .insert(obj@widget, toolkit, value, where, font.attr,do.newline,...)
          })

##' dispatch with toolkit
##' @alias insert
setGeneric(".insert",function(obj, toolkit,value, where=c("end","beginning","at.cursor"),
                              font.attr=NULL, do.newline=TRUE,...)
           standardGeneric(".insert"))


