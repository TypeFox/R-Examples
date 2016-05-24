##' @include guiComponentWithItems.R

##' Combobox class
setClass("gCombobox",
         contains="guiComponentWithItems",
         prototype=prototype(new("guiComponentWithItems"))
         )

##' constructor for combobox
##'
##' @export
##' @param items Items to select from. A vector or a data frame. If a
##' data frame, then first column is values. Second is optional, but
##' can specify a stock icon name, third is optional and can be used
##' to specify a tooltip. These may not be supported in all toolkits.
##' @param selected integer. Which item (by index) is selected. Use -1 for no selection
##' @param editable logical. Is user allowed to edit value
##' @param coerce.with A function of function name to be called before
##' selected value is returned by \code{svalue}
##' @param handler Called when combobox value is changed.
##' @param action passed to handler
##' @param container parent container
##' @param ... passed to parent container's \code{add} method
##' @param toolkit toolkit
##' @return Returns an object of class \code{gCombobox} for which the following methods are overriden:
##' \enumerate{
##' \item \code{svalue} Return selected value by name or (if \code{index=TRUE} by index). The latter only if \code{editable=FALSE}.
##' 
##' \item \code{\svalue<-} Set the selected value by value or if \code{index=TRUE} by index.
##'
##' \item \code{[} return items to select from
##'
##' \item \code{[<-} Set items to select from.
##' }
##' @example ~/pmg/r-forge/gwidgets/pkg/gWidgets/inst/tests/ex-gcombobox.R
gcombobox =function(
  items, selected = 1, editable = FALSE, coerce.with=NULL, handler = NULL,      action = NULL, container = NULL, ... ,
  toolkit=guiToolkit()){

  ## adjust items. Pass in data frame
  if(!is.data.frame(items) && !is.matrix(items))
    items <- data.frame(items, stringsAsFactors=FALSE)

  if(!is.data.frame(items))
    items <- as.data.frame(items)

  
  widget =  .gdroplist (toolkit,
    items=items, selected=selected, editable=editable, coerce.with=coerce.with, handler=handler, action=action, container=container, ...
    )
  obj = new( 'gCombobox',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##'
##' @alias gcombobox
setGeneric( '.gdroplist' , function(toolkit,
                                    items, selected = 1, editable = FALSE, coerce.with = NULL, handler = NULL,      action = NULL, container = NULL, ... ) standardGeneric( '.gdroplist' ))

##' Alias for gcombobox. Deprecated
##' 
##' @alias gcombobox
gdroplist <- gcombobox 



##' svalue method for combobox
##'
##' Main property of checkbox group are the states. 
##' @param obj object
##' @param index If \code{TRUE} returns index of selected. Otherwise a character vector.
##' @param drop ignored
##' @return character or index. Some toolkits return \code{-1} when nothing is selected.
##' @exports
setMethod("svalue", signature(obj="gCombobox"),
          function(obj, index=NULL, drop=NULL, ... ) {
            .svalue(obj@widget, obj@toolkit, ...,index=index, drop=drop)            
          })




##' set state for combobox
##'
##' @param obj
##' @param index if \code{TRUE} then an index is expected. If \code{FALSE}, then name should match available unless \code{editable=TRUE} is specified
##' @param ... ignored
##' @param value numeric or character depending on \code{index}
##' @return void
##' @exports
setReplaceMethod("svalue", signature(obj="gCombobox"),
          function(obj, index=NULL, ...,value) {
            .svalue(obj@widget, obj@toolkit, index=index, ...) <- value
            return(obj)
          })

##' replacement method for combobox selection items
##'
##' Ensure that value is a data frame. One can pass a vector or a
##' one-column data frame to inidicate the possible values for
##' selection, a second column is used for an icons (if possible), a
##' third for a tooltip (if possible).
##' get_from_templace for guiWidget
setReplaceMethod("[",signature(x="gCombobox"),
          function(x,i,j,...,value) {
            
            ## adjust value. Pass in data frame
            if(!is.data.frame(value) && !is.matrix(value))
              value <- data.frame(value, stringsAsFactors=FALSE)
            
            if(!is.data.frame(value))
              value <- as.data.frame(value)

            callNextMethod(x, i, j, ..., value=value)
          })

  
