##' @include guiComponentWithItems.R

##' Radio button class
setClass("gRadio",
         contains="guiComponentWithItems",
         prototype=prototype(new("guiComponentWithItems"))
         )
##' Constructor for radio button widget
gradio = function(items,selected=1, horizontal=FALSE, handler=NULL,
  action=NULL, container=NULL, ...,
  toolkit=guiToolkit()) {

  ## check input
  if(length(x <- unique(items) ) != length(items))
    gwCat("Using unique items for selection values")
  
  radio = .gradio(toolkit, x, selected, horizontal, handler, action, container,...)
  obj = new("gRadio",widget=radio, toolkit=toolkit)
  return(obj)
}

##' method for dispatch into toolkit
setGeneric(".gradio",function(toolkit,
                              items, selected=1, horizontal=FALSE, handler=NULL, action=NULL,
                              container=NULL,
                              ...) standardGeneric(".gradio"))


##' replace selections
##'
##' Ensure values are unique
setReplaceMethod("[",signature(x="gRadio"),
          function(x, i, j,...,value) {
            ## check input
            if(length(value) != length(value <- unique(value)))
              gwCat("Using unique values for selection values")

            callNextMethod(x, i, j, ..., value=value)
          })
