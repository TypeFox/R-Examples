## build widget based on gcheckbox
setClass("gCheckboxgrouptcltk",
         representation = representation("gComponentR5tcltk",
           coercewith="NULLorFunction"),
         contains="gComponentR5tcltk",
         prototype=prototype(new("gComponentR5tcltk"))
         )

setMethod(".gcheckboxgroup",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   items, checked = FALSE,
                   horizontal=FALSE, use.table=FALSE,
                   handler = NULL, action = NULL, container = NULL, ...) {

            force(toolkit)
            
            if(missing(items) || length(items) == 0)
              stop("Need items to be a vector of items")

            if(is(container,"logical") && container)
              container = gwindow()
            if(!is(container,"guiWidget")) {
              warning(gettext("Container is not correct. No NULL containers possible\n" ))
              return()
            }
            

            if(use.table) {
              obj <- .gcheckboxgrouptable(toolkit, items=items, checked=checked,
                                   handler=handler, action=action, container=container, ...)
              return(obj)
            }


            theArgs = list(...)
            if(!is.null(theArgs$coerce.with)) {
              coerce.with = theArgs$coerce.with
            } else {
              if(is.numeric(items))
                coerce.with = as.numeric
              else
                coerce.with = as.character
            }
            if(is.character(coerce.with))
              coerce.with = get(coerce.with)


            tt = getWidget(container)

            cbg_widget <- getRefClass("CheckButtonGroup")$new(parent=tt, items=items,
                                                              selected=checked, horizontal=horizontal)

            obj <- new("gCheckboxgrouptcltk", block=cbg_widget$get_widget(), widget=cbg_widget$get_widget(),
                       R5widget=cbg_widget,
                       toolkit=toolkit, coercewith = coerce.with, e = new.env())

            svalue(obj) <- checked
            
            
            
            ## add to container
            add(container,  obj,...)
  
            ## add handler
            if(!is.null(handler))
              tag(obj, "handler.id") <- addhandlerchanged(obj, handler, action)

            invisible(obj)
          })


### methods
## setMethod(".svalue",
##           signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxgrouptcltk"),
##           function(obj, toolkit, index=NULL, drop=NULL, ...) {
            
##             cbg_widget <- obj@R5widget
##             index <- getWithDefault(index, FALSE)
##             if(index) {
##               return(cbg_widget$get_index())
##             } else {
##               val <- cbg_widget$get_value()
##               if(!is.null(obj@coercewith))
##                 return(obj@coercewith(val))
##               else
##                 return(val)
##             }
          
            

##           })

## toggles state to be T or F
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxgrouptcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {

                   cbg_widget <- obj@R5widget
                   index <- getWithDefault(index, FALSE)

                   if(index) {
                     cbg_widget$set_index(value)
                   } else if(is.logical(value)) {
                     n <- length(obj)
                     value <- rep(value, length.out=n)
                     cbg_widget$set_index(which(value))
                   } else {
                     cbg_widget$set_value(value)
                   }
                   return(obj)
                 })

## [ and [<- refer to the names -- not the TF values

## setMethod("[",
##           signature(x="gCheckboxgrouptcltk"),
##           function(x, i, j, ..., drop=TRUE) {
##             .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
##           })
## setMethod(".leftBracket",
##           signature(toolkit="guiWidgetsToolkittcltk",x="gCheckboxgrouptcltk"),
##           function(x, toolkit, i, j, ..., drop=TRUE) {
##             cbg_widget <- x@R5widget
##             items <- cbg_widget$get_items()
##             if(missing(i))
##               items
##             else
##               items[i]
        
##           })

## assigns names
## setReplaceMethod("[",
##                  signature(x="gCheckboxgrouptcltk"),
##                  function(x, i, j,..., value) {
##                    .leftBracket(x, x@toolkit, i, j, ...) <- value
##                    return(x)
##                  })

## setReplaceMethod(".leftBracket",
##           signature(toolkit="guiWidgetsToolkittcltk",x="gCheckboxgrouptcltk"),
##           function(x, toolkit, i, j, ..., value) {

##             cbg_widget <- x@R5widget
##             if(!missing(i)) {
##               items <- cbg_widget$get_items()
##               items[i] <- value
##               value <- items
##             }
##             cbg_widget$set_items(value)

##              return(x)
##           })


## setMethod(".length",
##           signature(toolkit="guiWidgetsToolkittcltk",x="gCheckboxgrouptcltk"),
##           function(x,toolkit) {

##             cbg_widget <- x@R5widget
##             cbg_widget$no_items()
## #            length(tag(x,"items"))
##           })


## ## inherited enabled isn't workgin                
## setReplaceMethod(".enabled",
##                  signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxgrouptcltk"),
##                  function(obj, toolkit, ..., value) {

##                    cbg_widget <- obj@R5widget                   
##                    cbg_widget$set_enabled(value)
##                    return(obj)
                  
##                  })


## This handler code is common to gradio and gcheckboxgroup. Should abstract out into a superclass.
## IF we do that, we should also use CheckButton bit
setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxgrouptcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {

            cbg_widget <- obj@R5widget            
            user.data=list(obj=obj, handler=handler, action=action)
##            id <- cbg_widget$add_handler("<ButtonRelease-1>",
            id <- cbg_widget$add_handler("command",
                                        handler=function(user.data) {
                                          h <- user.data[c("obj", "action")]
                                          user.data$handler(h)
                                  },
                                        user.data=user.data)
            invisible(id)
            
          })
## clicked is changed
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxgrouptcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addhandlerclicked(obj, toolkit, handler, action, ...)
          })

