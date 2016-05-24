setClass("gRadiotcltk",
         representation = representation("gComponentR5tcltk",
           coercewith="NULLorFunction"),
         contains="gComponentR5tcltk",
         prototype=prototype(new("gComponentR5tcltk"))
         )

## constructor
setMethod(".gradio",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   items, selected=1, horizontal=FALSE,
                   handler=NULL, action=NULL,
                   container=NULL,       
                   ...
                   ) {
            force(toolkit)

           

            if(is(container,"logical") && container)
              container = gwindow()
            if(!is(container,"guiWidget")) {
              warning(gettext("Container is not correct. No NULL containers possible\n" ))
              return()
            }

            tt = getWidget(container)

            ## use coerce with
            theArgs = list(...)
            if(!is.null(theArgs$coerce.with)) {
              coerce.with = theArgs$coerce.with
            } else {
              if(is.numeric(items))
                coerce.with = as.numeric
              else if(is.logical(items))
                coerce.with = as.logical
              else
                coerce.with = as.character
            }
            if(is.character(coerce.with))
              coerce.with = get(coerce.with)

            items <- as.character(items)


            rb_widget <- getRefClass("RadioButton")$new(parent=tt, items=items, horizontal=horizontal)
            
            obj = new("gRadiotcltk",block=rb_widget$get_widget(), widget=rb_widget$get_widget(),
              R5widget=rb_widget,
              toolkit=toolkit, ID=getNewID(), e = new.env(),
              coercewith = coerce.with)


            svalue(obj, index=TRUE) <- selected
            
            ## add to container
            add(container,  obj,...)
  
            ## add handler
            if(!is.null(handler))
              addhandlerchanged(obj, handler, action)

            
            invisible(obj)
          })

## methods
## setMethod(".svalue",
##           signature(toolkit="guiWidgetsToolkittcltk",obj="gRadiotcltk"),
##           function(obj, toolkit, index=NULL, drop=NULL, ...) {

##             rb_widget <- obj@R5widget
##             index <- getWithDefault(index, FALSE)
##             if(index) {
##               return(rb_widget$get_index())
##             } else {
##               val <- rb_widget$get_value()
##               if(!is.null(obj@coercewith))
##                 return(obj@coercewith(val))
##               else
##                 return(val)
##             }
          
              
##           })

## svalue<-
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gRadiotcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {

                   if(is.data.frame(value))
                     value <- value[,1, drop=TRUE]
                   
                   rb_widget <- obj@R5widget                   
                   index <- getWithDefault(index, FALSE)
                   if(index) {
                     rb_widget$set_index(value)
                   } else {
                     rb_widget$set_value(value)
                   }
                   return(obj)
                
                   
                   return(obj)
                 })


## setMethod(".leftBracket",
##           signature(toolkit="guiWidgetsToolkittcltk",x="gRadiotcltk"),
##           function(x, toolkit, i, j, ..., drop=TRUE) {
##             ## return(items)
##             rb_widget <- x@R5widget
##             items <- rb_widget$get_items()
##             if(missing(i))
##               items
##             else
##               items[i]
            
         
##           })
            
## setMethod("[",
##           signature(x="gRadiotcltk"),
##           function(x, i, j, ..., drop=TRUE) {
##             .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
##           })


## ## This sets the labels for the buttons
## ## add in markup here.
## setReplaceMethod(".leftBracket",
##           signature(toolkit="guiWidgetsToolkittcltk",x="gRadiotcltk"),
##           function(x, toolkit, i, j, ..., value) {

##             rb_widget <- x@R5widget
            
##             if(!missing(i)) {
##               items <- rb_widget$get_items()
##               items[i] <- value
##               value <- items
##             }
##             rb_widget$set_items(value)
##             return(x)
            
##           })

## setReplaceMethod("[",
##                  signature(x="gRadiotcltk"),
##                  function(x, i, j,..., value) {
##                    .leftBracket(x, x@toolkit, i, j, ...) <- value
##                    return(x)
##                  })

## setMethod(".length",
##           signature(toolkit="guiWidgetsToolkittcltk",x="gRadiotcltk"),
##           function(x,toolkit) {
##             rb_widget <- x@R5widget
##             rb_widget$no_items()
##             ##length(tag(x,"items"))
##           })

## ## inherited enabled isn't workgin                
## setReplaceMethod(".enabled",
##                  signature(toolkit="guiWidgetsToolkittcltk",obj="gRadiotcltk"),
##                  function(obj, toolkit, ..., value) {
##                    rb_widget <- obj@R5widget
##                    rb_widget$set_enabled(value)
##                    return(obj)
                  
##                  })


##################################################
## handlers


##' only one handler per widget
##'
##' This could be changed, but only if asked ...
##' This does not get called by svalue -- it should?
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gRadiotcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            rb_widget <- obj@R5widget            
            user.data=list(obj=obj, handler=handler, action=action)
##            id <- rb_widget$add_handler("<ButtonRelease-1>",
            id <- rb_widget$add_handler("command",
                                        handler=function(user.data) {
                                          h <- user.data[c("obj", "action")]
                                          user.data$handler(h)
                                  },
                                        user.data=user.data)
            invisible(id)
          })       
            

            

## click and changed the same
setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gRadiotcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addhandlerchanged(obj,toolkit,handler,action,...)
          })


