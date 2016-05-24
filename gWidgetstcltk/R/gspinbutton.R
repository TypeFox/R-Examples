## Could make spinbutton slider, subclass as methods are identical
## setClass("gSpinbuttontcltk",
##          contains="gComponenttcltk",
##          prototype=prototype(new("gComponenttcltk"))
##          )

setClass("gSpinbuttontcltk",
         representation = representation("gComponentR5tcltk"),
         contains="gComponentR5tcltk",
         prototype=prototype(new("gComponentR5tcltk"))
         )

setMethod(".gspinbutton",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   from=0,to=10,by=1,value=from,digits=0,
                   handler=NULL, action=NULL,
                   container=NULL, ...) {

            force(toolkit)

            
            if(is(container,"logical") && container)
              container = gwindow()
            if(!is(container,"guiWidget")) {
              warning("Container is not correct. No NULL containers possible\n" )
              return()
            }

            tt = getWidget(container)
            sp_widget <-  getRefClass("SpinButton")$new(parent=tt, from=from, to=to, by=by)

            obj = new("gSpinbuttontcltk",block=sp_widget$get_widget(), widget=sp_widget$get_widget(),
              R5widget=sp_widget,
              toolkit=toolkit, ID=getNewID(), e = new.env())

            ## add to container
            add(container,  obj,...)
            
            ## add handler
            if(!is.null(handler))
              addhandlerchanged(obj, handler, action)

            invisible(obj)
            
            ## ## no spinbutton in the tcltk
            ## vals =  as.character(seq(from,to,by=by))
            

            ## tt = getWidget(container)
            ## gp = ttkframe(tt)
            
            ## sb = tkwidget(gp, "spinbox", from=from, to=to, increment=by)
            ## tcl(sb,"set",value)
            ## tkpack(sb, expand=TRUE, fill="both")
            
            ## obj = new("gSpinbuttontcltk",block=gp, widget=sb,
            ##   toolkit=toolkit, ID=getNewID(), e = new.env())
            
            ## add(container, obj,...)
            
            ## if (!is.null(handler))  {
            ##   id = addhandlerchanged(obj, handler, action)
            ## }
            
            ## invisible(obj)
          })

## ### methods
## setMethod(".svalue",
##           signature(toolkit="guiWidgetsToolkittcltk",obj="gSpinbuttontcltk"),
##           function(obj, toolkit, index=NULL, drop=NULL, ...) {
##             sp_widget <- obj@R5widget
##             sp_widget$get_value()

##             ## sb = getWidget(obj)
##             ## val = as.numeric(tcl(sb,"get"))
##             ## return(val)
##           })

## setReplaceMethod(".svalue",
##                  signature(toolkit="guiWidgetsToolkittcltk",obj="gSpinbuttontcltk"),
##                  function(obj, toolkit, index=NULL, ..., value) {
##                    sp_widget <- obj@R5widget

##                    sp_widget$set_value(value)
                   
##                    ## sb = getWidget(obj)
##                    ## tcl(sb,"set",value)
##                    return(obj)
##                  })

## enabled -- use tkconfigure, not tcl
setReplaceMethod(".enabled",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gSpinbuttontcltk"),
                 function(obj, toolkit, ..., value) {
                   if(as.logical(value))
#                     tcl(getWidget(obj),"state","!disabled")
                     tkconfigure(getWidget(obj),state="normal")
                   else
#                     tcl(getWidget(obj),"state","disabled")
                     tkconfigure(getWidget(obj),state="disabled")
                   return(obj)
                 })



## ## Method to replace values of spin button
## setReplaceMethod("[",
##                  signature(x="gSpinbuttontcltk"),
##                  function(x, i, j,..., value) {
##                    .leftBracket(x, x@toolkit, i, j, ...) <- value
##                    return(x)
##                  })

## setReplaceMethod(".leftBracket",
##           signature(toolkit="guiWidgetsToolkittcltk",x="gSpinbuttontcltk"),
##           function(x, toolkit, i, j, ..., value) {
##             obj <- x
##             sp_widget <- obj@R5widget
##             sp_widget$set_items(value)


##             ## widget <- getWidget(obj)

##             ## ## check that value is a regular sequence
##             ## if(length(value) <=1) {
##             ##   warning("Can only assign a vector with equal steps, as produced by seq")
##             ##   return(obj)
##             ## }
##             ## if(length(value) > 2 &&
##             ##    !all.equal(diff(diff(value)), rep(0, length(value) - 2))) {
##             ##   warning("Can only assign a vector with equal steps, as produced by seq")
##             ##   return(obj)
##             ## }
##             ## ## get current value, increment
##             ## curValue <- svalue(obj)
##             ## inc <- head(diff(value), n=1)

##             ## tkconfigure(widget, from=min(value), to =max(value), increment=inc)
##             ## tcl(widget, "set", curValue)

##             ## all done
##             return(obj)
##           })



## size has no height
setReplaceMethod(".size", 
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gSpinbuttontcltk"),
                 function(obj, toolkit, ..., value) {
                   width <- ceiling(value[1]/widthOfChar)
                   tkconfigure(getWidget(obj), width=width)
                   return(obj)
                 })

### handlers
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gSpinbuttontcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {

            sp_widget <- obj@R5widget            
            user.data=list(obj=obj, handler=handler, action=action)
            ##            id <- rb_widget$add_handler("<ButtonRelease-1>",
            id <- sp_widget$add_handler("command",
                                        handler=function(user.data) {
                                          h <- user.data[c("obj", "action")]
                                          user.data$handler(h)
                                        },
                                        user.data=user.data)
            invisible(id)



            ##                             #.addhandlerclicked(obj, toolkit, handler, action,...)

            ## changeHandler <- handler

            ## ## need a pause
            ## addhandler(obj,toolkit, signal="<Button-1>",
            ##            action=action, 
            ##            handler = function(h,...) {
            ##                tcl("after",150,function(...) {
            ##                  changeHandler(h,...) ## need to pause
            ##                })
            ##              })

            ## addhandler(obj,toolkit, signal="<Return>",
            ##            action=action, 
            ##            handler = function(h,...) {
            ##              tcl("after",150,function(...) {
            ##                changeHandler(h,...) ## need to pause
            ##              })
            ##            })

          })
