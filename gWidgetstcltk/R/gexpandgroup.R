## expander group, like a group, only expands, contracts if requested
## inherits from ggroup, see ggroup's arguments: horizontal, spacing, container
setClass("gExpandgrouptcltk",
         contains="gGrouptcltk",
         prototype=prototype(new("gGrouptcltk"))
         )


setMethod(".gexpandgroup",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   text="", markup=FALSE,horizontal=TRUE,
                   handler=NULL, action=NULL,
                   container = NULL, ...){

            force(toolkit)


            ## the ... arguments get passed as follows
            ## expand, container, anchor go for outergroup
            ## horizontal, spacing, use.scrollwindow pass to inner group

            oggroup <- function(container,spacing, user.scrollwindow, ...)
              ggroup(container=container, horizontal=FALSE, ...)
            iggroup <- function(container, horizontal, expand, anchor, ...)
              ggroup(container=container, horizontal=horizontal, ...)

            
##             theArgs = list(...)
##             groupArgs = list()
##             for(i in c("spacing","use.scrollwindow")) {
##               if(!is.null(theArgs[[i]])) {
##                 groupArgs[[i]] = theArgs[[i]]
##                 theArgs[[i]] = NULL
##               }
##             }

##             theArgs$horizontal = FALSE
##             theArgs$container = container

            cg = oggroup(container, ...)
#            cg = do.call("ggroup",theArgs)

            labelGroup = ggroup(horizontal=TRUE, container=cg)

            rightArrow = system.file("images","1rightarrow.gif",package="gWidgets")
            downArrow = system.file("images","1downarrow.gif",package="gWidgets")

            icon = gimage(downArrow,container=labelGroup)
            label = glabel(text, container=labelGroup)

            ## we need this so that getBlock doesn't find cg's block
            eg1 = ggroup(container=cg, expand=TRUE,horizontal)

##             groupArgs$container=eg1
##             groupArgs$expand=TRUE
##            eg = do.call("ggroup", groupArgs)
            eg <- iggroup(container=eg1, horizontal=horizontal, ...)

#            obj = new("gExpandgrouptcltk",block = eg1, widget = eg,

            obj = new("gExpandgrouptcltk",block = cg, widget = eg, horizontal=horizontal,
              toolkit = toolkit, ID = getNewID(), e = new.env())


            tag(obj, "containerGroup") <- cg
            tag(obj, "expandGroup") <- eg1
            tag(obj, "icon") <- icon
            tag(obj, "label") <- label
            tag(obj, "state") <- FALSE
            tag(obj, "rightArrow") <- rightArrow
            tag(obj, "downArrow") <- downArrow
            tag(obj, "height") <- tkcget(getWidget(obj), "-height")
              
            changeState = function(h,...) {
              if((state <- tag(obj,"state"))) {
                visible(obj) <- FALSE
              } else {
                visible(obj) <- TRUE
              }
            }

            addHandlerClicked(icon, handler=changeState)
            addHandlerClicked(label, handler=changeState)


            visible(obj) <- FALSE       # initial state
            
            ## must take care of closing/opening
            if(!is.null(handler)) {
              addHandlerChanged(obj, handler=handler, action=action)
            }


            invisible(obj)
          })



## methods
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gExpandgrouptcltk", value = "gWidgettcltk"),
          function(obj, toolkit, value,  ...) {
            ## add value to expandgroup
            add(obj@widget, value, ...)
          })
setMethod(".addSpace",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gExpandgrouptcltk"),
          function(obj, toolkit, value,  ...) {
            ## add value to expandgroup
            addSpace(obj@widget, value, ...)
          })
setMethod(".addSpring",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gExpandgrouptcltk"),
          function(obj, toolkit,  ...) {
            ## add value to expandgroup
            addSpring(obj@widget, ...)
          })
## push onto label
setReplaceMethod(".font",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gExpandgrouptcltk"),
          function(obj, toolkit,  ..., value) {
            ## add value to expandgroup
            font(tag(obj,"label")) <- value
            return(obj)
          })


## Should make
## a) svalure refer to padding, ala ggroup padding
## b) names refer to label
## c) font refer to font of label
## d) visible refer to state

## value refers to padding
## FOr svalue<- we still accept non-numeric for setting lable
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gExpandgrouptcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            svalue(tag(obj,"label"))
          })

setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gExpandgrouptcltk",
                           value = "numeric"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   svalue(obj@widget, value)
                   return(obj)
                 })
## set name, but is deprecated
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gExpandgrouptcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   gwCat("Using names<- to set label value")
                   svalue(tag(obj,"label")) <- as.character(value)
                   return(obj)
                 })


## visible method
setMethod(".visible",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gExpandgrouptcltk"),
          function(obj, toolkit, set=TRUE,...) {
            tag(obj,"state")
          })

## control expand/close with logical
setReplaceMethod(".visible",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gExpandgrouptcltk"),
                 function(obj, toolkit, ..., value) {
                   W <- getWidget(obj)
                   ## cg = tag(obj,"containerGroup")
                   ## eg = tag(obj,"expandGroup")
                   if( (value <- as.logical(value)) ) {
                    ## true, expand
                     ##                     add(cg, eg, expand=TRUE)
                     tkpack("propagate", W, TRUE)
                     tkconfigure(W, height=tag(obj, "height"))
                     svalue(tag(obj,"icon")) <- tag(obj,"downArrow")
                   } else {
                     ##                     delete(cg,eg)
                     tag(obj, "height") <- tkwinfo("height", W)
                     tkpack("propagate", W, FALSE)
                     tkconfigure(W, height=1)
                     svalue(tag(obj,"icon")) <- tag(obj,"rightArrow")
                   }
                   tag(obj,"state") <-value

                   return(obj)
                 })


## names refers to label
setMethod(".names",
          signature(toolkit="guiWidgetsToolkittcltk",x="gExpandgrouptcltk"),
          function(x, toolkit) {
            svalue(tag(x,"label"))
          })

setReplaceMethod(".names",
                 signature(toolkit="guiWidgetsToolkittcltk",
                           x="gExpandgrouptcltk"),
                 function(x, toolkit, value) {
                   svalue(tag(x,"label")) <- as.character(value)
                   return(x)
                 })


setReplaceMethod(".font",
                 signature(toolkit="guiWidgetsToolkittcltk",
                           obj="gExpandgrouptcltk"),
                 function(obj, toolkit, value) {
                   font(tag(obj, "label")) <- value
                   return(obj)
                 })


## handlers
## putonto expander button
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gExpandgrouptcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addHandlerChanged(tag(obj,"icon"), handler, action,...)
            addHandlerChanged(tag(obj,"label"), handler, action,...)
          })
