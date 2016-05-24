## Constants
setBufferFonts <- function(textview, font.attr) {
  font.attr <- unlist(font.attr)
  nms <- names(font.attr)   # a vector -- not alist
  out <- ""
  if("style" %in% nms)
    out <- paste(out, toupperFirst(font.attr['style']))
  if("weight" %in% nms)
    out <- paste(out, toupperFirst(font.attr['weight']), sep=" ")
  if("size" %in% nms) {
    sz <- fontSizes[font.attr['size']]
    sz <- ceiling(12*sz)    # to font size
  } else {
    sz <- 12
  }
  out <- paste(out, sz, sep=" ")
  font <- pangoFontDescriptionFromString(out)
  textview$modifyFont(font)
  ## now for color
  if("color" %in% nms) {
    color <- font.attr['color']
    textview$modifyText(GtkStateType['normal'], color)
  }
}



setClass("gTextRGtk",
#         representation(tags="list"),
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

setMethod(".gtext",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   text=NULL,
                   width=NULL, height=300,
                   font.attr = NULL, wrap = TRUE,
                   handler = NULL, action=NULL,
                   container=NULL, ...) {

            force(toolkit)
            

            ## make textview
            textview = gtkTextViewNew()
            textview$SetLeftMargin(10)
            textview$SetRightMargin(10)
            if(wrap)
              textview$SetWrapMode(GtkWrapMode['word'])
            else
              textview$SetWrapMode(GtkWrapMode['none'])
            
            ## pack in a scrollwindow
            sw = gtkScrolledWindowNew()
#            group = ggroup()
#            add(group, sw, expand=TRUE)
            sw$SetPolicy("GTK_POLICY_AUTOMATIC","GTK_POLICY_AUTOMATIC")
            if(!is.null(width))
              sw$SetSizeRequest(width,height)


            sw$Add(textview)
            textview$Show()
            

#            obj = new("gTextRGtk", block=group, widget=textview, tags=tags, toolkit=toolkit)
#            obj = new("gTextRGtk", block=sw, widget=textview, tags=tags, toolkit=toolkit)

            obj <- as.gWidgetsRGtk2(textview, block=sw)
            
            ##   ## Handle attributes
            ##   if(!is.null(font.attr))
            ##     font(obj) <- font.attr


            ## font.attr specifies text properties for the entire buffer (gWidgets 0.0-39)
            if(!is.null(font.attr)) {
              .font(textview, toolkit) <- font.attr
#              setBufferFonts(textview, font.attr)
            }
            
            if(!is.null(text)) {
              add(obj, text, do.newline=FALSE)
            }
            
  
            ## attach to container
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE)
              add(container, obj,...)
            }
            
            if (!is.null(handler)) {
              id = addhandlerkeystroke(obj, handler, action)
            }
            return(obj)
          })


as.gWidgetsRGtk2.GtkTextView <- function(widget, ...) {

  theArgs <- list(...)
  if(!is.null(theArgs$block))
    block <- theArgs$block
  else
    block <- widget
  
  obj <- new("gTextRGtk", block=block, widget=widget,
             toolkit=guiToolkit("RGtk2"))

  ## ## add tags if not there
  ## if(is.null(tag(obj,"tags"))) {
  ##   buffer <- widget$GetBuffer()
  ##   tags <- .addTags(buffer)
  ##   tag(obj,"tags") <- tags
  ## }
  
  return(obj)
}

## add tags to buffer
## return tags
.addTags <- function(buffer) {

  ## weights
  fontWeights = names(PangoWeight)
  fontWeights = fontWeights[fontWeights != "normal"] # also in Styles

  tagtbl <- buffer$getTagTable()
  
  for(i in fontWeights)  
    if(is.null(tagtbl$lookup(i)))
      buffer$createTag(i, weight = PangoWeight[i])
  
  ## styles
  fontStyles = names(PangoStyle)
  for(i in fontStyles)
    if(is.null(tagtbl$lookup(i)))
      buffer$createTag(i, style = PangoStyle[i])
  ## family
  buffer$createTag("monospace",family = "monospace")

  for(i in names(fontSizes))
    if(is.null(tagtbl$lookup(i)))
      buffer$createTag(i, scale = fontSizes[i])
  
  ## colors -- 
  ##             fontColors = c("black","blue","red","yellow","brown","green","pink")
  ##             for(i in fontColors) {
  ##               buffer$createTag(i,foreground = i)
  ##               buffer$createTag(Paste(i,".background"),background = i)
  ##             }
  fontColors <-  colors()
  lapply(colors(), function(i) {
    if(is.null(tagtbl$lookup(i))) 
      buffer$createTag(i,foreground = i)
    if(is.null(tagtbl$lookup(Paste(i,".background"))))
    buffer$createTag(Paste(i,".background"),background = i)
  })
  
  
  
  tags = list(
    styles = fontStyles,
    family = "monospace",
    weights = fontWeights,
    sizes = names(fontSizes),
    foreground.colors = fontColors,
    background.colors = paste(fontColors,".background", sep="")
    )

  return(tags)
}


### methods

## drop=TRUE to get only mouse selected text
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTextRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            ## grab all text
            buffer = obj@widget$GetBuffer()
            if(is.null(drop) || drop == FALSE) {
              start = buffer$GetStartIter()$iter
              end = buffer$GetEndIter()$iter
            } else {
              ## return only **selected** text
              ## if drop==TRUE
              bounds = buffer$GetSelectionBounds()
              if(bounds$retval == FALSE) return("")
              start = bounds$start
              end = bounds$end
            }
            val <- buffer$GetText(start,end)
            return(val)
            })
          
##  svalue<-() replaces text
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTextRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   textbuffer = obj@widget$GetBuffer()
                   if(length(value) > 1)
                     value = paste(value, collapse="\n")
                   textbuffer$SetText(value)
                   return(obj)
                 })


## clear all text in buffer
setMethod("dispose",signature(obj="gTextRGtk"),
          function(obj,...)  {
            .dispose(obj, obj@toolkit, ...)
          })
setMethod(".dispose",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTextRGtk"),
          function(obj, toolkit,  ...) {
            buffer = obj@widget$GetBuffer()
            startiter =  buffer$GetStartIter()$iter
            enditer =  buffer$GetEndIter()$iter
            buffer$Delete(startiter, enditer)
          })



### Add method is a workhorse for this class. Value can be
## * a line of text
## * a vector of lines of text
## * an gWidget instance
## need to do where value of "point"
## add, as a method, needs to have a consistent signature. I'

## add text
setMethod(".insert",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj = "gTextRGtk"),
          function(obj, toolkit, value, where = c("end","beginning","at.cursor"),
                   font.attr = NULL,
                   do.newline = TRUE, ...) {
            ## just call add
            where = match.arg(where)
            .add(obj, toolkit, value, where=where, font.attr=font.attr,
                 do.newline=do.newline, ...)
          })
## add does all the work
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTextRGtk",value="character"),
          function(obj, toolkit, value,  ...) {
            theArgs = list(...)                      # look for font.attr, do.newline, where

            do.newline = ifelse(is.null(theArgs$do.newline), TRUE, as.logical(theArgs$do.newline))
            markup = theArgs$font.attr

            if(!is.null(markup)) {
              if(is.null(tag(obj, "tags")))
                tag(obj, "tags") <- .addTags(obj@widget$getBuffer())
              markup = markup[markup %in% unlist(tag(obj,"tags"))] # only some markup
            }
            where <- getWithDefault(theArgs$where, "end")

            
            view <- obj@widget
            buffer = view$GetBuffer()
            iter = switch(where,
              "end"=buffer$GetEndIter()$iter,
              "beginning"=buffer$GetStartIter()$iter,
              {gwCat("Only end, beginning implemented")
               buffer$GetEndIter()$iter
             })
            
            for(i in 1:length(value) ) {
              if(is.null(markup)) {
                buffer$Insert(iter, value[i])
              } else {
                lst = list(object=buffer, iter=iter, text=value[i])
                for(key in names(markup)) {
                  if(is.list(markup))
                    lst[[key]] <- markup[[key]]
                  else
                    lst[[key]] <- markup[key]
                }
                do.call("gtkTextBufferInsertWithTagsByName",lst)
              }
              if(do.newline) buffer$Insert(iter,"\n")
            }

            ## scroll to end -- if appended to end
            if(where == "end") {
              gdkWindowProcessAllUpdates()
              while (gtkEventsPending())
                gtkMainIterationDo(blocking=FALSE)

             end <- buffer$getEndIter()$iter
              view$scrollToIter(end, within.margin = 0,
                                use.align=TRUE)
            }
          })

## add a widget
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTextRGtk",value="guiWidget"),
          function(obj, toolkit, value,  ...) {
            .add(obj,toolkit, value@widget, ...)
          })

setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTextRGtk",value="gWidgetRGtk"),
          function(obj, toolkit, value,  ...) {
            theArgs = list(...)                      # look for font.attr, do.newline, where
            
            do.newline = ifelse(is.null(theArgs$do.newline), TRUE, as.logical(theArgs$do.newline))

            where = ifelse(is.null(theArgs$where), "end",theArgs$where)
            buffer = obj@widget$GetBuffer()
            iter = switch(where,
              "end"=buffer$GetEndIter()$iter,
              "beginning"=buffer$GetStartIter()$iter,
              {gwCat("Only end, beginning implemented")
               buffer$GetEndIter()$iter
             })
            

            anchor = buffer$CreateChildAnchor(iter)
            getWidget(obj)$AddChildAtAnchor(getWidget(value), anchor)
            if(do.newline) buffer$Insert(iter,"\n")
            
            })


## set the font for the selected area of the gtext object
setReplaceMethod(".font",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTextRGtk"),
                 function(obj, toolkit, ..., value) {

                   textview <- getWidget(obj)
                   buffer <- textview$buffer
                   tagtbl <- buffer$getTagTable()

                   ## get bounds.
                   bounds <- buffer$GetSelectionBounds()
                   if(bounds$retval == FALSE) {
                     ## if no text selected, we set for entire buffer
                     ## change entire buffer -- new as of 0.64
                     start <- buffer$GetStartIter()$iter
                     end <- buffer$GetEndIter()$iter
                     buffer$removeAllTags(start, end)
                     ## now set font by calling again
                   } else {
                     start <- bounds$start
                     end <- bounds$end
                   }

                   tags <- sapply(value, function(i) i, simplify=FALSE) # to list

                   
                   ## we have family, style, weight, size, color
                   weights <- RGtk2::PangoWeight
                   if(!is.null(wt <- tags$weight) && wt %in% names(weights)) {
                     if(is.null(tagtbl$lookup(wt)))
                       buffer$createTag(wt, weight=weights[wt])
                     buffer$ApplyTagByName(wt, start, end)
                   }

                   
                   ## style
                   styles <- RGtk2::PangoStyle
                   if(!is.null(style <- tags$style) && style %in% names(styles)) {
                     if(is.null(tagtbl$lookup(style)))
                       buffer$createTag(style, style=styles[style])
                     buffer$ApplyTagByName(style, start, end)
                   }
                   
                   ## family
                   families <- c("normal","sans", "serif", "monospace")
                   if(!is.null(family <- tags$family) && family %in% families) {
                     if(is.null(tagtbl$lookup(family)))
                       buffer$createTag(family, family=family)
                     buffer$ApplyTagByName(family, start, end)
                   }

                   ## size
                   ## Pango Scale for converting between name and numeric value
                   
                   if(!is.null(size <- tags$size)) {
                     if(is.character(size)) 
                       size <- fontSizes[size]
                     else
                       size <- size/12
                     if(is.null(tagtbl$lookup(size)))
                       buffer$createTag(size, scale=size)
                     buffer$ApplyTagByName(size, start, end)
                   }
                   

                   
                   
                   ## color
                   if(!is.null(color <- tags$color) && color %in% colors()) {
                     if(is.null(tagtbl$lookup(color)))
                       buffer$createTag(color, foreground=color)
                     ## do we need to remove colors?
                     sapply(colors(), function(i) {
                       if(!is.null(tagtbl$lookup(i)))
                         buffer$RemoveTagByName(i, start, end)
                     })
                     buffer$ApplyTagByName(color, start, end)
                   }
                   
                   ## how to modify background color?
                     ## bg <- sprintf("%s.background", color)
                     ## if(is.null(tagtbl$lookup(bg)))
                     ##   buffer$createTag(bg, background=bg)

                       
                       



                     
                     ## get tags that are known
                     ## tags = value
                     ## tags = tags[tags %in% unlist(tag(obj,"tags"))]
                     ## if(length(tags) == 0) {
                     ##   cat(gettext("Invalid font specification\n"))
                     ##   return(obj)
                     ## }
                   

                     ## for(i in tags) {
                     ##   ## color is special
                     ##   if(length(names(i)) && names(i)[1] == "color") {
                     ##     sapply(colors(), function(j)
                     ##            buffer$RemoveTagByName(j, bounds$start, bounds$end))
                     ##   }
                       
                     ## buffer$ApplyTagByName(i, bounds$start, bounds$end)
                     ## }
                   
                   return(obj)
                 })



### handlers
setMethod(".addhandlerkeystroke",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTextRGtk"),
          function(obj,toolkit, handler=NULL, action=NULL,...) {
            widget <- getWidget(obj)
            ID <-
              gSignalConnect(widget,signal = "key-release-event", # or key-press-event
                             f = function(d,widget,event,...) {
                               h <- list(obj=d$obj,action=d$action)
                               key <- event$GetString()
                               h$key <- key
                               ## for modifiers
                               state <- event$getState()
                               if(state == 0)
                                 modifier <- NA
                               else 
                                 modifier <- gsub("-mask$", "",names(GdkModifierType)[GdkModifierType == state])
                               h$modifier <- modifier
                               
                               if(!is.null(d$handler) &&
                                  is.function(d$handler))
                                 d$handler(h,...)
                               return(FALSE) # propogate
                             },
                             user.data.first = TRUE,
                             data = list(obj=obj,handler=handler, action=action)
                             )
            invisible(ID)
          })


## generic
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTextRGtk"),
          function(obj,toolkit, handler=NULL, action=NULL,...) {
            .addhandlerkeystroke(obj,toolkit,handler,action,...)
          })

##################################################
##################################################
## testing This is copied from pygtk tutorial
addWidgetAtPoint = function(obj, value) {
  evb = gtkEventBoxNew()
  evb$SetVisibleWindow(FALSE)
  evb$SetBorderWidth(15)
  evb$AddEvents(c(GdkEventMask["button-press-mask"],
                  GdkEventMask["button-release-mask"],
                  GdkEventMask["button-motion-mask"],
                  GdkEventMask["button-motion-hint-mask"]))
  
  widget = value@widget
  if(is(widget,"gContainer") || is(widget,"gComponent"))
    widget = widget@widget                 # for instance, ggroup
  
  evb$Add(widget)
  evb$ShowAll()

  ## connect move handler?
  gtktry(connectSignal(evb,
                signal = "button-press-event",
                f = movableWidgetButtonPressHandler,
                data = list(obj=obj@widget),
                user.data.first = TRUE),
      silent=TRUE)
  gtktry(connectSignal(evb,
                signal = "button-release-event",
                f = movableWidgetButtonReleaseHandler,
                data = list(obj=obj@widget),
                user.data.first = TRUE),
      silent=TRUE)
  gtktry(connectSignal(evb,
                signal = "motion-notify-event",
                f = movableWidgetButtonMotionNotifyHandler,
                data = list(obj=obj@widget),
                user.data.first = TRUE),
      silent=TRUE)
  
  
  ## get xpos, ypos
  ptr = obj@widget$GetPointer()
  xpos = ptr$x; ypos = ptr$y
  xpos = 1; ypos = 1

  buffer = obj@widget$GetBuffer()
  iter = buffer$GetEndIter()$iter
  anchor = buffer$CreateChildAnchor(iter)
  obj@widget$AddChildAtAnchor(evb, anchor)

  return()
  obj@widget$AddChildInWindow(evb, GtkTextWindowType['widget'],
               xpos, ypos)                 

}

movableWidgetButtonPressHandler = function(h, widget, event, ...) {

  textview = h$obj
  info = widget$GetData("moveable-widget-data")

  if(is.null(info)) {
    info = list("start_x"= NA, "start_y"=NA, button = NA)
    widget$SetData("moveable-widget-data", info)
  }

  if(!is.list(info[['button']]) || is.na(info[['button']])) {
    info$button = event
    allocation = widget$GetAllocation()
    info[['start_x']] = allocation$x
    info[['start_y']] = allocation$y
    info[['click_x']] = allocation$x + event$GetX()
    info[['click_y']] = allocation$y + event$GetY()
    widget$SetData("moveable-widget-data", info)
  }
  return(FALSE)
}

movableWidgetButtonReleaseHandler = function(h, widget, event, ...) {
  info = widget$GetData("moveable-widget-data")
  if(!is.list(info[['button']]) || is.na(info[['button']])) {
    gwCat("relase handler failed\n")
    return(FALSE)
  }
  
  info = widget$GetData("moveable-widget-data")
  
  x = info[['start_x']] + event$GetX() + widget$GetAllocation()$x - info[['click_x']]
  y = info[['start_y']] + event$GetY() + widget$GetAllocation()$y - info[['click_y']]
  
  widget$SetData("moveable-widget-data", NULL)
  
  h$obj$MoveChild(widget, x,y) 
  
  return(FALSE)
}

movableWidgetButtonMotionNotifyHandler = function(h, widget, event, ...) {
  
  info = widget$GetData("moveable-widget-data")
  
  ptr = widget$GetPointer()
  allocation = widget$GetAllocation()
  x = ptr$x + allocation$x
  y = ptr$y + allocation$y
  
  x = info[['start_x']] + (x - info[['click_x']])
  y = info[['start_y']] + (y - info[['click_y']])
  
  h$obj$MoveChild(widget, x,y)   
  return(FALSE)
}
