## TODO
## * FONTS
## some common function

## does gtext object have a selection
hasSelection = function(obj) tclvalue(tktag.ranges(getWidget(obj),"sel")) != ""

## begin
setClass("gTexttcltk",
           representation(tags="list"),
           contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )

setMethod(".gtext",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   text=NULL,
                   width=NULL, height=200,
                   font.attr = NULL, wrap = TRUE,
                   handler = NULL, action=NULL,
                   container=NULL, ...) {


            force(toolkit)
            theArgs <- list(...)
            ladd <- function(...,bg) add(...) # don't pass bg to add -- if present

            if(is(container,"logical") && container)
              container <- gwindow()
            if(!is(container,"guiWidget")) {
              warning("Container is not correct. No NULL containers possible\n" )
              return()
            }

            ## options
            ## wrap
            if(wrap) wrap <- "word" else wrap <- "none"
            ## background color
            bg <- if(!is.null(theArgs$bg)) theArgs$bg else "white"
            
            tt <- getWidget(container)
            gp <- ttkframe(tt)
            
            xscr <- ttkscrollbar(gp, orient="horizontal",
                                 command=function(...)tkxview(txt,...))
            yscr <- ttkscrollbar(gp, 
                                 command=function(...)tkyview(txt,...))
            
            txt <- tktext(gp,
                          bg=bg, setgrid=FALSE, #font="courier",
                          undo = TRUE,                  # undo support
                          xscrollcommand=function(...)tkset(xscr,...),
                          yscrollcommand=function(...)tkset(yscr,...),
                          wrap=wrap)
            
            ## pack into a grid
            ## see tkFAQ 10.1 -- makes for automatic resizing
            tkgrid(txt,row=0,column=0, sticky="news")
            tkgrid(yscr,row=0,column=1, sticky="ns")
            tkgrid(xscr, row=1, column=0, sticky="ew")
            tkgrid.columnconfigure(gp, 0, weight=1)
            tkgrid.rowconfigure(gp, 0, weight=1)

            ## from tcltk2 package, this package is installed

            if(windowingsystem() != "aqua") {
              tcl("autoscroll::autoscroll", xscr)
              tcl("autoscroll::autoscroll", yscr)
            }

            
            ## set point
            tkmark.set(txt,"insert","0.0")
            
            obj <- new("gTexttcltk", block=gp, widget=txt, tags=list(),
              toolkit=toolkit,ID=getNewID(), e = new.env())

            ## font.attr sets text properties for entire buffer
            if(!is.null(font.attr)) {
              font(obj) <- font.attr
            }


            ## add initial text
            if(!is.null(text)) {
              add(obj, text)
            }

            ## set height if requested
            if(!is.null(width))
              ## width height in terms of characters
              size(obj) <- c(width,height)
            
##            adddropsource(obj)
            adddroptarget(obj)

            
            ## attach to container
            add(container, obj,...)

            ## add handler
            if (!is.null(handler)) {
              obj@e$handler.id <- addhandlerkeystroke(obj, handler, action)
            }
            return(obj)
          })

## drop=TRUE to get only mouse selected text
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gTexttcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {

            ## rongui request, if INDEX = TRUE return selected text
            ## by index in the buffer
            if(!is.null(index) && index == TRUE) {
              ## get the selected text from gtext,
              ## return the index instead of text.
              if(hasSelection(obj))
                ## row.column: row 1-based, column 0-based
                val <- as.character(tktag.ranges(getWidget(obj),"sel"))
              else
                val <- c(0,0)
              return(as.numeric(val))
            }

            ## otherwise we return text
            ## if drop=FALSE or NULL grab all text
            ## if drop=TRUE, get selected text only
            if(is.null(drop) || drop == FALSE) {
              val <-  tclvalue(tkget(getWidget(obj),"0.0","end"))
              ## strip off last "\n"'s
              val <- gsub("\n*$","",val)
            } else {
              range <- as.numeric(tktag.ranges(getWidget(obj),"sel"))
              ## range is numeric(0) if none
              if(length(range) > 0)
                val <- tclvalue(tkget(getWidget(obj),"sel.first","sel.last"))
              else
                val <- ""
            }
                                        ## val = unlist(strsplit(val,"\n"))
            return(val)
          })

##  svalue<-() replaces text
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gTexttcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {

                   ## how to clear out old text
                   tkdelete(getWidget(obj),"0.0","end")

                   if(length(value) > 1)
                     value <- paste(value, collapse="\n")
                   tkinsert(getWidget(obj),"end",value)

                   tksee(getWidget(obj),"0.0")
                   
                   return(obj)
                 })


## clear all text in buffer
setMethod("dispose",signature(obj="gTexttcltk"),
          function(obj,...)  {
            .dispose(obj, obj@toolkit, ...)
          })
setMethod(".dispose",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gTexttcltk"),
          function(obj, toolkit,  ...) {
            svalue(obj) <- ""
          })


### insert (was add) method is a workhorse for this class. Value can be
## * a line of text
## * a vector of lines of text
## need to do where value of "at.cursor"

## add text
setMethod(".insert",
          signature(toolkit="guiWidgetsToolkittcltk",obj = "gTexttcltk"),
          function(obj, toolkit, value, where = c("end","beginning","at.cursor"),
                   font.attr = NULL,
                   do.newline = TRUE, ...) {
            ## just call add
            where <- match.arg(where)
            .add(obj, toolkit, value, where=where, font.attr=font.attr,
                 do.newline=do.newline, ...)
          })

## should be .insert, but legacy
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gTexttcltk",value="character"),
          function(obj, toolkit, value,  ...) {
            theArgs <- list(...)                      # look for font.attr, do.newline, where

            do.newline <- ifelse(is.null(theArgs$do.newline), TRUE, as.logical(theArgs$do.newline))



            
            where <- ifelse(is.null(theArgs$where), "end", theArgs$where)
            where <- switch(where,
                            "at.cursor"="insert",
                            "beginning"="0.0",
                            "end")
                           

            txt <- getWidget(obj)
            
            value <- paste(value,collapse="\n")
            if(do.newline)
              value = paste(value,"\n",sep="")

            ### Handle markup here
            markup <- theArgs$font.attr
            
            if(!is.null(markup)) {
              ## bit of a hack to set font
              fname <- paste(as.character(date()),rnorm(1), sep="") ## some random string
              fontList <- fontlistFromMarkup(markup, fname)
              do.call("tkfont.create", fontList)

              tkmark.set(txt, "left","insert"); tkmark.gravity(txt,"left","left")
              tkmark.set(txt, "right","insert"); tkmark.gravity(txt,"right","right")
              tkinsert(txt, where, value)
              tktag.add(txt, fname, "left","right")
              tktag.configure(txt, fname, font=fname)

              if("color" %in% names(markup))
                tktag.configure(txt, fname, foreground=markup['color'])
            } else {
              ## no markup
              tkinsert(getWidget(obj),where,value)
            }

            ## does this place the cursor? TK FAQ 10.6
            tksee(getWidget(obj),where) # where = "end" or "0.0"
            
          })

## add a widget
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gTexttcltk",value="guiWidget"),
          function(obj, toolkit, value,  ...) {
            .add(obj,toolkit, value@widget, ...)
          })

setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gTexttcltk",value="gWidgettcltk"),
          function(obj, toolkit, value,  ...) {

            message("gtext: implement adding a widget to text area\n")
            return()
            })


## set the font for the selected area of the gtext object
setReplaceMethod(".font",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gTexttcltk"),
                 function(obj, toolkit, ..., value) {
                   ## if a selection, set font, else set font for buffer
                   widget <- getWidget(obj)
                   if(hasSelection(obj)) {
                     selected <- as.character(tktag.ranges(getWidget(obj),"sel"))
                     fname <- paste(as.character(date()),rnorm(1), sep="") ## some random string
                     ## make font, tag in buffer, configure tag
                     fontList <- fontlistFromMarkup(value)
                     do.call("tkfont.create", merge(list(fname), fontList))
                     tktag.add(widget, fname, selected[1], selected[2])
                     tktag.configure(widget, fname, font=fname)
                     if("color" %in% names(value))
                       tktag.configure(widget, fname, foreground=value['color'])
                   } else {
                     ## clear out old tags -- we are resetting
                     tagNames <- as.character(tktag.names(widget))
                     sapply(tagNames, function(i) tktag.delete(widget, i))
                     .font(widget, toolkit, ...) <- value
                   }
                   return(obj)
                 })



setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gTexttcltk"),
          function(obj,toolkit, handler=NULL, action=NULL,...) {
            .addhandlerkeystroke(obj,toolkit,handler,action)
          })

