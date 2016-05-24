## TODO:
## * drag and drop onto tabs, raise on motion,
## * the [ method is not working
## * dispose is forgetting -- not hiding -- the child widget
## * where is my delete method?

setMethod(".gnotebook",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   tab.pos = 3,                          # same as pos= in text
                   closebuttons = FALSE,
                   dontCloseThese = NULL,                 # integer of tabs not to close
                   container=NULL,                           # add to this container
                   ...) {
            
            force(toolkit)

            if(is(container,"logical") && container)
              container = gwindow()
            if(!is(container,"guiWidget")) {
              warning("Container is not correct. No NULL containers possible\n" )
              return()
            }
            
            tt <- getWidget(container)
            gp <- ttkframe(tt)
            nb <- ttknotebook(gp)
            tkpack(nb, expand=TRUE, fill="both")                  # pack into gp, gp packed during add(.)
            
            ## tabpos
            if(tab.pos !=3)
              gwCat(gettext("tab.pos is not implemented\n"))

            ## create gnotebook object
            obj = new("gNotebooktcltk", block=gp, widget=nb,
              toolkit=toolkit,ID=getNewID(),e = new.env(),
              closebuttons = as.logical(closebuttons),
              dontCloseThese = ifelse(is.null(dontCloseThese),0,dontCloseThese))


            ## add to container
            add(container, obj, ...)

            
            invisible(obj)
          })

### methods
## return the current tab number
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gNotebooktcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {

            nb <- getWidget(obj)
            curTab <- tclvalue(tcl(nb,"select"))
            allTabs <- unlist(strsplit(tclvalue(tcl(nb,"tabs")),"\\s+"))
            which(curTab == allTabs) # 1-based
          })

## set the current tab to value
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gNotebooktcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {


                   nb <- getWidget(obj)
                   n = length(obj)

                   value <- max(1,min(value,n))
                   tcl(nb,"select",value - 1) # 0 -based
                   
                   return(obj)
                 })


## remove the current tab
## this should be called delete -- which is used to remove objects
setMethod(".dispose",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gNotebooktcltk"),
          function(obj, toolkit,  ...) {

            nb <- getWidget(obj)

            theArgs = list(...)
            to.right=ifelse(!is.null(theArgs$to.right), theArgs$to.right,FALSE)
            dontCloseThese = obj@dontCloseThese
            if(dontCloseThese == 0) dontCloseThese = NULL
            deleteOK = function(i) {
              if(is.null(dontCloseThese)) return(TRUE)
              if(i %in% dontCloseThese) return(FALSE)
              return(TRUE)
            }
            cur.pageno = svalue(obj)

            ## we clear out the current page unless there is more!
            inds = 0
            if(to.right) {
              n = length(obj)
              no.right = n - cur.pageno
              if(no.right > 0) 
                inds = no.right:0      # must work from last backwards
            }
            ## clear out would like "hide" here, as then we can
            ## readd. Not working here? why not?

            children <- obj@e$childComponents
            for(i in inds) {
              j = cur.pageno + i
              if(deleteOK(j)) {
#                tcl(nb,"hide", j - 1)
                tcl(nb,"forget", j - 1)
                children[[j]] <- NULL
              }
            }
            obj@e$childComponents <- children

            if(cur.pageno > 0 && length(children)) {        # error if no pages
              if(cur.pageno <= length(obj))
                svalue(obj) <- cur.pageno
              else
                svalue(obj) <- length(obj)
            }
          })


setMethod(".delete",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gNotebooktcltk", widget="gWidgettcltk"),
          function(obj, toolkit, widget, ...) {
            nb <- getWidget(obj)
            childWidget <- getBlock(widget)

            tcl(nb,"forget",childWidget)
            ## remove from childComponents
            children <- obj@e$childComponents
            ind <- sapply(children, function(i) digest(i) == digest(widget))
            if(any(ind))
              children <- children[!ind]
            obj@e$childComponents <- children
          })



### add() is a workhorse method here. Several args available in ...
#add.gNotebook = functionf(obj, value,
#  label="", markup = FALSE,
#  index = NULL, override.closebutton = FALSE, ...) {
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gNotebooktcltk",
                    value="guiWidget"),
          function(obj, toolkit, value,  ...) {
            .add(obj, toolkit, value@widget, ...)
          })

setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gNotebooktcltk",
                    value="gWidgettcltk"),
          function(obj, toolkit, value,  ...) {

            ## add parent, children
            childComponents <- obj@e$childComponents
            if(is.null(childComponents))
              childComponents <- list()
            obj@e$childComponents <- c(childComponents, value)
            value@e$parentContainer <- obj


            ## in ... we have many possibilies
            ## label -- for setting label  (also look for name)
            ## index for setting the index of page to add
            ## markup -- markup label
            ## override.closebutton -- to not put closebutton even if set in constructor

            nb <- getWidget(obj)
            widget <- getBlock(value)
            
            ## process ...
            theArgs = list(...)                      # for making generic
            if(!is.null(theArgs$label)) {
              label = theArgs$label
            } else if(!is.null(theArgs$name)) {
              label = theArgs$name
            } else {
              label = id(obj)
              if(is.null(label))
                label = "unnamed"
            }
            
            index = if (is.null(theArgs$index)) NULL else theArgs$index
            if(!is.null(theArgs$pageno)) index = theArgs$pageno # also called paegno
            markup = if (is.null(theArgs$markup)) FALSE  else theArgs$markup

            override.closebutton =
              if (is.null(theArgs$override.closebutton))
                FALSE
              else
                as.logical(theArgs$override.closebutton)


            packingOptions = list()
            packingOptions$sticky = if(is.null(theArgs$anchor))
              "nw"
            else
              xyToAnchor(theArgs$anchor)         
            
            if(!is.null(theArgs$expand) && theArgs$expand) 
              packingOptions$sticky="news"

            
            ## label -- a string in tcltk
            if(!is.character(label))
              label = svalue(label)
            
            
            ## closebutton
            ## closebutton
            file <- system.file("images/cancel.gif",package="gWidgets")
            closeb <- tcl("image","create","photo",file=file)
            
            doCloseButton <- FALSE
            
            if(!is.null(obj@closebuttons) &&
               as.logical(obj@closebuttons) &&
               !override.closebutton) {
              doCloseButton <- TRUE
              gwCat(gettext("gnotebook: close buttons are not active\n"))
            } 

            ## where
#            if(!is.null(index)) index = max(1,min(index,n))

            ## add drop motion for labels


            ##
            ## Can't do close buttons until we can identify when a tab click is on the icon
            ## We should be able to: nb identify element x y should work, as in
            ## tkbind(nb, "<Button-1>", function(x,y) {
            ##      tcl(nb, "identify", "element", x, y)
            ## })
            ## however this failes, only identify works and we can't sort out label from image
            ## if(doCloseButton) {
            ##   packingOptions$image=closeb
            ##   packingOptions$compound = "right"
            ## }
            
            if(is.null(index)) {
              f <- function(...)
                tcl(nb,"add", widget, text=label,...)
              do.call(f,packingOptions)
            } else {
              if(doCloseButton)
                f <- function(...) 
                  tcl(nb,"add",index-1, widget, text=label,...)
              do.call(f,packingOptions)
            }
          })
            
## Regular R methods treat gnotebook like a vector

## find out number of pages
setMethod(".length",
          signature(toolkit="guiWidgetsToolkittcltk",x="gNotebooktcltk"),
          function(x, toolkit) {
            nb <- getWidget(x)
            as.numeric(tclvalue(tcl(nb,"index","end")))
          })

## return tabnames
setMethod(".names",signature(toolkit="guiWidgetsToolkittcltk",x="gNotebooktcltk"),
          function(x, toolkit) {
            nb <- getWidget(x)
            n <- length(x)
            if(n > 0)
              vals <- sapply(1:n, function(i) tclvalue(tcl(nb,"tab",i - 1, "-text")))
            else
              vals <- NA
            return(vals)
          })

## can assigne with names(x) <-x or even names(x)[i] <- "single name"
setReplaceMethod(".names",
                 signature(toolkit="guiWidgetsToolkittcltk",x = "gNotebooktcltk"),
                 function(x,toolkit, value) {
                   nb <- getWidget(x)
                   n = length(x)
                   
                   if(length(value) != n)
                     stop(gettext("New names for notebook must have proper length"))
                   
                   lapply(1:n, function(i) {
                     tcl(nb,"tab",i-1, text=value[i])
                   })

                   return(x)
                 })


## return widget contained in notebook page i as a  list or single widget
setMethod("[",
          signature(x="gNotebooktcltk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gNotebooktcltk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            
            if(missing(i))
              i = 1:length(x)

            children <- x@e$childComponents[i]
            if(length(children) == 1)
              children <- children[[1]]

            return(children)
          })


## Puts widget into a position
setReplaceMethod("[",
                 signature(x="gNotebooktcltk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gNotebooktcltk"),
          function(x, toolkit, i, j, ..., value) {
            ##
            message(gettext("Can't add widget via [<-\n"))
            return()

            nb <- getWidget(x)
            widget <- getBlock(value)

            if(missing(i))
              stop(gettext("Missing value for i"))
            ## works, but parent is all messed up
            tcl(nb,"add",i - 1, value)
            
          })


### handlers
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gNotebooktcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addhandlerexpose(obj,toolkit, handler,action)
          })


setMethod(".addhandlerexpose",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gNotebooktcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addHandler(obj,toolkit,"<<NotebookTabChanged>>",handler, action)
          })

