setClass("gFrametcltk",
         contains="gGrouptcltk",
         prototype=prototype(new("gGrouptcltk"))
         )

## add a frame for packing. subclass of gGroup
setMethod(".gframe",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   text = "", markup=FALSE,
                   pos = 0, ## pos in [0,1] 0 for left, (.01,.99) center, 1 for right
                   horizontal=TRUE,
                   container=NULL,
                   ...) {

            force(toolkit)
            

            ## we can't do any markup here. Font() could be used
            if(markup) {
              gwCat(gettext("HTML markup not supported for title. \n"))
              text = gsub("<[^>]*>","",text)    # strip off HTML
            }

            ## where to put
            labAnchor = "nw"
            if(.33 < pos  && pos < .66)
              labAnchor = "n"
            else if(.66 <= pos)
              labAnchor = "ne"

            
            theArgs <- list(...)

            if(is(container,"logical") && container)
              container = gwindow()
            if(!is(container,"guiWidget")) {
              warning("Container is not correct. No NULL containers possible\n" )
              return()
            }

            tt <- getWidget(container)
            f <- tkwidget(tt, "ttk::labelframe", text=text, labelanchor=labAnchor)

            ## put in some padding. Adjust with svalue
            if(!is.null(theArgs$spacing))
              spacing <- theArgs$spacing
            else
              spacing <- 5
            tkconfigure(f,"padding"=spacing)

            ## XXX -- not sure this is supposed to be here
            ## ## handle expand and anchor arguments for packing frame
            ## argList = list(f)
            
            ## if(!is.null(theArgs$expand) && theArgs$expand) {
            ##   argList$expand = TRUE
            ##   argList$fill = "both"
            ## }
            ## if(is.null(theArgs$anchor))
            ##   theArgs$anchor= c(-1,1)
            
            ## argList$anchor = xyToAnchor(theArgs$anchor)

            ## do.call("tkpack",argList)

            
            obj = new("gFrametcltk",
              block=f, widget=f, toolkit=toolkit,
              horizontal=horizontal,
              ID=getNewID(), e = new.env())

            tag(obj,"title") <- text

            ## attach to container if there
            if(!is.null(container)) {
              add(container, obj,...)
            }

            return(obj)
          })

## methods
## should be same as from ggroup:
## svalue<- padding
## names<- name

## sets the padding. Same as ggroup
setReplaceMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gFrametcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ..., value) {
            ## adds some breathing room to object
            ## value is pixels
            widget <- getWidget(obj)
            tcl(widget,"configure","padding"=as.numeric(value))
            return(obj)
          })

## should put in a names argument to change label value
## return label name
setMethod(".names",signature(toolkit="guiWidgetsToolkittcltk",
                             x="gFrametcltk"),
          function(x, toolkit) {
            tag(x,"title")
          })


setReplaceMethod(".names",
                 signature(toolkit="guiWidgetsToolkittcltk",x = "gFrametcltk"),
                 function(x,toolkit,value) {

                   f <- x@widget
                   ## XXX What to put here?
                   tkconfigure(f,text=as.character(value))
                   tag(x,"title") <- value
                   return(x)
                 })

