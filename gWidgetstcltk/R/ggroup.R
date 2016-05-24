## class in aaaClasses.R
## constructor
setMethod(".ggroup",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   horizontal = TRUE, spacing = 5,
                   use.scrollwindow = FALSE, 
                   container = NULL, ... 
                   ) {

            force(toolkit)
            
            theArgs = list(...)                   # raise.on.dragmotion
            
            if(is.null(spacing))
              spacing = 0


            if(is.null(container)) {
              message(gettext("No NULL containers in tcltk. Creating a new window\n"))
              container=gwindow()
            } else if(is.logical(container) && container) {
              container = gwindow()
            }

            if(!is(container,"guiWidget")) {
              container = gwindow()
            }

            tt <- getWidget(container)

            ## implement scrollbars if asked.
            ## XXX Not quite working as desired...
            if(use.scrollwindow == TRUE && windowingsystem() != "aqua") {
              ## cf http://mail.python.org/pipermail/python-list/1999-June/005180.html
              block <- ttkframe(tt)

              ## put in one direction only
              widget <- tkcanvas(block)
              
              if(horizontal) {
                xscr <- ttkscrollbar(block, orient="horizontal",
                                   command=function(...)tkxview(widget,...))
                tkconfigure(widget, xscrollcommand = function(...) tkset(xscr,...))

                ## Pack into a grid
                ## see tkFAQ 10.1 -- makes for automatic resizing
                tkgrid(widget,row=0,column=0, sticky="news")
                tkgrid(xscr, row=1, column=0, sticky="ew")
                tkgrid.rowconfigure(block, 0, weight=1)
                
                tcl("autoscroll::autoscroll", xscr)

              } else {
                yscr <- ttkscrollbar(block, 
                                     command=function(...)tkyview(widget,...))
                tkconfigure(widget, yscrollcommand = function(...) tkset(yscr,...))
              
                ## Pack into a grid
                ## see tkFAQ 10.1 -- makes for automatic resizing
                tkgrid(widget,row=0,column=0, sticky="news")
                tkgrid(yscr,row=0,column=1, sticky="ns")
                tkgrid.columnconfigure(block, 0, weight=1)
                
                tcl("autoscroll::autoscroll", yscr)
              }
              
              ## Set up frame
              gp <- ttkframe(widget)
              gpID <- tcl(widget,"create","window",0,0,anchor="nw",window=gp)
              tkgrid.columnconfigure(widget,0,weight=1)
              tkgrid.rowconfigure(widget,0,weight=1)
              

                
              ## give an initial size
#              gpwidth <- getWithDefault(theArgs$width, 300)
#              gpheight <- getWithDefault(theArgs$height, 300)
#              if(horizontal)
#                tkitemconfigure(widget, gpID, height=gpheight)
#              else
#                tkitemconfigure(widget, gpID, width=gpwidth)
              
              tcl("update","idletasks")


              ## tkbind(widget,"<Configure>",function() {
              ##   bbox <- tcl(widget,"bbox","all")
              ##   tcl(widget,"config",scrollregion=bbox)
              ## })

              tkbind(block, "<Map>", function() {
                if(horizontal) {
                  width <- tkwinfo("width", block)
                  tkconfigure(widget, width=width)
                } else {
                  height <- tkwinfo("height", block)
                  tkconfigure(widget, height=height)
                }
              })

              tkbind(gp, "<Map>", function() {
                if(horizontal) {
                  tkconfigure(widget, height=tkwinfo("height", gp))
                } else {
                  tkconfigure(widget, width=tkwinfo("width", gp))
                }
              })
              
              tkbind(gp,"<Configure>",function() {
                bbox <- tcl(widget,"bbox","all")
                tcl(widget,"config",scrollregion=bbox)
              })

              
              tkbind(widget, "<Configure>", function(W) {
                width <- as.numeric(tkwinfo("width", W))
                height <- as.numeric(tkwinfo("height", W))
                
#                gpwidth <- as.numeric(tkwinfo("width", gp))
#                gpheight <- as.numeric(tkwinfo("height", gp))
                gpwidth <- as.numeric(tkwinfo("width", block))
                gpheight <- as.numeric(tkwinfo("height", block))
                
                if(gpwidth < width && !horizontal)
                  tkitemconfigure(widget, gpID, width=width)
                if(gpheight < height && horizontal)
                  tkitemconfigure(widget, gpID, height=height)
              })
            } else {
              gp <- ttkframe(tt)
              tkconfigure(gp, borderwidth=0) # XXX
              block <- gp
              widget <- NULL      # for later
            }

            tkconfigure(gp, padding=spacing)
            

            if(!is.null(theArgs$debug)) {
              theArgs$debug <- NULL
              tkconfigure(gp,borderwidth=4, relief="solid")
            }
            
            obj = new("gGrouptcltk",block=block, widget=gp,
              horizontal=horizontal,
              e = new.env(), ID=getNewID(), toolkit=toolkit
              )

            
            ## to move widget when scrolling
            ## if(!is.null(widget <- tag(value,"scrollable.widget"))) {
            ##  tkxview.moveto(widget,1)
            ##  tkyview.moveto(widget,1)
            ## }
            .tag(obj,toolkit, i="scrollable.widget") <- widget
            obj@e$i <- widget

            ## attach to container if there
            if(!is.null(container)) {
              theArgs$obj <- container
              theArgs$value <- obj
              do.call("add", theArgs)
            }

            ## raise if we drag across
            if(!is.null(theArgs$raise.on.dragmotion)) {
#              tkbind(gp, "<Motion>", function(W) {
#                tkfocus(W)
#              })
                     
            }
            return(obj)
          })


##################################################
## methods

setReplaceMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gGrouptcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ..., value) {
            ## adds some breathing room to object
            ## value is pixels
            gp <- getWidget(obj)
#            tkcofigure(gp,padx=value,pady=value)
            tkconfigure(gp,padding = value)            

            return(obj)
          })

##################################################
## handlers
