## some dialogs for R
## dialogs don't get windows, they make them
## dialogs are modal
## dialogs return their value -- not an object. so source(gfile()) should work

setMethod(".gmessage",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   message,
                   title = "message",
                   icon = c("info","warning","error","question"),
                   parent=NULL,
                   handler = NULL,
                   action = NULL,
                   ...
                   ) {
            
            force(toolkit)

            icon = match.arg(icon)
            
            icon = Paste("GTK_MESSAGE_",toupper(match.arg(icon)))
            button = "GTK_BUTTONS_OK"

            ## parent
            if(!is.null(parent)) {
              parent <- getBlock(parent)
              if(!is(parent,"GtkWindow"))
                parent <- parent$GetWindow()
              if(!is(parent,"GtkWindow"))
                parent <- NULL          # give up
            }
              
            
            ## use message dialog for Gtk
            dlg = gtkMessageDialogNew(
              parent = parent,
              flags = 0,
              buttons = button,
              type=icon,
              message[1]
              )

            ## secret bit. Needs API! If message has length more than
            ## 1, use rest for secondary text.
            if(length(message) > 1)
              dlg['secondary-text'] <- paste(message[-1], collapse = "\n")
            
            dlg$SetTitle(title)
            dlg$GrabFocus()
            dlg$GetWindow()$Raise()
            dlg$setDefaultResponse(GtkResponseType["ok"])
            
            ## run in modal mode
            response = dlg$Run()
            h = list(obj=dlg, ref=dlg, action=action)
            if(response == GtkResponseType["cancel"] ||
               response == GtkResponseType["close"] ||
               response == GtkResponseType["delete-event"]) {
              dlg$Destroy()
              invisible(FALSE)
            } else if(response == GtkResponseType["ok"]) {
              if(!is.null(handler)) handler(h)
              dlg$Destroy()
              invisible(TRUE)
            } else {
              gwCat("Don't know this response")
              print(response)
              dlg$Destroy()
              invisible(NA)
            }
          })

## if OK then run handler, else not
setMethod(".gconfirm",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   message,
                   title = "Confirm",
                   icon = c("info", "warning", "error", "question"),
                   parent=NULL,
                   handler = NULL,
                   action = NULL,
                   ...
                   ) {
            
            if(missing(icon)) icon="question"
            icon = match.arg(icon)
            icon = Paste("GTK_MESSAGE_",toupper(match.arg(icon)))
            
            icon = "GTK_MESSAGE_QUESTION"
            buttons = "GTK_BUTTONS_OK_CANCEL"
            
            ## parent
            if(!is.null(parent)) {
              parent <- getBlock(parent)
              if(!is(parent,"GtkWindow"))
                parent <- parent$GetWindow()
              if(!is(parent,"GtkWindow"))
                parent <- NULL          # give up
            }
              


            dlg = gtkMessageDialogNew(
              parent = parent,
              flags = 0,
              buttons = buttons,
              type=icon,
              message[1]
              )
            ## secret bit. Needs API! If message has length more than
            ## 1, use rest for secondary text.
            if(length(message) > 1)
              dlg['secondary-text'] <- paste(message[-1], collapse = "\n")
            
            dlg$SetTitle(title)            
            dlg$GrabFocus()
            dlg$GetWindow()$Raise()
            dlg$setDefaultResponse(GtkResponseType["ok"]) # fails -- need to use gtkDialog directly
            
            ## add callback to close
            close.handler = function(h,...) h$obj$Destroy()
            
            ## run in modal mode
            response = dlg$Run()
            h = list(obj=dlg, action=action)
            if (response == GtkResponseType["close"] ||
                response == GtkResponseType["delete-event"] ||
                response == GtkResponseType["cancel"]) {
              dlg$Destroy()
              invisible(FALSE)
            } else if(response == GtkResponseType["ok"]) {
              if(!is.null(handler)) handler(h)
              dlg$Destroy()
              invisible(TRUE)
            } else {
              gwCat("Don't know this response")
              print(response)
              dlg$Destroy()
              invisible(NA)
            }
            
          })


## Add input to the above
## h,... in handler has componets action, input (for value)
setMethod(".ginput",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   message,
                   text="",
                   title = "Input",
                   icon = c("info", "warning", "error", "question"),
                   parent=NULL,                   
                   handler = NULL,
                   action = NULL,
                   ...
                   ) {
            
            icon = Paste("GTK_MESSAGE_",toupper(match.arg(icon)))

            ## parent
            if(!is.null(parent)) {
              parent <- getBlock(parent)
              if(!is(parent,"GtkWindow"))
                parent <- parent$GetWindow()
              if(!is(parent,"GtkWindow"))
                parent <- NULL          # give up
            }
            

            ## use message dialog for Gtk
            dlg = gtkMessageDialogNew(
              parent = parent,
              flags = 0,
              buttons = "GTK_BUTTONS_OK_CANCEL",
              type=icon,
              message[1]
              )
            dlg$SetTitle(title)
            dlg$setDefaultResponse(GtkResponseType["ok"])
            
            ## secret bit. Needs API! If message has length more than
            ## 1, use rest for secondary text.
            if(length(message) > 1)
              dlg['secondary-text'] <- paste(message[-1], collapse = "\n")
            dlg$GrabFocus()
            dlg$GetWindow()$Raise()

            
            group = ggroup(horizontal=FALSE)
#            glabel(message, container=group)
            input = gedit(text,container=group)
            
            
            ## find the area to pack the entry widget
##            dlg$GetVbox()[[1]]$PackStart(getBlock(group))
            dlg$GetVbox()$PackStart(getBlock(group)) 
            ##  dlg$GetVbox()[[2]]$GetWidget()$PackStart(group$ref) 
            ##  dlg$GetVbox()$PackStart(group$ref)

            ## set as default
            widget <- getWidget(input)
            widget['can-default'] <- TRUE
            widget$grabFocus()
            widget$grabDefault()

            
            ## run in modal mode
            response = dlg$Run()
            h = list(obj=dlg, ref=dlg, action=action, input=svalue(input))
            if(response == GtkResponseType["cancel"] ||
               response == GtkResponseType["close"] ||
               response == GtkResponseType["delete-event"]) {
              dlg$Destroy()
              invisible(NA)
            } else if(response == GtkResponseType["ok"]) {
              if(!is.null(handler)) handler(h)
              val = svalue(input)
              dlg$Destroy()
              ## input is widget, return value of widget
              invisible(val)
            } else {
              gwCat("Don't know this response")
              print(response)
              dlg$Destroy()
              invisible(NA)
            }
            
          })

## add a widget to the dialog. This is modal
setMethod(".gbasicdialog",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   title = "Dialog",
                   widget,
                   parent=NULL,                   
                   handler = NULL,
                   action = NULL,
                   ...
                   ) {
            
            
            ## parent
            if(!is.null(parent)) {
              parent <- getBlock(parent)
              if(!is(parent,"GtkWindow"))
                parent <- parent$GetWindow()
              if(!is(parent,"GtkWindow"))
                parent <- NULL          # give up
            } else {
              parent <- gtkWindowNew(show=FALSE)
            }
            
            ## button order
            theArgs <- list(...)
            buttons <- getWithDefault(theArgs$buttons, c("ok","cancel"))

            buttonMap <- function(name) {
              if(name == "ok")
                list("gtk-ok", GtkResponseType["ok"])
              else if(name =="yes")
                list("gtk-yes", GtkResponseType["ok"])
              else if(name == "cancel")
                list("gtk-cancel", GtkResponseType["cancel"])
              else if(name == "close")
                list("gtk-close", GtkResponseType["close"])
              else if(name =="no")
                list("gtk-no", GtkResponseType["cancel"])
              else
                list("gtk-yes", GtkResponseType["ok"])
            }

            l <- list(title=title, parent=parent, flags=c("modal"))
            for(i in buttons) {
              m <- buttonMap(i)
              l[[length(l) + 1]] <- m[[1]]
              l[[length(l) + 1]] <- m[[2]]
            }

            dlg <- do.call("gtkDialog", l)
              
            ## dlg = gtkDialog(title,
            ##   parent=parent,
            ##   c("modal"),
            ##   "gtk-ok", GtkResponseType["ok"],
            ##   "gtk-cancel", GtkResponseType["cancel"]
            ##   )
            ## dlg$SetTitle(title)
            dlg$GrabFocus()
            dlg$GetWindow()$Raise()
            
            
            tag(widget,"dlg") <- dlg
            
            ## group to pack widget in
            group = ggroup()
            add(group, widget, expand=TRUE)
            
            ## find the area to pack the entry widget
            dlg$GetVbox()$PackStart(getBlock(group))
            
            ## run in modal mode
            response = dlg$Run()
            h = list(obj=widget, action=action)
            if(response == GtkResponseType["cancel"] ||
               response == GtkResponseType["close"] ||
               response == GtkResponseType["delete-event"]) {
              ## cancel action
              dlg$Destroy()
              return(FALSE)
            } else if(response == GtkResponseType["ok"]) {
              if(!is.null(handler))
                handler(h)
              dlg$Destroy()
              return(TRUE)              # was widget, but TRUE now
            } else {
              ## default action
              gwCat("Don't know this response")
              print(response)
              dlg$Destroy()
              invisible(NA)
            }
            
          })

## with no paret
setClass("gBasicDialogNoParentRGtk",
         contains="gContainerRGtk",
         prototype=prototype(new("gContainerRGtk"))
         )

setMethod(".gbasicdialognoparent",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   title = "Dialog",
                   parent=NULL,                   
                   handler = NULL,
                   action = NULL,
                   ## do.buttons=TRUE
                   ...
                   ) {

            ## parent
            if(!is.null(parent)) {
              parent <- getBlock(parent)
              if(!is(parent,"GtkWindow"))
                parent <- parent$GetWindow()
              if(!is(parent,"GtkWindow"))
                parent <- NULL          # give up
            } else {
              parent <- gtkWindowNew(show=FALSE)
            }
            



            theArgs <- list(...)
            buttons <- getWithDefault(theArgs$buttons, c("ok","cancel"))
            
            buttonMap <- function(name) {
              if(name == "ok")
                list("gtk-ok", GtkResponseType["ok"])
              else if(name =="yes")
                list("gtk-yes", GtkResponseType["ok"])
              else if(name == "cancel")
                list("gtk-cancel", GtkResponseType["cancel"])
              else if(name == "close")
                list("gtk-close", GtkResponseType["close"])
              else if(name =="no")
                list("gtk-no", GtkResponseType["cancel"])
              else
                list("gtk-yes", GtkResponseType["ok"])
            }

            l <- list(title=title, parent=parent, flags=c("modal"), show=FALSE)
            for(i in buttons) {
              m <- buttonMap(i)
              l[[length(l) + 1]] <- m[[1]]
              l[[length(l) + 1]] <- m[[2]]
            }

            ## do buttons?
            do.buttons <- getWithDefault(theArgs$do.buttons, TRUE)
            if(do.buttons) {
              dlg <- do.call("gtkDialog", l)
            } else {
              dlg <- gtkDialogNew(show=FALSE)
              dlg$setTransientFor(parent)
              ## hide separator and button box. Uses internals -- bad idea if widget changes
              sapply(dlg$getChildren()[[1]]$getChildren(), gtkWidgetHide)
            }
            ## dlg = gtkDialog(title,
            ##   parent=parent,
            ##   flags = 0,
            ##   "gtk-ok", GtkResponseType["ok"],
            ##   "gtk-cancel", GtkResponseType["cancel"],
            ##   show=FALSE)
            dlg$SetTitle(title)
            dlg$setDefaultResponse(GtkResponseType["ok"])
            
            obj <- new("gBasicDialogNoParentRGtk",
                       block=dlg, widget=dlg, toolkit=guiToolkit("RGtk2"))
            tag(obj,"handler") <- handler
            tag(obj,"action") <- action

            return(obj)
          })

setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",
                    obj="gBasicDialogNoParentRGtk", value="guiWidget"),
          function(obj, toolkit, value, ...) {
            .add(obj, toolkit, value@widget, ...)
          })
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",
                    obj="gBasicDialogNoParentRGtk", value="gWidgetRGtk"),
           function(obj, toolkit, value, ...) {
             tag(obj,"widget") <- value
             
             ## group to pack widget in
             group <- gtkHBox(spacing=5)
             group$PackStart(getBlock(value))
            
             ## find the area to pack the entry widget
             dlg <- getWidget(obj)
             dlg$GetVbox()$PackStart(group)
            
             dlg$GrabFocus()
#             dlg$GetWindow()$Raise()

            
          })

##' dispose of dialog
setMethod(".dispose", signature(toolkit="guiWidgetsToolkitRGtk2",
                                obj="gBasicDialogNoParentRGtk"),
          function(obj, toolkit,...) {
            dlg <- getWidget(obj)
            dlg$Destroy()
            return(TRUE) 
          })

setMethod(".visible",
                 signature(toolkit="guiWidgetsToolkitRGtk2",
                           obj="gBasicDialogNoParentRGtk"),
                 function(obj, toolkit, set=NULL, ...) {

                   if(as.logical(set)) {

                     dlg <- getWidget(obj)
                     handler <- tag(obj,"handler")
                     action <- tag(obj,"action")
                     widget <- tag(obj,"widget")
                     
                     ## run in modal mode
                     dlg$show()
                     response = dlg$Run()


                     h = list(obj=widget, action=action)

                     if(response == GtkResponseType["cancel"] ||
                        response == GtkResponseType["close"] ||
                        response == GtkResponseType["delete-event"]) {
                       ## cancel action
                       dlg$Destroy()
                       return(FALSE)
                     } else if(response == GtkResponseType["ok"]) {
                       if(!is.null(handler))
                         handler(h)
                       dlg$Destroy()
                       return(TRUE)              # was widget, but TRUE now
                     } else if(response == GtkResponseType["delete-event"]) {
                       ## window manager close
                       return(invisible(FALSE))
                     } else if(response == GtkResponseType["none"]) {
                       ## dispose() call
                       return(invisible(FALSE))
                     } else {
                       ## default action
                       message("Don't know this response: ", response)
                       dlg$Destroy()
                       return(invisible(NA))
                     }
                   } else {
                     gwCat("gbasicdialog: call visible(obj,set=TRUE) to see.\n")
                     return(invisible(NA))
                   }
                 })
                 


setMethod(".galert",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   message,
                   title = "message",
                   delay = 3,
                   parent=NULL,
                   ...
                   ) {
            force(toolkit)

            ## insert an info bar here?
            if(is(parent, "gWindow")) {
              ib <- gtkInfoBar(show=FALSE)
              ib$setNoShowAll(TRUE)
              ib$setMessageType("warning")
              ca <- ib$getContentArea()
              ca$packStart(gtkLabel(message))
              b <- gtkButton(stock.id="cancel")
              ca$packEnd(b, FALSE)

              ## add to ibg
              ibg <- tag(parent, "infobargroup")
              ourgroup <- ggroup(container=ibg, expand=TRUE, horizontal=FALSE)
              add(ourgroup, ib, expand=TRUE)
              visible(ibg) <- TRUE
              ib$show()

              closeBar <- function(...) {
                ib$hide()
                visible(ibg) <- FALSE
                delete(ibg, ourgroup)
                FALSE
              }
              gSignalConnect(b, "clicked", closeBar)
              gTimeoutAdd(delay*1000, closeBar)
            } else {
              w <- gwindow(title, width=250, height=100, parent = parent)
              g <- ggroup(container =w)
              l <- gbutton("  ", container =g)
              getToolkitWidget(l)$modifyBg(GtkStateType['normal'], color="red")
              label <- glabel(message, container =g, expand=TRUE)
              font(label) <- c("weight"="bold")
              gimage(filename="close",dirname="stock", container =g, handler = function(h,...) dispose(w))
              
              addHandlerIdle(w, handler = function(h,...) dispose(w),
                             interval = as.numeric(delay)*1000)
            }

            invisible(w)
          })
 
