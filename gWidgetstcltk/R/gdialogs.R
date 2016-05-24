## some dialogs for R
## dialogs don't get windows, they make them
## dialogs are modal
## dialogs return their value -- not an object. so source(gfile()) should work

## we don't implement gbasiddialog. -- how to do so not clear?

## TODO:

## used to create all three dialogs
tcltkDialog = function(
  message,
  text = "",
  title = "Input",
  icon = c("info","warning","error","question"),
  type = c("message","confirm","input"),
  parent = NULL,
  handler = NULL,
  action = NULL,
  ...
  ) {

  ## top level widnow
  dlg <- tktoplevel()
  f <- ttkframe(dlg, padding=3)
  tkpack(f, expand=TRUE, fill="both")
  
  if(!is.null(parent)) {
    parent <- getBlock(parent) ## needs to be top level window
    parent <- getTopParent(parent)
    curgeo <- tclvalue(tkwm.geometry(parent))
    ## widthXheight+xpos+ypos
    pos <- unlist(strsplit(curgeo, "\\+"))
    sz <- unlist(strsplit(pos[1],"x"))
    xpos = as.numeric(pos[2]); ypos=as.numeric(pos[3])
    tkwm.geometry(dlg,paste("+",xpos+10,"+",ypos+10,sep="")) # shift
    
    tkwm.transient(dlg, parent) # set transient
    tkbind(parent,"<Destroy>",function(...) tkdestroy(dlg))
  }

      
 
  
  ## set up dlg window
  tkwm.deiconify(dlg)
#  tkgrab.set(dlg) ## was giving errors
  tkfocus(dlg)
  tkwm.title(dlg,title)
  tkwm.resizable(dlg, FALSE, FALSE)

  dlgframe <- ttkframe(f, padding=3)
  tkpack(dlgframe, expand=TRUE, fill="both")
                       

   
  ## set up icon
  ## These are stupid icons!!!
  icon = match.arg(icon)
  allIcons = getStockIcons()
  iconFile = switch(icon,
    "warning"=allIcons$alert,
    "error" = allIcons$error,
    "question" = allIcons$help,
    allIcons$ok
    )
  imageID = paste("gdialogs",as.character(runif(1)),sep="")
  tcl("image","create","photo",imageID,file=iconFile)  
  icon = ttklabel(dlgframe,image=imageID)
  tkgrid(icon,row=0,column=0)

  ## set up label
  if(missing(message) || is.null(message))
    message <- ""
  l <- ttklabel(dlgframe, text = paste(as.character(message), sep="\n"))
  tkgrid(l, row=0, column = 1, stick ="nw", padx=25, pady=5)


  
  ## entry widget for input
  if(type == "input") {
    textEntryVarTcl <- tclVar(text)
    textEntryWidget <-
      ttkentry(dlgframe,
              width=max(25,as.integer(1.3*nchar(text))),
              textvariable=textEntryVarTcl)
    tkgrid(textEntryWidget,row = 1, column=1,stick="nw", padx=5,pady=5)
  }
  
  ## what to return? TRUE or FALSE or string for ginput
  ReturnVal <- FALSE
  
  
  onOK <- function() {
    if(type == "input") 
      ReturnVal <<- tclvalue(textEntryVarTcl)
    else
      ReturnVal <<- TRUE
    
    ## call handler if asked
    if(!is.null(handler)) 
      handler(list(obj=NULL, action=action, input=ReturnVal))
    
    tkgrab.release(dlg)
    tkdestroy(dlg)
  }
  onCancel <- function(){
    if(type == "input")
      ReturnVal <<- NA
    else
      ReturnVal <<- FALSE
    tkgrab.release(dlg)
    tkdestroy(dlg)
  }
  
  gp <- ttkframe(f)

  OK.but     <-ttkbutton(gp,text="   OK   ",command=onOK, state="active")
  Cancel.but <-ttkbutton(gp,text=" Cancel ",command=onCancel)

  tkpack(gp, fill="y")
  if(type == "confirm" || type == "input")
    tkpack(Cancel.but,side="left")
  tkpack(OK.but,side="left")
  
  if(type == "input")
    tkfocus(textEntryWidget)            # set focus
  else
    tkfocus(OK.but)

  tkbind(dlg, "<Destroy>", function() {
    tkgrab.release(dlg)
  })
  if(type == "input")
    tkbind(textEntryWidget, "<Return>", onOK)

  tkwait.window(dlg)

  invisible(ReturnVal)
}




setMethod(".gmessage",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   message,
                   title = "message",
                   icon = c("info","warning","error","question"),
                   parent = NULL,
                   handler = NULL,
                   action = NULL,
                   ...
                   ) {

            icon = match.arg(icon)
            l <- list(icon=icon,
                      message=gettext(message[1]),
                      title = title,
                      type="ok")
            if(length(message) > 1)
              l$detail=gettext(message[2])
            
            if(!is.null(parent))
              l$parent <- getWidget(parent)

            out <- do.call("tkmessageBox",l)

            
            if(is.logical(out) && out  && !is.null(handler)) {
              h = list()
              h$obj=NULL; h$action=action
              handler(h)
            }
            
            return(out)
            
            ## ## old
            ## return(tcltkDialog(
            ##                    message,
            ##                    title=title,
            ##                    icon=icon,
            ##                    type="message",
            ##                    parent = parent,
            ##                    handler=handler,
            ##                    action=action,
            ##                    ...))

##             icon = match.arg(icon)
            
##             ret = tkmessageBox(
##               message=message,
##               title=title,
##               icon=icon)
##             if(as.character(ret) == "ok")
##               TRUE
##             else
##               FALSE
          })
  
## if OK then run handler, else not
setMethod(".gconfirm",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   message,
                   title = "Confirm",
                   icon = c("info", "warning", "error", "question"),
                   parent = NULL,
                   handler = NULL,
                   action = NULL,
                   ...
                   ) {
            icon = match.arg(icon)
            
            l <- list(icon=icon,
                      message=gettext(message[1]),
                      title = title,
                      type="yesno")

            if(length(message) > 1)
              l$detail=gettext(message[2])
            
            if(!is.null(parent))
              l$parent <- getWidget(parent)
            out <- do.call("tkmessageBox",l)
            val <- switch(as.character(out),
                          "yes"=TRUE,
                          "no" = FALSE,
                          FALSE)

            if(val && !is.null(handler)) {
              h = list()
              h$obj=NULL; h$action=action
              handler(h)
            }
            
            return(val)

            ## return(tcltkDialog(
            ##                    message,
            ##                    title=title,
            ##                    icon=icon,
            ##                    type="confirm",
            ##                    parent = parent,
            ##                    handler=handler,
            ##                    action=action,
            ##                    ...))

##             icon = match.arg(icon)

##             ret = tkmessageBox(
##               message=message, 
##               title=title,
##               icon=icon,
##               type="yesnocancel"
##               )

##             val = switch(as.character(ret),
##               "yes"=1,
##               "no"=0,
##               "cancel"=-1)

##             if(!is.null(handler)) {
##               h = list()
##               h$obj=NULL; h$action=action
##               handler(h)
##             }
              
            
##             return(val)

          })

 
## Add input to the above
## h,... in handler has componets action, input (for value)
setMethod(".ginput",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   message,
                   text = "",
                   title = "Input",
                   icon = c("info","warning","error","question"),
                   parent = NULL,
                   handler = NULL,
                   action = NULL,
                   ...
                   ) {

            return(tcltkDialog(
                               message,
                               text = text,
                               title=title,
                               icon=icon,
                               type="input",
                               parent = parent,
                               handler=handler,
                               action=action,
                               ...))

          })

## add a widget to the dialog. This is modal
## see next one for one that gets called here
setMethod(".gbasicdialog",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   title = "Dialog",
                   widget,
                   parent = NULL,
                   handler = NULL,
                   action = NULL,
                   ...
                   ) {
            message(gettext("gbasiddialog isn't implemented in tcltk"),"\n")
            return()
          })

## with no paret
setClass("gBasicDialogNoParenttcltk",
         contains="gWindowtcltk",
         prototype=prototype(new("gContainertcltk"))
         )

setMethod(".gbasicdialognoparent",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   title = "Dialog",
                   parent=NULL,                   
                   handler = NULL,
                   action = NULL,
                   ...
                   ) {
            
            dlg <- gwindow(title, parent=parent, visible=FALSE)
            tt <- dlg@widget@widget
            
            g <- ggroup(container = dlg, horizontal=FALSE, expand=TRUE)
            
            obj <- new("gBasicDialogNoParenttcltk",
                       block=dlg, widget=g, toolkit=guiToolkit("tcltk"))
            tag(obj,"handler") <- handler
            tag(obj,"action") <- action
            tag(obj,"tt") <- tt

            args <- list(...)
            tag(obj, "do.buttons") <- getWithDefault(args$do.buttons, TRUE)
            
            tkbind(tt, "<Destroy>", function() {
              tkgrab.release(tt)
            })

            return(obj) 
          })


setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",
                    obj="gBasicDialogNoParenttcltk", value="guiWidget"),
          function(obj, toolkit, value, ...) {
            .add(obj, toolkit, value@widget, ...)
          })

setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",
                    obj="gBasicDialogNoParenttcltk", value="gWidgettcltk"),
           function(obj, toolkit, value, ...) {
             add(obj@widget, value, ...)
             ## keep these around
             tag(obj,"widget") <- value
          })

##' close window
setMethod(".dispose",
                 signature(toolkit="guiWidgetsToolkittcltk",
                           obj="gBasicDialogNoParenttcltk"),
                 function(obj, toolkit, ...) {
                   flag <- tag(obj, "flag")
                   tclvalue(flag) <- "destroy"
                 })

setMethod(".visible",
                 signature(toolkit="guiWidgetsToolkittcltk",
                           obj="gBasicDialogNoParenttcltk"),
                 function(obj, toolkit, set=NULL, ...) {

                   if(as.logical(set)) {

                     handler <- tag(obj,"handler")
                     action <- tag(obj,"action")
                     widget <- tag(obj,"widget")
                     tt <- tag(obj,"tt")

                     dlg <- obj@block
                     g <- obj@widget



                     
                     ## we use tclwait.variable, rather than window
                     ## with window, we need to destroy widget before returning loop
                     ## and then widget is destroyed before we can use it.
                     flag <- tclVar("")
                     tag(obj, "flag") <- flag
                     
                     ## bind to destroy event
                     tkwm.protocol(dlg@widget@block, "WM_DELETE_WINDOW", function() {
                       tclvalue(flag) <- "destroy"
                     })

                     
                     ans <- FALSE

                     if(tag(obj, "do.buttons")) {
                       buttonGroup = ggroup(container=g, expand=TRUE, fill="x") ## just x XXX
                       addSpring(buttonGroup)
                       OKbutton = gbutton("OK",container=buttonGroup,action = tt,
                         handler=function(h,...) {
                           ans <<- TRUE
                           tkgrab.release(h$action)
                           tclvalue(flag) <- "destroy"
                         })
                       addSpace(buttonGroup, 10)
                       Cancelbutton = gbutton("Cancel",container=buttonGroup, action=tt,
                         handler=function(h,...) {
                           ans <<- FALSE
                           tkgrab.release(h$action)
                           tclvalue(flag) <- "destroy"
                         })
                     }

                     ## make window visible and on top of stack
                     visible(dlg) <- TRUE
                     focus(dlg) <- TRUE
                     ## make modal
                     tkwait.variable(flag)
                       
                     ## process response
                     if(ans) {
                       ## yes
                       if(!is.null(handler)) {
                         handler(list(obj=widget,action=action))
                       }
                       dispose(dlg)
                       return(ans)
                     } else {
                       ## no
                       dispose(dlg)
                       return(ans)
                     }
                   } else {
                     ## nothing
                     dispose(dlg)                     
                     return(NA)
                   }
                 })




## alert
setMethod(".galert",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   message,
                   title = "message",
                   delay = 3,
                   parent=NULL,
                   ...
                   ) {
            force(toolkit)
            
            w <- gwindow(title, width=250, height=50, parent = parent)
            g <- ggroup(container = w)
            l <- glabel("  ", container = g)
            label <- glabel(message, container = g, expand=TRUE)
            font(label) <- c("weight"="bold")
            gimage(filename="dismiss",dirname="stock", container = g, handler = function(h,...) dispose(w))
            
            addHandlerMouseMotion(label, handler = function(h,...) dispose(w))
            
            w
          })
