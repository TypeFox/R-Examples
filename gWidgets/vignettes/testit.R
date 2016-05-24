FROM = "gWidgetsRGtk <gWidgetsRGtk@gmail.com>"
buddyList = c("My Friend <myfriend@gmail.com>","My dog <mydog@gmail.com>")

Rmail = function(draft = NULL, ...) {
  ## Define main widgets, store in a list for ease of use
  widgets = list()
  widgets$to = gdroplist(c(), editable=TRUE)
  widgets$from = glabel(FROM, editable=TRUE)
  widgets$subject = gedit()
  widgets$text = gtext()
  
  ## Handle drafts. Either a list or a filename to source")
  ## The generic svalue() method makes setting values easy")
  if(!is.null(draft)) {
    if(is.character(draft)) {
      print(draft)
      sys.source(draft,envir=environment())
      print("smudget")
      print(draft)
    }
    if(is.list(draft)) {
      print("draft is list")
      sapply(c("to","from","subject","text"), function(i)  {
        cat("assign",i,"value of",draft[[i]],"\n")
        svalue(widgets[[i]]) <- draft[[i]]
      })
    }
  }
  
  ## Helper functions
  sendIt = function(...) {
    tmp = tempfile()
    
    cat("To:", svalue(widgets$to),"\n",file = tmp, append=TRUE)
    cat("From:", svalue(widgets$from),"\n", file=tmp, append=TRUE)
    cat("Subject:", svalue(widgets$subject),"\n", file=tmp, append=TRUE)
    cat("Date:", format(Sys.time(),"%d %b %Y %T %Z"),"\n", file=tmp, append=TRUE)
    cat("X-sender:", "R", file=tmp, append=TRUE)
    cat("\n\n", file=tmp, append=TRUE)
    cat(svalue(widgets$text), file=tmp, append=TRUE)
    cat("\n", file=tmp, append=TRUE)
    
    ## Use UNIX sendmail to send message
    system(paste("sendmail -t  <", tmp))
    ## Add To: to buddyList
    if(exists("buddyList"))
    assign("buddyList", unique(c(buddyList,svalue(widgets$to))), inherits=TRUE)
    
    ## Close window, delete file
    unlink(tmp)
    dispose(window)
  }

  ## Function to save a draft to the file draft.R
  saveDraft = function(...) {
    draft = list()
    sapply(c("to","from","subject","text"), function(i) 
    draft[[i]] <<- svalue(widgets[[i]])
    )
    dump("draft","draft.R")
    cat("Draft dumped to draft.R\n")
  }
  
  ## A simple dialog
  aboutMail = function(...) gmessage("Sends a message")
  
  ## Make main window from top down
  
  window = gwindow("Compose mail")
  group = ggroup(horizontal=FALSE, spacing=0, container = window)
  ## Remove border
  svalue(group) <- 0
  
  ## Menubar is defined by a list
  menubarlist = list()
  menubarlist$File$Save$handler = saveDraft
  menubarlist$File$Send$handler = sendIt
  menubarlist$File$Quit$handler = function(...) dispose(window)
  menubarlist$File$Quit$icon = "quit"
  menubarlist$Help$About$handler = aboutMail
  add(group, gmenu(menubarlist))
  
  ## Toolbar is also defined by a list
  toolbarlist = list()
  toolbarlist$Send$handler = sendIt
  toolbarlist$Send$icon = "connect"
  toolbarlist$Save$handler = saveDraft
  toolbarlist$Save$icon = "save"
  add(group, gtoolbar(toolbarlist))
  
  
  ## Put headers in a glayout() container
  tbl = glayout(container = group)
  
  ## To: field. Looks for buddyList
  tbl[1,1] = glabel("To:")
  tbl[1,2] = widgets$to
  if(exists("buddyList")) widgets$to[] <- buddyList
  
  ## From: field. Click to edit value
  tbl[2,1] = glabel("From:")
  tbl[2,2] = widgets$from
  
  ## Subject: field. Handler updates window title
  tbl[3,1] = glabel("Subject:")
  tbl[3,2] = widgets$subject
  addHandlerKeystroke(widgets$subject, handler = function(h,...)
  svalue(window) = paste("Compose mail:",svalue(h$obj),collapse=""))
  
  ## Layout needs to be finalized
  visible(tbl) <- TRUE
  
  ## Add text box for message, but first some space
  addSpace(group, 5)
  add(group, widgets$text, expand=TRUE)
  
  ## That's it.
}
