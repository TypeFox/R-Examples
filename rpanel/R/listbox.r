w.listbox <- function(parent, title = NA, labels, rows = length(labels), initval = labels[1], 
                      action=I, pos=NULL, foreground=NULL, background="white", font=NULL,
                      sleep = 0.01) {
  if (is.na(title))
     widget <- w.createwidget(parent, pos, background)
  else
     widget <- w.createwidget(parent, pos, background, title)
  widget$.type <- "listbox"

  if (rows < length(labels)) { 
     widget$.widget <- handshake(tklistbox, parent$.handle, height=rows, selectmode = "single",
                                 yscrollcommand = function(...) handshake(tkset, scr, ...))
     scr <- handshake(tkscrollbar, parent$.handle, repeatinterval = 5,
                      command=function(...) handshake(tkyview, widget$.widget,...)) 
     w.appearancewidget(widget, font, foreground, background, scr)    
  }
  else {
     widget$.widget <- handshake(tklistbox, parent$.handle, height=rows, selectmode = "single")
     w.appearancewidget(widget, font, foreground, background)
  }

  selection <- 1
  
#  tkinsert(widget$.widget, "end", "test")

  for (i in (1:length(labels))) {
    Sys.sleep(sleep)
    handshake(tkinsert, widget$.widget, "end", as.character(labels[[i]]))
    if (labels[[i]] == initval) selection <- i - 1
  } 
 
  handshake(tkselection.set, widget$.widget, selection)
  
  f <- function(...) action(labels[as.numeric(handshake(tkcurselection, widget$.widget)) + 1])

  handshake(tkbind, widget$.widget, "<ButtonRelease-1>", f)
   
  invisible(widget)
}

rp.listbox <- function(panel, variable, vals, labels = vals, rows = length(labels),
         initval = vals[1], pos = NULL, title = deparse(substitute(variable)), action = I, 
         foreground = NULL, background = NULL, font = NULL,
         parentname = deparse(substitute(panel)), sleep = 0.01,
         name = paste("listbox", .nc(), sep = ""), ...) {
   	
  if (!exists(panel$panelname, .rpenv, inherits = FALSE)) { # if the panelname is not set then
     panelname <- deparse(substitute(panel)) # the panel name should be the panel deparse subst'ed
     # 13/03/2012 these lines are not commented out in previous version
     #    panel <- rp.control.get(panelname, panel) # now get the panel
     #    panel$panelname = panelname # now set the panelname properly
     #    assign(panelname, panel, envir=.rpenv) # now send back the panel
  } 
  else 
    panelname <- panel$panelname 
    # 13/03/2012 these lines are not commented out in previous version
    #    panel <- rp.control.get(panelname, panel) # now get the panel
  
  varname <- deparse(substitute(variable))
  if (!rp.isnull(panelname, varname))
     variable = rp.var.get(panelname, varname)
  else
     variable = initval; rp.var.put(panelname, varname, variable)

  if (is.null(pos) && length(list(...)) > 0) pos <- list(...)

  f <- function(val) {
     rp.var.put(panelname, varname, val)
     panel <- rp.control.get(panelname)
     panel <- action(panel)  
     rp.control.put(panelname, panel)
  }

  if (rp.widget.exists(panelname, parentname))
     parent <- rp.widget.get(panelname, parentname)
  else
     parent <- panel
  if (is.list(pos) && !is.null(pos$grid)) parent <- rp.widget.get(panelname, pos$grid)
  
  widget <- w.listbox(parent, title, labels, rows, initval = variable, action=f, pos,
                      foreground, background, font, sleep)
  rp.widget.put(panelname, name, widget)
  if (.rpenv$savepanel) rp.control.put(panelname, panel) # put the panel back into the environment
  invisible(panelname)
}
