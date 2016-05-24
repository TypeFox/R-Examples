w.textentry <- function(parent, label, text, action=I, pos=NULL, foreground=NULL, background=NULL, font=NULL, 
  width=20, keydown=FALSE)
{
  widget <- w.createwidget(parent, pos, background)
  widget$.type = "textentrys"
  widget$.var <- c()
  widget$.label <- list()   
  widget$.text <- list()
  
  f <- function(...) 
  { 
    variable <- c()
    for (j in (1:length(label))) 
    { 
      variable[j] <- as.character(handshake(tclvalue, widget$.var[[j]]))
    }
    if (length(label) == 1) { action(as.character(variable)) } else { action(variable) }
  }

  widget$.f <- f

  for (i in 1:length(label))
  { 
    widget$.var[i] <- list(handshake(tclVar, text[[i]]))
    labeller <- w.text(widget, label[[i]], NA, 
      pos=list(
        column=0, 
        row=i-1, 
        sticky="w", 
        cweight=1,
        width=as.integer(handshake(.Tcl, paste('font measure systemfont "', label[[i]], '"', sep="") )),
        height=as.integer(handshake(.Tcl, 'font metrics systemfont -linespace'))
      ), foreground, background, font)
    entry <- w.createwidget(widget, 
      pos=list(
        column=1, 
        row=i-1, 
        sticky="ew", 
        cweight=100,
        width=width*as.integer(handshake(.Tcl, paste('font measure systemfont "1"', sep="") )),
        height=as.integer(handshake(.Tcl, 'font metrics systemfont -linespace'))
      ), background)
    entry$.type = "textentry"
    entry$.widget <- handshake(tkentry, widget$.handle, width=width, textvariable=widget$.var[[i]])
    if (keydown) { handshake(tkbind, entry$.widget, "<KeyRelease>", f) } # run the function for every keypress - inadvisable
    handshake(tkbind, entry$.widget, "<Key-Return>", f) # run the function when enter is pressed
    w.appearancewidget(entry, font, foreground, background) # this applies to the textentry only
    widget$.label[i] <- list(labeller)
    widget$.entry[i] <- list(entry)
  }
  invisible(widget)
}
  
rp.textentry <- function(panel, variable, action = I, labels = NULL, names = labels,
                         title = NULL, initval = rep(NA, length(labels)),
                         pos = NULL, foreground = NULL, background = NULL, 
                         font = NULL, width = 20, keydown = FALSE, 
                         parentname = deparse(substitute(panel)), 
                         name = paste("textentry", .nc(), sep=""), ...) {
                         	
  if (!exists(panel$panelname, .rpenv, inherits = FALSE)) { # if the panelname is not set then
     panelname = deparse(substitute(panel)) # the panel name should be the panel deparse subst'ed
# 13/03/2012 these lines are not commented out in previous version
#    panel <- rp.control.get(panelname, panel) # now get the panel
#    panel$panelname = panelname # now set the panelname properly
#    assign(panelname, panel, envir=.rpenv) # now send back the panel
  } 
  else { 
    panelname = panel$panelname 
# 13/03/2012 these lines are not commented out in previous version
#    panel <- rp.control.get(panelname, panel) # now get the panel
  }
  
  varname <- deparse(substitute(variable))
  if (!rp.isnull(panelname, varname))
     variable <- rp.var.get(panelname, varname)
  else { 
    variable        <- initval
    names(variable) <- labels
    rp.var.put(panelname, varname, variable, labels)
  } 
  
  if ((is.null(title))  && (is.null(labels))) labels <- varname
  if ((!is.null(title)) && (is.null(labels))) labels <- title
  if (is.null(pos) & length(list(...)) > 0)   pos    <- list(...)
  
  f <- function(val) {
    valexisting <- rp.var.get(panelname, varname)
    names(val) <- names(valexisting)
    rp.var.put(panelname, varname, val)
    panel <- rp.control.get(panelname)
    panel <- action(panel)
    rp.control.put(panelname, panel)
  }
   
  if (rp.widget.exists(panelname, parentname))
     parent <- rp.widget.get(panelname, parentname)
  else 
     parent <- panel
  if (is.list(pos) && !is.null(pos$grid))
     parent <- rp.widget.get(panelname, pos$grid)
  
  widget <- w.textentry(parent, labels, text = variable, action = f, pos, 
                        foreground, background, font, width, keydown)
  rp.widget.put(panelname, name, widget)

  if (.rpenv$savepanel) rp.control.put(panelname, panel) # put the panel back into the environment
  invisible(panelname)
}
