w.text <- function(parent, text, action = NULL,
                   pos = NULL, foreground = NULL, background = NULL,
                   font = NULL, width = NULL) {

  widget <- w.createwidget(parent, pos, background)
  widget$.type = "text"
  
  widget$.var <- handshake(tclVar, text)
  f <- function(...) { action() }
  if (is.null(width))
  	 widget$.widget <- handshake(tklabel, parent$.handle,
  	                      text=handshake(tclvalue, widget$.var))

  else
     widget$.widget <- handshake(tklabel, parent$.handle,
                          text=handshake(tclvalue, widget$.var), width = width)
  handshake(tkconfigure, widget$.widget, textvariable = widget$.var)
  if (!is.null(action)) handshake(tkbind, widget$.widget, "<Button-1>", f)
  w.appearancewidget(widget, font, foreground, background)
  invisible(widget)
}

w.text.change <- function(widget, text) {
# not really possible to handshake this
  tclvalue(widget$.var) <- text
}

rp.text <- function(panel, text, pos = NULL, action = I,
                    foreground = NULL, background = NULL,
                    font = NULL, width = NULL, parentname=deparse(substitute(panel)), 
                    name = paste("text", .nc(), sep = ""), ...) {
  if (!exists(panel$panelname, .rpenv, inherits = FALSE))
    panelname = deparse(substitute(panel))
# 13/03/2012 these lines are not commented out in previous version
#    panel <- rp.control.get(panelname, panel) # now get the panel
#    panel$panelname = panelname # now set the panelname properly
#    assign(panelname, panel, envir=.rpenv) # now send back the panel
  else 
    panelname = panel$panelname 
# 13/03/2012 these lines are not commented out in previous version
#    panel <- rp.control.get(panelname, panel) # now get the panel

  if (is.null(pos) & length(list(...)) > 0) pos <- list(...)

  f <- function() {
# 13/03/2012 this next line was not commented out in the previous version
#    panel <- rp.control.get(panelname)
    panel <- action(panel)
    rp.control.put(panelname, panel)    
  }  

  if (rp.widget.exists(panelname, parentname))
  	 parent <- rp.widget.get(panelname, parentname)
  else
     parent <- panel
  if (is.list(pos) && !is.null(pos$grid)) parent <- rp.widget.get(panelname, pos$grid)
  
  widget <- w.text(parent, text, action=f, pos, foreground, background, font, width)
  rp.widget.put(panelname, name, widget)

  if (.rpenv$savepanel) rp.control.put(panelname, panel)
  invisible(panelname)
}

rp.text.change <- function(panel, name, text) {
  if (!exists(panel$panelname, .rpenv, inherits = FALSE))
    panelname <- deparse(substitute(panel))
# 13/03/2012 these lines are not commented out in previous version
#    panel <- rp.control.get(panelname, panel) # now get the panel
#    panel$panelname = panelname # now set the panelname properly
#    assign(panelname, panel, envir=.rpenv) # now send back the panel
  else 
    panelname <- panel$panelname 
# 13/03/2012 these lines are not commented out in previous version
#    panel <- rp.control.get(panelname, panel) # now get the panel
  w.text.change(rp.widget.get(panelname, name), text)
  invisible()
}
