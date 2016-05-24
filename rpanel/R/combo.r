w.combo <- function(parent, prompt=NULL, values, pos=NULL, action=I, foreground=NULL, background=NULL, font=NULL, editable=FALSE)
{
  widget <- w.createwidget(parent, pos, background)
  widget$.type = "combobox"
  f <- function(...) 
  { 
    val <- values[as.numeric(handshake(tcl, widget$.combo[[1]]$.widget,"getvalue"))+1]
    action(val)
  }
               
  handshake(tclRequire, "BWidget")
  
  maxlen <- function(values) {
    len <- 0
    for(i in 1:length(values)) {
#      if (nchar(values[i]) > len) {
#        len <- nchar(values[i])
#      }
      if (as.integer(handshake(.Tcl, paste('font measure systemfont "', values[i], '"', sep="") )) > len) {
        len <- as.integer(handshake(.Tcl, paste('font measure systemfont "', values[i], '"', sep="") ))
#        print(len)
      }
    }
    invisible(len)  
  }

  maxlen2 <- function(values) {
    len <- 0
    for(i in 1:length(values)) {
      if (nchar(values[i]) > len) {
        len <- nchar(values[i])
      }
#      if (as.integer(handshake(.Tcl, paste('font measure systemfont "', values[i], '"', sep="") )) > len) {
#        len <- as.integer(handshake(.Tcl, paste('font measure systemfont "', values[i], '"', sep="") ))
#        print(len)
#      }
    }
    invisible(len)  
  }

  label <- w.text(widget, prompt, NA, 
    pos=list(column=0, row=0, sticky="w", cweight=1, 
      width=as.integer(handshake(.Tcl, paste('font measure systemfont "', prompt, '"', sep="") )),
      height=as.integer(handshake(.Tcl, 'font metrics systemfont -linespace'))
    ), foreground, background, font)

  combo <- w.createwidget(widget, 
    pos=list(column=1, row=0, sticky="ew", cweight=100, 
      width=maxlen(values), 
      height=as.integer(handshake(.Tcl, 'font metrics systemfont -linespace'))
    ), background)
    
  combo$.type <- "combo"
  combo$.widget <- handshake(tkwidget, parent$.handle, "ComboBox", editable=FALSE, values=values, modifycmd=f, width=maxlen2(values),
     textvariable = tclVar(values[1]))
  w.appearancewidget(combo, font, foreground, background)
  widget$.label <- list(label)
  widget$.combo <- list(combo)

  invisible(widget)
}

rp.combo <- function(panel, variable, prompt=NULL, vals, initval=vals[1], pos=NULL, action=I, foreground=NULL, 
  background=NULL, font=NULL, editable=FALSE, parentname=deparse(substitute(panel)), name=paste("combo", .nc(), sep=""), ...) 
{
  if (!exists(panel$panelname, .rpenv, inherits = FALSE)) # if the panelname is not set then
  { 
    panelname = deparse(substitute(panel)) # the panel name should be the panel deparse subst'ed
# 13/03/2012 these lines are not commented out in previous version
#    panel <- rp.control.get(panelname, panel) # now get the panel
#    panel$panelname = panelname # now set the panelname properly
#    assign(panelname, panel, envir=.rpenv) # now send back the panel
  } 
  else 
  { 
    panelname = panel$panelname 
# 13/03/2012 these lines are not commented out in previous version
#    panel <- rp.control.get(panelname, panel) # now get the panel
  }
  
  varname = deparse(substitute(variable))
  if (!rp.isnull(panelname, varname)) { variable = rp.var.get(panelname, varname) } 
  else { variable = initval; rp.var.put(panelname, varname, variable); } 

  if (is.null(pos)) { if (length(list(...)) > 0) { pos <- list(...) } }

  f <- function(val)
  {
    rp.var.put(panelname, varname, val)
    panel <- rp.control.get(panelname)
    panel <- action(panel)
    rp.control.put(panelname, panel)
  }

  if (rp.widget.exists(panelname, parentname)) { parent <- rp.widget.get(panelname, parentname) }
  else { parent <- panel }
  if (is.list(pos)) { if (!is.null(pos$grid)) { parent <- rp.widget.get(panelname, pos$grid) } }
  
  widget <- w.combo(parent, prompt, vals, pos, action=f, foreground, background, font, editable)
  rp.widget.put(panelname, name, widget)
  if (.rpenv$savepanel) { rp.control.put(panelname, panel) } # put the panel back into the environment
  invisible(panelname)
}
