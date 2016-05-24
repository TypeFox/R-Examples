w.radiogroup <- function(parent, title, vals, labels, initval = vals[1], action = I, pos = NULL,
                         foreground = NULL, background = "white", font = NULL)
{
  widget       <- w.createwidget(parent, pos, background, title)
  widget$.type <- "radiogroup"
  widget$.var  <- handshake(tclVar, initval)
  widget$.cb   <- list()
  # cb is used so that same disposal method can be used for both radiogroup and checkgroup

  for (i in (1:length(labels))) {
    f  <- function(...) action(handshake(tclvalue, widget$.var))
    rb <- w.createwidget(widget, 
      pos=list(
        column=0, 
        row=i-1, 
        sticky="w",
        width=as.integer(handshake(.Tcl, paste('font measure systemfont "', labels[[i]], '"', sep="") )),
        # 06/08 Offset of 2 is needed for unknown reasons - note 1 does not work
        height=2+as.integer(handshake(.Tcl, 'font metrics systemfont -linespace'))        
      ), background)
    rb$.type = "radiobutton"
    rb$.widget <- handshake(tkradiobutton, widget$.handle, command=f, text = labels[[i]],
                            variable = widget$.var)
    handshake(tkconfigure, rb$.widget, variable = widget$.var, value = vals[[i]])    
    w.appearancewidget(rb, font, foreground, background)    
    widget$.cb[i] <- list(rb)
  } 
  invisible(widget)
}

rp.radiogroup <- function(panel, variable, vals, labels = NULL, initval = vals[1],
      pos = NULL, title = deparse(substitute(variable)), 
      action = I, foreground = NULL, background = NULL, font = NULL, 
      parentname = deparse(substitute(panel)), name = paste("radiogroup", .nc(), sep = ""), ...) {

  # The following line appears to be redundant, but is not. There is a problem with 
  # deparse substitute and rtcltk where the wrong value is taken into tcltk. Printing 
  # or assigning the title gets round this problem.
  title <- paste(title)

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
  
# vals is not presently in use
  if (is.null(labels)) labels <- vals   

  varname <- deparse(substitute(variable))
  if (!rp.isnull(panelname, varname))
     variable <- rp.var.get(panelname, varname)
  else {
     variable <- initval
     rp.var.put(panelname, varname, variable)
# This is commented out but shows an interesting avenue - this returns the altered panel back to the calling environment.
# This was developed to deal with the rp.do problem - that it uses a local/global copy of panel and not the true version
# actually stored in .rpenv - this has been changed so that rp.do uses the .rpenv copy of the panel, not the global/local
# copy.
#    paa <- eval(parse(text=panelname), envir=.rpenv); 
#    assign(panelname, paa, envir=parent.frame()) 
  } 

  if (is.null(pos) & length(list(...)) > 0) pos <- list(...)
  if (is.null(labels)) labels <- varname
   
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
  if (is.list(pos)) { if (!is.null(pos$grid)) { parent <- rp.widget.get(panelname, pos$grid) } }
  
  widget <- w.radiogroup(parent, title, vals, labels, initval = variable, action=f,
                         pos, foreground, background, font)

  rp.widget.put(panelname, name, widget)
  rp.control.get(panelname)  
  
  if (.rpenv$savepanel) rp.control.put(panelname, panel) # put the panel back into the environment
  invisible(panelname)
}
