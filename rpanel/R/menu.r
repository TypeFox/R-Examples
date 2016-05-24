w.menu <- function(panel, labels, action=I, foreground=NULL, background=NULL, font=NULL) 
{
  f <- function(option) { function(...) { action(option) } } # function must be declared inside function
  
  menu <- list(.type="menu")
  menu$.widget <- handshake(tkmenu, panel$.handle)
  fileMenu <- list(.type="menu")
  handshake(tkconfigure, panel$.handle, menu=menu$.widget)
  w.appearancewidget(menu, font, foreground, background)
  for (i in (1:length(labels)))
  {
    submenu <- unlist(labels[i])
    fileMenu$.widget <- handshake(tkmenu, menu$.widget,tearoff=FALSE)
    for (j in (2:length(submenu)))
    {
      handshake(tkadd, fileMenu$.widget, "command", label=submenu[j], command=eval(parse(text=paste("f('", submenu[j], "')", sep="")))) # The eval is necessary.
    }
    handshake(tkadd, menu$.widget,"cascade",label=submenu[1],menu=fileMenu$.widget)
    w.appearancewidget(fileMenu, font, foreground, background)
  }
  
  invisible(menu)
}

rp.menu <- function(panel, variable, labels, initval=NULL, action=I, foreground=NULL, background=NULL, font=NULL, 
  name=paste("menu", .nc(), sep="")) 
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
  # variable then isn't used - but could be for a default initial value

  f <- function(val)
  {
    rp.var.put(panelname, varname, val)
    panel <- rp.control.get(panelname)
    panel <- action(panel)
    rp.control.put(panelname, panel)
  }
    
  widget <- w.menu(panel, labels, action=f, foreground, background, font)
  rp.widget.put(panelname, name, widget)
  if (.rpenv$savepanel) { rp.control.put(panelname, panel) } # put the panel back into the environment
  invisible(panelname)
}
